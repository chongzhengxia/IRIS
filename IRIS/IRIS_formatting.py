import sys, numpy, argparse, os
#                                          []                  []                 1文件或者2目录    {}
# loadSamplelist(args.rmats_sample_order, sample_fin_list,  sample_header,    sample_name_field, sample_size)
def loadSamplelist(fin_samples,           sample_fin_list,  sample_header,    sample_name_field, sample_size):
	for l in open(fin_samples):   # 
		ls=l.strip()
		sample_fin_list.append(ls)
		for r in open(ls):
			#     既可以从bam文件解析文件名，也可以从bam文件的目录解析文件名 sample_name_field==2
			rs=map(lambda x:x.split('/')[-sample_name_field].split('.bam')[0].split('.aln')[0],r.strip().strip(',').split(','))
			# rs sample1
			if sample_name_field==2:
				sn_list=r.strip().strip(',').split(',') # sn_list [/path/sample1.aln/Aligned.sortedByCoord.out.bam]
				for e,sn in enumerate(rs):
					if len(sn)==0:
						rs[e]=sn_list[e].split('/')[-1].split('.')[0]
			sample_header+=rs
			sample_size[ls]=len(r.split(','))
	return sample_fin_list, sample_header, sample_size

def parseEventRowSE(line_split):
	return line_split[1].strip('"')+'\t'+line_split[2].strip('"')+'\t'+'\t'.join(line_split[3:7]+line_split[8:10])

def parseEventRow(line_split):
	return line_split[1].strip('"')+'\t'+line_split[2].strip('"')+'\t'+'\t'.join(line_split[3:11])

def loadGTF(gtf):
	exon_start_dict={}
	exon_end_dict={}
	for l in open(gtf):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if ls[2]=='exon':
			chrom=ls[0]
			if chrom.startswith('chr')==False:
				chrom='chr'+chrom
			exon_start_dict[ls[6]+':'+chrom+':'+ls[3]]=''
			exon_end_dict[ls[6]+':'+chrom+':'+ls[4]]=''
	return exon_start_dict, exon_end_dict

def checkNovelSS(head, ls, splicing_event_type, exon_start_dict, exon_end_dict):# This is a conservative def of novelSS than rMATS4.1 (0.4% events less- complex cases )
	ld=dict(zip(head, ls))
	strand=ld['strand']
	chrom=ld['chr']
	if splicing_event_type=='SE':
		check1=strand+':'+chrom+':'+str(int(ld['exonStart_0base'])+1) not in exon_start_dict
		check2=strand+':'+chrom+':'+str(int(ld['downstreamES'])+1) not in exon_start_dict
		check3=strand+':'+chrom+':'+ld['exonEnd'] not in exon_end_dict	
		check4=strand+':'+chrom+':'+ld['upstreamEE'] not in exon_end_dict
	elif splicing_event_type=='A5SS':
		check1=strand+':'+chrom+':'+str(int(ld['shortES'])+1) not in exon_start_dict
		check2=strand+':'+chrom+':'+str(int(ld['flankingES'])+1) not in exon_start_dict
		check3=strand+':'+chrom+':'+ld['longExonEnd'] not in exon_end_dict	
		check4=False
		#check4=strand+':'+chrom+':'+str(int(ld['shortES'])+1) not in exon_end_dict
	elif splicing_event_type=='A3SS':
		check1=strand+':'+chrom+':'+str(int(ld['longExonStart_0base'])+1) not in exon_start_dict
		check2=strand+':'+chrom+':'+str(int(ld['shortES'])+1) not in exon_start_dict
		check3=strand+':'+chrom+':'+ld['flankingEE'] not in exon_end_dict	
		check4=False

	elif splicing_event_type=='RI':
		check1,check2,check3,check4=[False,False,False,False]
	else:
		exit('choose AS type.')
	return check1, check2, check3, check4


def mergeEvents(events_fin_list, splicing_event_type, novelSS, exon_start_dict, exon_end_dict): 
	# 遍历每个文件，然后将unique事件全部存储到字典当中并返回
	parseRow=parseEventRow
	if splicing_event_type=='SE':
		parseRow=parseEventRowSE
	total_event_dict={}
	for i, events_fin in enumerate(events_fin_list):
		# 遍历文件的每一行，将ENSG编号，基因名称，染色体号，正负链，前中后剪接信息返回为制表符分隔的文本信息，如果这些信息在所有的
		# 事件字典中存在，则继续遍历
		head=[]
		for index,event_l in enumerate(open(events_fin)):
			if index==0:
				head=event_l.strip().split('\t')
				continue
			event_ls=event_l.strip().split('\t')
			if novelSS:
				check1, check2, check3, check4= checkNovelSS(head, event_ls, splicing_event_type, exon_start_dict, exon_end_dict)
				novel=True if check1 or check2 or check3 or check4 else False
				if novel==False: # if no novel, will not parse the row and save
					continue
			events_cord=parseRow(event_ls)
			# 'ENSG00000156261\tCCT8\tchr21\t-\t29066715\t29067071\t29066715\t29066777\t29066890\t29067071'
			if events_cord in total_event_dict:
				continue
			total_event_dict[events_cord]=''
	return total_event_dict

def writeMergedEvents(events_fin_list, splicing_event_type, cov_cutoff, data_name, fout_path, novelSS, exon_start_dict, exon_end_dict):
	total_event_dict=mergeEvents(events_fin_list, splicing_event_type, novelSS, exon_start_dict, exon_end_dict)
	novelss_tag=''
	if novelSS:
		novelss_tag='.novelSS'
	total_event_list=sorted(total_event_dict.keys())
	fout=open(fout_path+'/prefilter_events.splicing_matrix.'+splicing_event_type+novelss_tag+'.cov'+str(cov_cutoff)+'.'+data_name+'.txt','w')
	for e in total_event_list:
		fout.write(e.strip()+'\n')
	fout.close()

	return total_event_list
    # mergeMatrixInBatch(fin_list, events_fin_list, sample_fin_list, cov_cutoff, data_name, splicing_event_type, sample_header, sample_size, total_event_list, file_batch_list, batch, fout_path, individual_filter, novelSS)
def mergeMatrixInBatch(fin_list, events_fin_list, sample_fin_list, cov_cutoff, data_name, splicing_event_type, sample_header, sample_size, total_event_list, file_batch_list, batch, fout_path, individual_filter, novelSS):
	parseRow=parseEventRow
	if splicing_event_type=='SE':
		parseRow=parseEventRowSE
		
	for b in range(0,len(total_event_list),batch):
		Intercep_Matrix={}      #{'事件'： {来源:(psi, cov)}}
		print '[INFO] Merging in progress. Working on batch ',b 
		batch_event_list= total_event_list[b:min(b+batch,len(total_event_list))]
		batch_event_dict= dict.fromkeys(batch_event_list, 0)
		# 将batch_event_dict 当作完全从 total_env_list 里初始化得来的字典，简化逻辑
		
		for n,fin in enumerate(fin_list):  # 通过n将每个样本的JC.raw.input.A3SS.txt与fromGTF.A3SS.txt一一对应起来

			eventID={}
			for index,event_l in enumerate(open(events_fin_list[n])):
				if index==0:
					continue
				event_ls=event_l.strip().split('\t')    
				event_cord=parseRow(event_ls)
	                       #事件字符串                # {事件：0}
				if event_cord in batch_event_dict:  
					eventID[event_ls[0]]=event_cord # {0：事件， 1：事件， 2：事件}如果从文件中获得的事件在切片排序后的字典中，则在eventID字典中重构这个
					                                # 事件，以该事件在该文件的ID号为key,以事件的信息为value
			print '[INFO] Merging file: ', fin, len(eventID)
			for index,r in enumerate(open(fin)):
				if index==0:
					continue
				rs=r.strip().split('\t')
				if rs[0] not in eventID:               # 如果JC.txt与fromGTF.txt的ID能够对应上，则计算psi和cov
					continue
				Incl=map(float,rs[1].split(','))
				Skip=map(float,rs[2].split(','))
				Cov=[num+Skip[o] for o,num in enumerate(Incl)]
				psi_values=[]
				for i,I in enumerate(Incl):
					if individual_filter: # individual_filter. Use Cov[i] for each individual sample
						if  Cov[i]< cov_cutoff:
							psi_values.append('NaN')
						else:
							psi_values.append(str(round(I/int(rs[5])/(I/int(rs[5])+Skip[i]/int(rs[6])),4)))
					else:
						if int(I)+int(Skip[i])==0:
							psi_values.append('NaN')
						else:
							psi_values.append(str(round(I/int(rs[5])/(I/int(rs[5])+Skip[i]/int(rs[6])),4)))

				if eventID[rs[0]] not in Intercep_Matrix:
					Intercep_Matrix[eventID[rs[0]]]={}
				if sample_fin_list[n] not in Intercep_Matrix[eventID[rs[0]]]:
					Intercep_Matrix[eventID[rs[0]]][sample_fin_list[n]]=(psi_values,Cov)   {'事件'： {来源:(psi, cov), 来源:(psi, cov)}}
				if len(psi_values)!=sample_size[sample_fin_list[n]]:
					exit('[Abort] Sample number does not match observations in JC file.')

		novelss_tag=''
		if novelSS:
			novelss_tag='.novelSS'
		file_path_name=fout_path+'/splicing_matrix/splicing_matrix.'+splicing_event_type+novelss_tag+'.cov'+str(cov_cutoff)+'.'+data_name+'.txt.batch_'+str(b)+'.txt'
		file_batch_list.append(file_path_name)
		fout=open(file_path_name,'w')
		header_line='AC\tGeneName\tchr\tstrand\texonStart\texonEnd\tupstreamEE\tdownstreamES\t'+'\t'.join(sample_header)
		if splicing_event_type=='A5SS':
			header_line='AC\tGeneName\tchr\tstrand\tlongExonStart\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\t'+'\t'.join(sample_header)
		if splicing_event_type=='A3SS':
			header_line='AC\tGeneName\tchr\tstrand\tlongExonStart\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\t'+'\t'.join(sample_header)
		if splicing_event_type=='RI':	
			header_line='AC\tGeneName\tchr\tstrand\triExonStart\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\t'+'\t'.join(sample_header)
		fout.write(header_line+'\n')
		for k in sorted(Intercep_Matrix.keys()):
			psi_value_all=[]
			cov_all=[]
			for sample in sample_fin_list:
				if sample in Intercep_Matrix[k]:
					psi_value_all+=Intercep_Matrix[k][sample][0]
					cov_all+=Intercep_Matrix[k][sample][1]
				else:
					psi_value_all+=['NaN']*sample_size[sample]
			#if individual_filter==False: #if filter by group and cov < cov_cutoff, skip this event K. Otherwise, psi_vallue_all is being wrote to output.
			#	mean=numpy.mean(cov_all)
			#	if mean < cov_cutoff: 
			#		continue
			if set(psi_value_all)==set(['NaN']): #remove full NaN events 2020
				continue
			fout.write(k+'\t'+'\t'.join(psi_value_all)+'\n')    # example 的txt文件

		fout.close()
	return file_batch_list   # 列表里只有一个文件地址，简化

def mergeMatrixInOne(file_batch_list, cov_cutoff, data_name, splicing_event_type, fout_path, novelSS):
	novelss_tag=''
	if novelSS:
		novelss_tag='.novelSS'
	fout_merge=open(fout_path+'/splicing_matrix/splicing_matrix.'+splicing_event_type+novelss_tag+'.cov'+str(cov_cutoff)+'.'+data_name+'.txt','w')
	header=0
	for file_batch in file_batch_list:
		for j,l in enumerate(open(file_batch)):
			if j==0:
				if header==0:
					header+=1
					fout_merge.write(l)
				continue
			fout_merge.write(l)
	fout_merge.close()
	os.system('rm '+fout_path+'/splicing_matrix/splicing_matrix.'+splicing_event_type+novelss_tag+'.cov'+str(cov_cutoff)+'.'+data_name+'.txt.batch_*.txt')	
	return 'splicing_matrix.'+splicing_event_type+novelss_tag+'.cov'+str(cov_cutoff)+'.'+data_name+'.txt'
	
# index_PsiMatrix(fout_path+'/splicing_matrix/'+merged_file_name,   fout_path+'/splicing_matrix',            '\t',               splicing_event_type)
def index_PsiMatrix(fn,                                              outdir,                                 delim,              splicing_event_type):
	out_fp = outdir+'/'+fn.split('/')[-1]+'.idx'
	line_formatter =  "{id}\t{offset}\n"
	offset = 0	
	col_index=10
	if splicing_event_type=='SE':#handle SE and other types of AS events
		col_index=8
	with open(fn, 'r') as fin:
		with open(out_fp, 'w') as fout:
			offset += len(fin.readline())    # 自动跳过了头行
			for line in fin:
				ele = line.strip().split(delim)
				eid = ':'.join([ele[0].split('_')[0].split('.')[0]]+ele[1:col_index])
				fout.write( line_formatter.format(id=eid, offset=offset) )
				offset += len(line)
	return

def main(args):
	cov_cutoff=args.cov_cutoff
	data_name=args.data_name
	sample_name_field=args.sample_name_field            # 1代表是文件的绝对路径，2代表的是文件的目录的绝对路径
	splicing_event_type=args.splicing_event_type        # 不包含互斥外显子
	individual_filter= args.sample_based_filter
	novelSS= args.novelSS
	exon_start_dict={}
	exon_end_dict={}
	if novelSS:
		gtf=args.gtf
		exon_start_dict, exon_end_dict= loadGTF(gtf)
	
	if sample_name_field==1:
		print '[INFO] Sample name parsed from bam file.  (alternatively can be parsed from up level folder)'
	if sample_name_field==2:
		print '[INFO] Sample name parsed from folder name above the bam file.  (alternatively can be parsed from bam file)'
	db_dir=args.iris_db_path.rstrip('/')             # 储存经过格式化以及索引的ASmatrix
	#prepare files/folders in IRIS db directory
	os.system('mkdir -p '+db_dir+'/'+data_name+' '+db_dir+'/'+data_name+'/splicing_matrix')
	fout_path=db_dir+'/'+data_name
	print '[INFO] output path: '+fout_path
	fin_list=[]
	sample_fin_list=[]
	events_fin_list=[]
	sample_size={}
	sample_header=[]
	file_batch_list=[]


	#PARSING INPUT FILE LISTS
	# rmats_mat_path_manifest参数可以输入多个rMATS输出文件的目录，fin_list列表包含了所有的JC.raw.input文件的绝对路径
	fin_list=[l.strip().rstrip('/')+'/JC.raw.input.'+splicing_event_type+'.txt' for l in open(args.rmats_mat_path_manifest)]
	# 列表包含了fromGTF文件的绝对路径 可以包含多个样本的相同事件
	events_fin_list=[l.strip().rstrip('/')+'/fromGTF.'+splicing_event_type+'.txt' for l in open(args.rmats_mat_path_manifest)]
	sample_fin_list, sample_header, sample_size= loadSamplelist(args.rmats_sample_order,sample_fin_list, sample_header,sample_name_field, sample_size)

	#MAKING MERGED EVENTS LIST
	total_event_list= writeMergedEvents(events_fin_list, splicing_event_type, cov_cutoff, data_name, fout_path, novelSS, exon_start_dict, exon_end_dict)
	if args.merge_events_only:
		exit('[INFO] Done merging events only.')
	print '[INFO] Done loading file dir. Total events:', len(total_event_list)

	#START MERGING MATRICES IN BATCH MODE FOLLOWING EVENTS LIST GENERATED.
	batch=20000
	file_batch_list=mergeMatrixInBatch(fin_list, events_fin_list, sample_fin_list, cov_cutoff, data_name, splicing_event_type, sample_header, sample_size, total_event_list, file_batch_list, batch, fout_path, individual_filter, novelSS)
	print '[INFO] Done merging matrices by batch.'
	merged_file_name=mergeMatrixInOne(file_batch_list, cov_cutoff, data_name, splicing_event_type, fout_path, novelSS)
	print '[INFO] Done merging matrices: '+merged_file_name

	
	#create index in IRIS db directory
	index_PsiMatrix(fout_path+'/splicing_matrix/'+merged_file_name,fout_path+'/splicing_matrix','\t', splicing_event_type)
	print '[INFO] Finished. Created matrix: '+fout_path

if __name__ == '__main__':
	main()
