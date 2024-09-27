import sys, csv, glob, os
from . import config
                         # writeShell(rMATS_path, fin_name, folder_name,          bam_dir, str(read_length), gtf, novelSS, task_name, task_dir)
            # rMATS程序路径       bam文件绝对路径  bam文件dirname    一级目录   
def writeShell(rMATS_path,              fin_name,    folder_name,     bam_dir, read_length_argument, gtf, novelSS, task_name, task_dir):
	fout_local=open(folder_name+'/bam_list.txt','w')  # 在bam文件的同一级目录中打开个文件
	fout_local.write(fin_name)                        # 在bam文件的同一级目录写入bam文件的绝对路径到一个文件里
	fout_local.close()

	sample_name=folder_name.split('/')[-1].split('.')[0]
	task_script_base = 'rMATS_prep.{}.sh'.format(sample_name)
	task_script = os.path.join(task_dir, task_script_base)
	fout=open(task_script,'w')
	fout.write('#!/bin/bash\n')
	novelSS_str=''
	if novelSS:
		novelSS_str='--novelSS '
	# TODO the '|| true' at the end of this command ignores a
	# failure return code from the python command.
	# rMATS produces the desired output file despite the error return.
	# A future version of rMATS may fix this behavior.
	fout.write('python {} --b1 {}/bam_list.txt --od {} --tmp {}/{}.RL{}/{}.tmp --anchorLength 1 --readLength {} --gtf {} -t paired --task prep --nthread 8 --statoff {}|| true\n'.format(rMATS_path, folder_name, folder_name, bam_dir, task_name, read_length_argument, sample_name, read_length_argument, gtf, novelSS_str))
	fout.close()
        # 每一个bam文件运行prep关闭statoff统计信息
                             #      找到的所有的star日志
def organizeReadLength(rMATS_path, file_list_mapping, gtf, novelSS, bam_prefix, task_name, task_dir):
	rl_dict={}           {folder_name:RL}
	folder_names={}      {RL:''}
	for fin_name in file_list_mapping:
		for l in open(fin_name):
			if l.find('Average input read length |')!=-1:
				map_rl=int(round(float(l.split('Average input read length |')[-1].strip())/2,0))  # 通过日志找到的RL
				rl_dict['/'.join(fin_name.split('/')[:-1])]=map_rl   {folder_name:RL}
				folder_names[map_rl]=''                               {RL:''}
				break
	bam_dir='/'.join(file_list_mapping[0].split('/')[:-2]) # bam_dir
	for folder_name in folder_names:
		os.system('mkdir -p '+bam_dir+'/'+task_name+'.RL'+str(folder_name))
	for folder_name in rl_dict:
		writeShell(rMATS_path, folder_name+'/'+bam_prefix+'.bam', folder_name, bam_dir, str(rl_dict[folder_name]), gtf, novelSS, task_name, task_dir)

def main(args):
	gtf=args.gtf
	task_name=args.data_name
	rMATS_path=args.rMATS_path         # rMATS的代码的路径
	bam_dir= args.bam_dir.rstrip('/')  # process_rna 的一级输出目录，二级目录是star的输出目录，是process_rna的 -p参数 --sampleID-outdir
	bam_prefix=args.bam_prefix         # bam文件的前缀，一般是Aligned.sortedByCoord.out
	novelSS=args.novelSS               # 是否开启rMATS 的novelSS选项，一般是关闭 
	task_dir=args.task_dir             # 任务脚本目录的路径
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)
	if args.read_length:              # read_length的设置，一般是从star的输出结果里获取
		read_length= int(args.read_length)
	
	print 'preparing rMATS-turbo prep directories'
	if args.read_length==False:
		mapping_log_file_list=glob.glob(bam_dir+'/*/Log.final.out')   # 找到star的输出结果从而获取read_length，应该是有个bug 应该加上Log.final.out的前缀
		organizeReadLength(rMATS_path, mapping_log_file_list, gtf, novelSS, bam_prefix, task_name, task_dir) #relocated based on the read length
	else:
		mapping_bam_list=glob.glob(bam_dir+'/*/'+bam_prefix+'.bam')  # 包含bam文件路径的列表
		os.system('mkdir -p '+bam_dir+'/'+task_name+'.RL'+str(read_length))
		for fin_name in mapping_bam_list:   #mapping_bam_list是列表，包含了全部sort后的bam文件的路径
			folder_name= '/'.join(fin_name.split('/')[:-1]) #bam文件的dirname
			# a/b/filename fin_name: a/b/filename folder_name: a/b
			
			writeShell(rMATS_path, fin_name, folder_name,          bam_dir, str(read_length), gtf, novelSS, task_name, task_dir)
                         # rMATS程序路径       bam文件绝对路径  bam文件dirname    一级目录    
 
if __name__ == '__main__':
	main()
