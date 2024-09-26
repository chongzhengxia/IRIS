import sys,glob,os
from . import config
        # out_dir/sample1   sample1全路径    这个函数在for 里被调用多次
def write_task_script(out_prefix, fin, label_string, starGenomeDir, gtf, task_dir):
	fq_dir=fin   #sample1的全路径
	fqs=glob.glob(fq_dir+'/*')    # [sample1R1的全路径，sample1R2的全路径]
	r1=[] #R1文件    一个样本的测序结果可能有多个lane以及文库，产生了多个R1文件
	r2=[] #R2文件
	for fq in fqs:
		if fq.find('1'+label_string+'f')!=-1:
			r1.append(os.path.abspath(fq))  #找到的R1文件全部放到r1列表里
		elif fq.find('2'+label_string+'f')!=-1:
			r2.append(os.path.abspath(fq))  #找到的R2文件全部放到r2列表里
	if len(r1)!=len(r2) or len(r1)==0 or len(r2)==0:
		print 'file name can not be recognize'
		return '','' 

	out_dir=out_prefix.rstrip('/')+'.aln'                            # out_dir/sample1.aln
	sample_name=out_dir.split('/')[-1].split('.')[0]                 # sample1
	task_script_base1 = 'STARmap.{}.sh'.format(sample_name)          # STARmap.sample1.sh
	task_script1 = os.path.join(task_dir, task_script_base1)         # task_dir/STARmap.sample1.sh  路径名字
	task_script_base2 = 'Cuffquant.{}.sh'.format(sample_name)        # Cuffquant.sample1.sh
	task_script2 = os.path.join(task_dir, task_script_base2)         # task_dir/Cuffquant.sample1.sh  路径名字
	fout1=open(task_script1,'w')                                     # 在task目录下创建了 STARmap.sample1.sh 文件
	fout2=open(task_script2,'w')
	fq_path=','.join(sorted(r1)+sorted(r2))                          # r1,r2排序，确保每个R1,R2一一对应，最终整合到一个列表里，列表里先是全部的R1,然后是全部的R2

	cmd1='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --mapping --sort -p '+out_dir+' '+fq_path  # 比对以及排序
	fout1.write('#!/bin/bash\n'+cmd1+'\n')
	cmd2='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --quant -p '+out_dir+' '+fq_path           # 定量以及检测AS事件
	fout2.write('#!/bin/bash\n'+cmd2+'\n')
	return task_script1,task_script2                  # 为什么需要返回值

def main(args):
	starGenomeDir=args.starGenomeDir
	gtf=args.gtf
	fastq_folder_dir=args.fastq_folder_dir.rstrip('/') #一级目录
	fastq_folder_list=glob.glob(fastq_folder_dir+'/*') #找出一级目录下的全部二级目录，每个二级目录里包含每个样本的R1，R2文件
	out_dir=args.outdir.rstrip('/')                    #bam文件的文件夹的输出目录
	task_name=args.data_name
	label_string=args.label_string
	os.system('mkdir -p '+out_dir)
	task_dir=args.task_dir				   #存放作业脚本的目录
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)
                       # 列表包含全部的二级目录
	for folder in fastq_folder_list:
		print folder                  # out_dir/sample1          sample1的全路径                                存放作业脚本的目录
		fn1,fn2=write_task_script(out_dir+'/'+folder.split('/')[-1], folder, label_string, starGenomeDir, gtf, task_dir)
		if fn1=='':
			continue


if __name__ == '__main__':
	main()
