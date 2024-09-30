import sys,glob,os,argparse
from . import config
                          #      task_name.RL150
def write_task_script(rMATS_path, bam_dir, task_name, gtf, novelSS, task_dir):
	#                task_name.RL150
	read_length=int(bam_dir.split('/')[-1].split('.')[-1][2:])   # 150
	dir_name=task_name+'_RL'+str(read_length)                    # task_name_RL150
	graphlist=glob.glob(bam_dir+'/*.tmp/*')
	os.system('mkdir -p '+bam_dir+'/'+dir_name+'.graph')
	os.system('cp '+bam_dir+'/*.tmp/* '+bam_dir+'/'+dir_name+'.graph/.')   # 将task_name.RL150/.tmp文件夹下的所有文件移动到 task_name.RL150/task_name_RL150.graph下
	print '[INFO] Done copy'
	                        task_name_RL150
	cmd='head -n1 -q '+bam_dir+'/'+dir_name+'.graph/*.rmats |paste -d, -s >'+bam_dir+'/'+dir_name+'_rmatspost_list.txt'
	print cmd
	os.system(cmd)

	task_script_base = 'rMATS_post.{}.sh'.format(dir_name)
	task_script = os.path.join(task_dir, task_script_base)
	fout=open(task_script,'w')

	novelSS_str=''
	if novelSS:
		novelSS_str='--novelSS '
	# TODO the '|| true' at the end of this command ignores a
	# failure return code from the python command.
	# rMATS produces the desired output files despite the error return.
	# A future version of rMATS may fix this behavior.
									# task_name_RL150_rmatspost_list.txt
	fout.write('#!/bin/bash\npython '+rMATS_path+' --b1 '+bam_dir+'/'+dir_name+'_rmatspost_list.txt --od '+bam_dir+'/'+dir_name+'.matrix --tmp '+bam_dir+'/'+dir_name+'.graph/ --anchorLength 1 --readLength '+str(read_length)+' --gtf '+gtf+' -t paired --nthread 8 --task post --statoff '+novelSS_str+'|| true\n')
		                                                                            # task_name.RL150              task_name_RL150.matrix
	fout.close()
	return 


def main(args):
	rMATS_path=args.rMATS_path
	prep_dir=args.bam_dir.rstrip('/')  # 一级目录，二级目录包含bam文件，
	gtf=args.gtf
	task_name=args.data_name
	novelSS=args.novelSS
	task_dir=args.task_dir
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)

	rl_bam_folders=glob.glob(prep_dir+'/'+task_name+'.RL*')  # 找到全部的rmats prep步骤的输出tmp目录的高一级目录
	#   task_name.RL150
	for bam_folder in rl_bam_folders:
		print bam_folder
		write_task_script(rMATS_path, bam_folder,task_name, gtf, novelSS, task_dir)


if __name__ == '__main__':
	main()
