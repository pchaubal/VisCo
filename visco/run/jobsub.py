# create a script to submit jobs
# create a method to postproces the jobs
# 
import os

snapshot_dir = '../snapshots'

def copycode(dir_name):
	print '\ncopying code...'
	os.system('mkdir ' + dir_name)
	copycommand = 'cp ../readsnap.f90 ../Makefile ../fields.f90 ../visco.f90 ' + dir_name
	os.system(copycommand)
	print 'Made folder',dir_name, 'and code copied.'

def edit_runvisco(snapname):
	# Read in the file
	with open('runvisco_template', 'r') as file :
		filedata = file.read()

	# Replace the target string
	# print filedata
	filedata = filedata.replace('snap001', snapname)

	# Write the file out again
	jobscript = 'runvisco_' + snapname
	with open(jobscript, 'w') as file:
		file.write(filedata)


# edit_runvisco('snap001')



def edit_visco(snapname):
	viscofile = '../parallelcode/' + snapname + '/visco.f90'
	with open(viscofile, 'r') as file :
		filedata =file.read()

	datafilepath = '../../data/' + snapname
	filedata = filedata.replace('./data',datafilepath)

	with open(viscofile,'w') as file:
		file.write(filedata)

def edit_readnap(snapname):
	readsnap = '../parallelcode/' + snapname + '/readsnap.f90'
	with open(readsnap, 'r') as file :
		filedata =file.read()

	datafilepath = '../snapshots/' + snapname
	filedata = filedata.replace('snapshot_021',datafilepath)

	with open(readsnap,'w') as file:
		file.write(filedata)


def submitjob(snapname):
	jobsub = 'bsub < '+'runvisco_' + snapname
	os.system(jobsub)
	print 'job summitted'

for file in os.listdir(snapshot_dir):
	codedir = '../parallelcode/' + file
	copycode(codedir)
	edit_visco(file)
	edit_readnap(file)
	edit_runvisco(file)
	submitjob(file)

# clean up the run folder after submitting all the jobs
os.system('rm runvisco_snap*')