# create a script to submit jobs
# create a method to postproces the jobs
# 
import os
import numpy as np


def copycode(dir_name):
	print '\ncopying code...'
	os.system('mkdir ' + dir_name)
	copycommand = 'cp ../readsnap.f90 ../Makefile ../fields.f90 ../visco.f90 ' + dir_name
	os.system(copycommand)
	print 'Made folder',dir_name, 'and code copied.'

def edit_runvisco(Om):
	# Read in the file
	with open('runvisco_template', 'r') as file :
		filedata = file.read()

	# Replace the target string
	# print filedata
	filedata = filedata.replace('snap001', 'varying_Om/Om='+str(Om))
	filedata = filedata.replace('snap', 'varying_Om/Om='+str(Om)[0:6])


	# Write the file out again
	jobscript = 'runvisco_' + 'Om='+str(Om)
	with open(jobscript, 'w') as file:
		file.write(filedata)


# edit_runvisco('snap001')
def makedatadir(Om):
	os.system('mkdir ../data/varying_Om/Om='+ str(Om)[0:6])
	


def edit_visco(Om):
	viscofile = '../parallelcode/varying_Om/Om=' + str(Om) + '/visco.f90'
	with open(viscofile, 'r') as file :
		filedata =file.read()

	datafilepath = '../../../data/varying_Om/Om=' + str(Om)[0:6]
	filedata = filedata.replace('./data',datafilepath)

	with open(viscofile,'w') as file:
		file.write(filedata)

def edit_readnap(Om):
	readsnap = '../parallelcode/varying_Om/Om=' + str(Om) + '/readsnap.f90'
	with open(readsnap, 'r') as file :
		filedata =file.read()

	datafilepath = '../../../../gadget2/lcdm/snapshots/Om_dm=' + str(Om) + '/snapshot_041'
	filedata = filedata.replace('snapshot_021',datafilepath)

	with open(readsnap,'w') as file:
		file.write(filedata)


def submitjob(Om):
	jobsub = 'bsub < '+'runvisco_Om=' + str(Om)
	os.system(jobsub)
	print 'job summitted'

for Om in np.linspace(0.2,0.3,100):
	codedir = '../parallelcode/varying_Om/Om=' + str(Om)
	copycode(codedir)
	makedatadir(Om)
	edit_visco(Om)
	edit_readnap(Om)
	edit_runvisco(Om)
	submitjob(Om)

# clean up the run folder after submitting all the jobs
os.system('rm runvisco_Om*')