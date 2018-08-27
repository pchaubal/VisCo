import numpy as np
import os
import itertools
import argparse





parser = argparse.ArgumentParser()
parser.add_argument("datadir", help="give the path to the data directory from the current directory")
datadir = parser.parse_args().datadir

fields_path = datadir+'/total.dat'
positionfile = datadir+'/positions.dat'
binsize = 0.5 # In the units of Mpc

print 'loading data from ' + datadir
data =[]
data = np.loadtxt(fields_path)
pos = np.loadtxt(positionfile)

field=[]
for i in range(len(data)):
	if not np.isnan(data[i]).any():
		elements = [pos[i]] + [a for a in data[i]]
		field.append(elements)

def combinations(source):
	print 'Calculating correlations ...'
	result = []
	for p1 in range(len(source)):
	        for p2 in range(p1+1,len(source)):
	        	r       = np.linalg.norm(source[p1][0]-source[p2][0])
	        	P_AD    = source[p1][1]*source[p2][3]
	        	P_AT    = source[p1][3]*source[p2][2]
	        	P_DD    = source[p1][1]*source[p2][1]
	        	P_DT    = source[p1][1]*source[p2][2]
	        	P_TT    = source[p1][2]*source[p2][2]
	        	d2_P_DD = source[p1][1]*source[p2][4]
	        	d2_P_DT = source[p1][4]*source[p2][2]
	        	d2_P_TT = source[p1][2]*source[p2][5]
	        	row     = [r,P_AD,P_AT,P_DD,P_DT,P_TT,d2_P_DD,d2_P_DT,d2_P_TT]
	        	result.append(row)
	return result

pairs = combinations(field)
np.savetxt(datadir+'/pairs.dat',pairs,fmt='%1.4e', delimiter='\t')

# a=[[1,2,1],[2,1,2],[3,2,5],[4,5,3]]
# bins = np.linspace(0,100,100)
# hist = np.histogram(pairs,bins=bins, normed=False)

# Now create the histograms to average the computed powers in a small range i.e. binwidth
BLOCK = 65536 # iteratting over block size to minimize memory
max_length = 100.0
hist=[]
for bin_no in np.arange(max_length/binsize):
	P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0
	N=0
	for i in np.arange(0,len(pairs),BLOCK):
		tmp_pairs= pairs[i:i+BLOCK]
		for j in np.arange(len(tmp_pairs)):
			keep  =  tmp_pairs[j][0] >= bin_no*binsize
			keep &=  tmp_pairs[j][0] <  (bin_no+1)*binsize
			if keep:
				P_AD    += tmp_pairs[j][1] 
				P_AT    += tmp_pairs[j][2] 
				P_DD    += tmp_pairs[j][3] 
				P_DT    += tmp_pairs[j][4] 
				P_TT    += tmp_pairs[j][5] 
				d2_P_DD += tmp_pairs[j][6] 
				d2_P_DT += tmp_pairs[j][7] 
				d2_P_TT += tmp_pairs[j][8]
				N       += 1
	if N!=0:
		# If bin is not empty, then normalize
		P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = P_AD/float(N), P_AT/float(N), P_DD/float(N), P_DT/float(N), P_TT/float(N), d2_P_DD/float(N), d2_P_DT/float(N), d2_P_TT/float(N) 
	else:
		# Otherwise set all correlation functions to zero
		P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0																																												
	normalized_corr = [(bin_no*binsize + binsize/2.0),P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT]
	hist.append(normalized_corr)
	print bin_no*binsize,(bin_no+1)*binsize,N

np.savetxt(datadir+'/correlation.dat',hist,fmt='%1.4e', delimiter='\t')
print 'Calculating sound velocities...'

rlist = []
cs2list     = []
cs2numlist  = []
cs2denlist  = []
cv2list     = []
cv2numlist  = []
cv2denlist  = []
for i in range(len(hist)):
	try:
		cs2 = (hist[i][2]*hist[i][7] - hist[i][1]*hist[i][8])/( (hist[i][7])**2 - hist[i][6]*hist[i][8] )
		cv2 = (hist[i][1]*hist[i][7] - hist[i][2]*hist[i][6])/( (hist[i][7])**2 - hist[i][6]*hist[i][8] )
	except ZeroDivisionError:
		# print "cs2 is infinity"
		cs2 = 0
		cv2 = 0

	cs2num = hist[i][2]*hist[i][7] - hist[i][1]*hist[i][8]
	cs2den = ( (hist[i][7])**2 - hist[i][6]*hist[i][8] )
	cv2num = (hist[i][1]*hist[i][7] - hist[i][2]*hist[i][6])
	cv2den = ( (hist[i][7])**2 - hist[i][6]*hist[i][8] )


	# print "cs2,cv2=",i,cs2,cv2
	rlist.append(hist[i][0])
	cs2list.append(cs2)
	cs2numlist.append(cs2num)
	cs2denlist.append(cs2den)
	cv2list.append(cv2)
	cv2numlist.append(cv2num)
	cv2denlist.append(cv2den)

np.savetxt(datadir+'/soundvels.dat',hist,fmt='%1.4e', delimiter='\t')
np.savetxt(datadir+'/rlist.dat',rlist,fmt='%1.4e', delimiter='\t')
np.savetxt(datadir+'/cs2.dat',cs2list,fmt='%1.4e', delimiter='\t')
np.savetxt(datadir+'/cv2.dat',cv2list,fmt='%1.4e', delimiter='\t')
