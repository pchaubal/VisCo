import numpy as np
import os
import itertools

fields_path = './data/total.dat'
positionfile = './data/positions.dat'
binsize = 0.1 # In the units of Mpc

print 'loading data ...'
data =[]
data = np.loadtxt(fields_path)
pos = np.loadtxt(positionfile)

field=[]
for i in range(len(data)):
	if not np.isnan(data[i]).any():
		elements = [pos[i]] + [a for a in data[i]]
		field.append(elements)

def combs(source):
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

pairs = combs(field)
np.savetxt('./postprocess/pairs.dat',pairs,fmt='%1.4e', delimiter='\t')


# print len(pairs), len(pairs[0])

N_bins = int(104/float(binsize))
hist = []

for bin_no in range(N_bins):
	P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0
	N=0
	for i in range(len(pairs)):
		if (bin_no*binsize < pairs[i][0] < (bin_no+1)*binsize):
			P_AD    = P_AD    + pairs[i][1] 
			P_AT    = P_AT    + pairs[i][2] 
			P_DD    = P_DD    + pairs[i][3] 
			P_DT    = P_DT    + pairs[i][4] 
			P_TT    = P_TT    + pairs[i][5] 
			d2_P_DD = d2_P_DD + pairs[i][6] 
			d2_P_DT = d2_P_DT + pairs[i][7] 
			d2_P_TT = d2_P_TT + pairs[i][8]
			N=N+1

	if N!=0:
		# If bin is not empty, then normalize
		P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = P_AD/float(N), P_AT/float(N), P_DD/float(N), P_DT/float(N), P_TT/float(N), d2_P_DD/float(N), d2_P_DT/float(N), d2_P_TT/float(N) 
	else:
		# Otherwise set all correlation functions to zero
		P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0																																												
	normalized_corr = [(bin_no*binsize + binsize/2.0),P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT]
	hist.append(normalized_corr)
	print bin_no*binsize,(bin_no+1)*binsize,N

np.savetxt('correlation.dat',hist,fmt='%1.4e', delimiter='\t')
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
		# cs2num = (hist[i][1]*hist[i][6] - hist[i][0]*hist[i][7])
		# cs2den = ( (hist[i][6])**2 - hist[i][5]*hist[i][7] )
		cv2 = (hist[i][1]*hist[i][7] - hist[i][2]*hist[i][6])/( (hist[i][7])**2 - hist[i][6]*hist[i][8] )
		# cv2num = (hist[i][0]*hist[i][3] - hist[i][1]*hist[i][2])
		# cv2den = ( (hist[i][6])**2 - hist[i][5]*hist[i][7] )

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

np.savetxt('./postprocess/soundvels.dat',hist,fmt='%1.4e', delimiter='\t')

import matplotlib.pyplot as plt
plt.title('cs2')
plt.plot(rlist,cs2list,'-o',markersize=3)
plt.savefig('cs2.jpg')
plt.show()
plt.clf()

plt.title('cs2num')
plt.plot(rlist,cs2numlist,'-o',markersize=3)
plt.savefig('cs2num.jpg')
plt.show()
plt.clf()

plt.title('cs2den')
plt.plot(rlist,cs2denlist,'-o',markersize=3)
plt.savefig('cs2den.jpg')
plt.show()
plt.clf()

plt.title('cv2')
plt.plot(rlist,cv2list,'-o',markersize=3)
plt.savefig('cv2.jpg')
plt.show()
plt.clf()

plt.title('cv2num')
plt.plot(rlist,cv2numlist,'-o',markersize=3)
plt.savefig('cv2num.jpg')
plt.show()
plt.clf()

plt.title('cv2den')
plt.plot(rlist,cv2denlist,'-o',markersize=3)
plt.savefig('cv2den.jpg')
plt.show()
plt.clf()
