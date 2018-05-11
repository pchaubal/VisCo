import numpy as np



stride = 0.5 # In the units of Mpc
max_r = 70.0



tf = np.asarray(np.loadtxt('./data/terf.dat'))

tf = np.sort(tf,axis=0)

print tf[0:10], '\t'


r,cs2,cv2=[],[],[]
for r_no in np.arange(1,max_r/stride + stride):
	P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0	
	print r_no																																											
	for i in range(int(2000*r_no+0),int(2000*r_no+2000),2):
		# print i
		# P = tf[i]*tf[i+1]
		P_AD    += tf[i][4]*tf[i+1][2]
		P_AT    += tf[i][4]*tf[i+1][1]
		P_DD    += tf[i][2]*tf[i+1][2]
		P_DT    += tf[i][2]*tf[i+1][3]
		P_TT    += tf[i][3]*tf[i+1][3]
		d2_P_DD += tf[i][2]*tf[i+1][5]
		d2_P_DT += tf[i][5]*tf[i+1][3]
		d2_P_TT += tf[i][3]*tf[i+1][6]

	# Normalize the tertiary fields now
	P_AD    = P_AD/2000.0    
	P_AT    = P_AT/2000.0    
	P_DD    = P_DD/2000.0    
	P_DT    = P_DT/2000.0    
	P_TT    = P_TT/2000.0    
	d2_P_DD = d2_P_DD/2000.0 
	d2_P_DT = d2_P_DT/2000.0 
	d2_P_TT = d2_P_TT/2000.0 

	cs_sq = (P_AT*d2_P_DT - P_AD*d2_P_TT)/(d2_P_DT**2 - d2_P_DD*d2_P_TT)
	cv_sq = (P_AD*d2_P_DT - P_AT*d2_P_DD)/(d2_P_DT**2 - d2_P_DD*d2_P_TT)

	r.append(r_no*stride)
	cs2.append(cs_sq)
	cv2.append(cv_sq)

import matplotlib.pyplot as plt
plt.title('cs2')
plt.plot(r,cs2,'-o',markersize=2)
plt.savefig('./postprocess/cs2.jpg')
plt.show()
plt.clf()

plt.title('cv2')
plt.plot(r,cv2,'-o',markersize=2)
plt.savefig('./postprocess/cv2.jpg')
plt.show()
plt.clf()


