import numpy as np



stride = 0.5 # In the units of Mpc
max_r = 70.0



tf = np.asarray(np.loadtxt('./data/terf.dat'))

tf = np.sort(tf,axis=0)

print tf[0:10], '\t'



for r_no in np.arange(1,max_r/stride + stride):
	P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0																																												
	for i in range(2000*r_no+0,2000*r_no+2000,2):
		# print i
		# P = tf[i]*tf[i+1]
		P_AD    += tf[i][2]*tf[i+1][2]
		P_AT    += tf[i][4]*tf[i+1][1]
		P_DD    += tf[i][2]*tf[i+1][2]
		P_DT    += tf[i][2]*tf[i+1][3]
		P_TT    += tf[i][3]*tf[i+1][3]
		d2_P_DD += tf[i][2]*tf[i+1][5]
		d2_P_DT += tf[i][5]*tf[i+1][3]
		d2_P_TT += tf[i][3]*tf[i+1][6]

