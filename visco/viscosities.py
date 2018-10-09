import numpy as np
import matplotlib.pylab as plt

stride = 0.5 # In the units of Mpc
max_r = 70.0 # max correlation length in Mpc

tf = np.asarray(np.loadtxt('/home/prakrut/codes/visco/terf.dat'))
tf = np.sort(tf,axis=0)
print tf.shape # There are 2000*140 paris = 560000 points
# The first 2000 pairs (=4000 points) correspond to a correlation length of 0.5 Mpc
# The next 2000 pairs to correlation length of 1.0 Mpc and so on

# Vectorized part
def moving_avg(list,n):
	avg = [np.mean(list[i:i+n]) for i in range(0,len(list),n)]
	return np.asarray(avg)

r = np.arange(stride,max_r+stride,stride)

# The first term i.e. tf[0::2,i] is taking every 2nd element starting from 0 of the list
# The second term i.e. tf[1::2,i] is taking every 2nd element but starting from 1 of the list
# So tf[0::2,i]*tf[1::2,j] is an array of length 280000 with 
# every 2000 elements corresponding to a unique correlation length
# The function moving_avg takes the mean of every 2000 elements of the list
P_AD_v    = moving_avg(tf[0::2,4]*tf[1::2,2],2000)
P_AT_v    = moving_avg(tf[0::2,4]*tf[1::2,3],2000)
P_DD_v    = moving_avg(tf[0::2,2]*tf[1::2,2],2000)
P_DT_v    = moving_avg(tf[0::2,2]*tf[1::2,3],2000)
P_TT_v    = moving_avg(tf[0::2,3]*tf[1::2,3],2000)
d2_P_DD_v = moving_avg(tf[0::2,2]*tf[1::2,5],2000)
d2_P_DT_v = moving_avg(tf[0::2,5]*tf[1::2,3],2000)
d2_P_TT_v = moving_avg(tf[0::2,3]*tf[1::2,6],2000)

P_AA_v = moving_avg(tf[0::2,4]*tf[1::2,4],2000)

cs2_v_num = (P_AT_v*d2_P_DT_v - P_AD_v*d2_P_TT_v)
cs2_v_den = (d2_P_DT_v**2.0 - d2_P_DD_v*d2_P_TT_v)
cs2_v = cs2_v_num/cs2_v_den

cv2_v_num = (P_AD_v*d2_P_DT_v - P_AT_v*d2_P_DD_v)
cv2_v_den = (d2_P_DT_v**2.0 - d2_P_DD_v*d2_P_TT_v)
cv2_v = cv2_v_num/cv2_v_den


# Same but looping
# cs2,cv2=[],[]
# for j in range(len(r)):
# 	P_AD, P_AT, P_DD, P_DT, P_TT, d2_P_DD, d2_P_DT, d2_P_TT = 0,0,0,0,0,0,0,0	
# 	for pt_ind in range(0,4000,2):
# 		# print i
# # 		# P = tf[i]*tf[i+1]
# 		P_AD    += tf[j*4000+ pt_ind][4] * tf[j*4000+ (pt_ind+1)][2]
# 		P_AT    += tf[j*4000+ pt_ind][4] * tf[j*4000+ (pt_ind+1)][3]
# 		P_DD    += tf[j*4000+ pt_ind][2] * tf[j*4000+ (pt_ind+1)][2]
# 		P_DT    += tf[j*4000+ pt_ind][2] * tf[j*4000+ (pt_ind+1)][3]
# 		P_TT    += tf[j*4000+ pt_ind][3] * tf[j*4000+ (pt_ind+1)][3]
# 		d2_P_DD += tf[j*4000+ pt_ind][2] * tf[j*4000+ (pt_ind+1)][5]
# 		d2_P_DT += tf[j*4000+ pt_ind][5] * tf[j*4000+ (pt_ind+1)][3]
# 		d2_P_TT += tf[j*4000+ pt_ind][3] * tf[j*4000+ (pt_ind+1)][6]

# # 	# Normalize the tertiary fields now
# 	P_AD    = P_AD/2000.0    
# 	P_AT    = P_AT/2000.0    
# 	P_DD    = P_DD/2000.0    
# 	P_DT    = P_DT/2000.0    
# 	P_TT    = P_TT/2000.0    
# 	d2_P_DD = d2_P_DD/2000.0 
# 	d2_P_DT = d2_P_DT/2000.0 
# 	d2_P_TT = d2_P_TT/2000.0 


# 	try:
# 		cs_sq = (P_AT*d2_P_DT - P_AD*d2_P_TT)/(d2_P_DT**2 - d2_P_DD*d2_P_TT)
# 		cv_sq = (P_AD*d2_P_DT - P_AT*d2_P_DD)/(d2_P_DT**2 - d2_P_DD*d2_P_TT)
# 	except ZeroDivisionError:
# 		cs_sq = 0
# 		cv_sq = 0

# 	cs2.append(cs_sq)
# 	cv2.append(cv_sq)

###############################################################################
# Plots
# Use the first two to compare the looped and vectorized versions
# Ad

# plt.title('cv2')
# plt.plot(r,cv2,'-o',markersize=2,label='Looping')
# plt.plot(r,cv2_v,'-o',markersize=2,label='Vectorized')
# plt.plot(r,cv2_v_num/r**2,label='cv2 numerator/r2')
# plt.plot(r,cv2_v_num*r**2,label='cv2 numerator*r2')
# plt.plot(r,cv2_v_num,label='numerator')
# plt.plot(r,cv2_v_den,label ='cv2 denominator')

# plt.plot(r[5:100],(cs2_v_num)[5:100],label ='cs2 numerator')
# plt.plot(r[5:100],(cs2_v_num*r**2.0)[5:100],label ='cs2 numerator*r2') # plotting only first 5 to 100 points
# plt.plot(r[5:100],(-cs2_v_num/r**2.0)[5:100],label ='cs2 numerator/r2')
# plt.plot(r[5:100],(cs2_v_den)[5:100],'-o',markersize=2,label ='cs2 denominator')


# plt.plot(r,P_AA_v, label='P_AA')
# plt.plot(r,P_AA_v/r**4, label='P_AA/r**4')
# plt.plot(r[20:100],(P_AA_v/r**2.6)[20:100], label='P_AA/r**4')
plt.plot(r[20:80],(d2_P_DT_v**2.0*r**3)[20:80], label='P_AA/r**4')

# plt.xlim(5,30)
plt.legend()
# plt.savefig('cv2.jpg')
plt.show()
plt.clf()
