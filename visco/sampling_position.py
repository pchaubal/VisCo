import numpy as np
import os

R = 10.0
Box_size = 60.0
stride = 0.5 # In the units of Mpc
max_r = 70.0

def generate_pair(distance):

	found_pair = False
	while (found_pair==False):
		r1 = np.random.rand(3)
		r2 = np.random.rand(3)
		d  = np.linalg.norm(r1-r2)

		r1 = (distance/d)*r1
		r2 = (distance/d)*r2

		if (np.max(r1)<Box_size) and (np.max(r2)<Box_size):
			found_pair = True
			# print np.linalg.norm(r1-r2)

	return ([r1,r2])

# print generate_pair(7.0)


r_list = []
r_pairs = []
positions=[]
for r_no in np.arange(1,max_r/stride + stride):
	print 'Generating pairs at:%s'%(r_no*stride)
	for i in range(2000):
		r = r_no*stride
		pairs = generate_pair(r)
		r_pairs.append(pairs)
		positions.append(pairs)
	r_list.append([r_no*stride,r_pairs])
	del r_pairs[:]

positions = np.asarray(positions).reshape(560000,3)
np.savetxt('sampling_point.dat',positions,fmt='%1.4e', delimiter='\t')
# print len(r_list), np.asarray(positions).reshape(560000,3)
