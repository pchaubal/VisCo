import numpy as np
import os

R = 10.0
Box_size = 60.0
binsize = 0.5 # In the units of Mpc
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

print generate_pair(7.0)


r_list = []
r_pairs = []
for bin_no in np.arange(max_r/binsize):
	print bin_no*binsize, (bin_no+1)*binsize
	for i in range(1000):
		r = 0.5*binsize*np.random.randn() + (bin_no*binsize +binsize*0.5)
		r_pairs.append(generate_pair(r))
	r_list.append([[bin_no*binsize, (bin_no+1)*binsize],r_pairs])
	del r_pairs[:]

print len(r_list)
