import numpy as np
import os

R = 10.0
Box_size = 60.0


def generate_pair(r):

	r1 = np.random.rand(3)
	r2 = np.random.rand(3)
	d  = np.linalg.norm(r1-r2)


	r1 = (r/d)*r1
	r2 = (r/d)*r2

	# print r1,r2
	# d= np.linalg.norm(r1-r2)
	print np.linalg.norm(r1-r2)
	return ([r1,r2])

print generate_pair(10.0)
# for i in range(len(r1)):
# 	if r[i]< Box_size:
# 		rlist.append(r1)
