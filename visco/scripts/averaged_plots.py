import numpy as np
import os
import re

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

r = np.loadtxt('../data/varying_Om/Om=0.2/rlist')

dirs = os.listdir('../data/varying_Om')
dirs.sort(key=natural_keys)

averaged_cs2 = []
averaged_cv2 = []
for drc in dirs:
	data = np.loadtxt('../data/varying_Om/'+drc+'/soundvels.dat')
	cs2_avg = np.loadtxt('../data/varying_Om/'+drc+'/cs2.dat')
	cs2_avg = np.average(cs2_avg[20:134])
	cv2_avg = np.loadtxt('../data/varying_Om/'+drc+'/cv2.dat')
	cv2_avg = np.average(cv2_avg[20:134])

	# print cs2_avg
	averaged_cs2.append(cs2_avg)
	averaged_cv2.append(cv2_avg)


np.savetxt('../postprocess/averaged_cs2.dat',averaged_cs2,fmt='%1.4e', delimiter='\t')
np.savetxt('../postprocess/averaged_cv2.dat',averaged_cv2,fmt='%1.4e', delimiter='\t')

import matplotlib.pyplot as plt
plt.title('cs2')
plt.plot(r,cs2list,'-o',markersize=3)
plt.savefig('cs2.jpg')
plt.show()
plt.clf()

plt.title('cs2num')
plt.plot(r,cs2numlist,'-o',markersize=3)
plt.savefig('cs2num.jpg')
plt.show()
plt.clf()