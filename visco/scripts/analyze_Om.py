import numpy as np
import os

for drc in os.listdir('../data/varying_Om'):
	print '\n',drc
	os.system('python ../soundvels.py ../data/varying_Om/'+ drc)
