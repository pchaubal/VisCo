import h5py
import numpy as np
import os
snapdir = 'snapdir_135' # Give here the name of the folder of hdf snapshots
savedir = 'binarysnaps'

position =np.empty((0,3))
velocity =np.empty((0,3))
for file in os.listdir(snapdir):
	print 'converting file', file
	f = h5py.File(snapdir+'/'+file, 'r')
	coor = f.get('PartType1/Coordinates').value
	vel = f.get('PartType1/Velocities').value

	position = np.concatenate((position,coor),axis=0)
	velocity = np.concatenate((velocity,vel),axis=0)


position = np.insert(position,0,position.shape[0],axis=0)
if not os.path.exists(savedir):
	os.makedirs(savedir)

fileobj1 = open(savedir+'/'+'positions.bin', mode='wb')
fileobj2 = open(savedir+'/'+'velocities.bin', mode='wb')
position.T.tofile(fileobj1)
velocity.tofile(fileobj2)
fileobj1.close()
fileobj2.close()

print position[0:3]