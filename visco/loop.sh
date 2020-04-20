#!/bin/bash
echo "Starting VisCo for multiple files ..."
make clean
make
for i in {0..9}
# rm -r ./data/data00$i
# rm ./data/terf00$i.dat
  do
mkdir ./data/data00$i
mpirun -np 6 ./execute ./snapshots-32/snapshot_00$i 
mv ./data/tf* ./data/data00$i/
cat ./data/data00$i/tf* > ./data/terf00$i.dat
echo "Calculation of snapshot_00$i is done"
cp ./data/terf00$i.dat ./data/terf.dat
python viscosities.py
mv ./postprocess/datafile.txt ./postprocess/datafile00$i.dat  
echo "Calculation of visocsity parameters from terf00$i is done and stored in POSTPROCESS"

done

for i in {10..41}
# rm -r ./data/data00$i
# rm ./data/terf00$i.dat
  do
mkdir ./data/data0$i
mpirun -np 6 ./execute ./snapshots-32/snapshot_0$i 
mv ./data/tf* ./data/data0$i/
cat ./data/data0$i/tf* > ./data/terf0$i.dat
echo "Calculation of snapshot_0$i is done"
cp ./data/terf0$i.dat ./data/terf.dat
python viscosities.py
mv ./postprocess/datafile.txt ./postprocess/datafile0$i.dat  
echo "Calculation of visocsity parameters from terf0$i is done and stored in POSTPROCESS"

done
