make clean 
rm *.so
rm *.mat
rm MTvec.txt
touch MTvec.txt

gfortran -c constants.f90 
make
mv newkry.* newkry.so
#python3 krylsolve.py