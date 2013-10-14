% compile.m 
% Script to compile the Fortran modules called by MATLAB for solving
% fingerprint PDEs. 

mex -f /home/kcoltin/Code/mexopts_mathpost.sh -g -largeArrayDims ...
	Finger.F Laplacian.F FDupdate.F -L/usr/lib -llapack -lblas ...
	 mxFDupdate.F -o FDupdate 

