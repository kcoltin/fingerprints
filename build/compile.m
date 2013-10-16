% compile.m 
% Script to compile the Fortran modules called by MATLAB for solving
% fingerprint PDEs. 

mex -f mexopts_fprints.sh -g -largeArrayDims Finger.F Laplacian.F ...
	FDupdate.F -L/usr/lib -llapack -lblas mxFDupdate.F -o FDupdate 

