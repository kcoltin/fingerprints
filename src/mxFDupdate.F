#include "fintrf.h" 
!===================================================================================================
#if 0
! mxFDupdate.F 
! .F file needs to be preprocessed to generate .f90 equivalent
! 
#endif
!
! mxFDupdate.f90
! Kevin Coltin 
!
! "Mex finite difference update." Contains functions for updating several steps of the finite 
! difference scheme for the Gierer-Meinhardt system for fingerprint pattern formation. This is 
! designed to be called by the MATLAB function fingerprintSim.m, to speed up the computation. 
! 
! INSTRUCTIONS FOR CALLING FROM MATLAB: 
! The syntax is 
! C = FDupdate (C, time, nt, dt, h, YMAX, XMIN, XMAX, a, b, d, gamma, theta);

! Gateway routine 
subroutine mexFunction (nlhs, plhs, nrhs, prhs) 
use FDupdate
use Finger 

implicit none 

mwPointer plhs(*), prhs(*) 
integer nlhs, nrhs 
mwSize, parameter :: NDIM = 3 ! number of dimensions 

! Type mwSize necessary to avoid seg fault
mwSize, parameter :: MW_ONE = 1 
mwSize, dimension(NDIM) :: dims ! dimensions of C 
mwSize csize, dsize, thetasize, nx, ny 

mwPointer C_pr_rhs, C_pr_lhs, t0_pr, nt_pr, dt_pr, h_pr, ymax_pr, xmin_pr, xmax_pr, &
               a_pr, b_pr, d_pr, gamma_pr, theta_pr 

double precision, dimension(:,:,:), allocatable :: C 
double precision :: t0, dt, h, a, b, gamma 
double precision, dimension(:,:), allocatable :: d 
integer :: nt, info, ysplit 
type(domain_params) :: domain 
double precision :: temp1 ! for holding real values until they can be converted to int
double precision, dimension(:), allocatable :: temp2
integer :: stat
character(len=100) :: msg 

! Declarations of built-in MEX functions used
mwPointer mxCreateNumericArray, mxGetPr, mxGetDimensions 
integer mxClassIDFromClassName 
  
! Check for correct number of arguments 
if (nrhs /= NARGIN) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate:nInput', 'Invalid number of inputs.')
elseif (nlhs /= 1) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate:nOutput', 'One output required.')
endif

! Get size of C, and allocate C and ymax 
call mxCopyPtrToInteger4 (mxGetDimensions (prhs(C_INDEX)), dims, NDIM)
csize = dims(1) * dims(2) * 2 ! third dimension is U-V, so =2  
nx = dims(1) 
ny = dims(2) 
dsize = dims(1) * dims(2) 

! Check dimensions 
if (dims(1) < 5) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'x-dimension of C is too small.') 
endif 
if (dims(2) < 5) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'y-dimension of C is too small.') 
endif 
if (dims(3) /= 2) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Third dimension of C must be of size 2.') 
endif 

! Allocate arrays
allocate (C(dims(1), dims(2), 2), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating C.')
endif 
allocate (d(dims(1), dims(2)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating d.')
endif 
allocate (domain%ymax(dims(1)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating domain%ymax.')
endif 
allocate (domain%xmin(dims(2)), stat=stat)
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating domain%xmin.')
endif 
allocate (domain%xmax(dims(2)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating domain%xmax.')
endif 
! note: theta and temp2 are allocated later. 

! Create 3-d array for the return argument 
plhs(1) = mxCreateNumericArray(NDIM, dims, mxClassIDFromClassName('double'), 0) 
C_pr_lhs = mxGetPr(plhs(1)) 

! Convert input arguments into correct fortran data types 
C_pr_rhs = mxGetPr(prhs(C_INDEX)) 
call mxCopyPtrToReal8 (C_pr_rhs, C, csize)  

t0_pr = mxGetPr(prhs(T0_INDEX)) 
call mxCopyPtrToReal8 (t0_pr, t0, MW_ONE) 

nt_pr = mxGetPr(prhs(NT_INDEX)) 
call mxCopyPtrToReal8 (nt_pr, temp1, MW_ONE) 
nt = int(temp1) 

dt_pr = mxGetPr(prhs(DT_INDEX)) 
call mxCopyPtrToReal8 (dt_pr, dt, MW_ONE) 

h_pr = mxGetPr(prhs(H_INDEX)) 
call mxCopyPtrToReal8 (h_pr, h, MW_ONE) 

ymax_pr = mxGetPr(prhs(YMAX_INDEX)) 
allocate (temp2(dims(1)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating temp2.')
endif 
call mxCopyPtrToReal8 (ymax_pr, temp2, nx) 
domain%ymax = int(temp2) 
deallocate(temp2, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating temp2.')
endif 

xmin_pr = mxGetPr(prhs(XMIN_INDEX)) 
allocate (temp2(dims(2)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating temp2.')
endif 
call mxCopyPtrToReal8 (xmin_pr, temp2, ny) 
domain%xmin = int(temp2) 
deallocate(temp2, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating temp2.')
endif 

xmax_pr = mxGetPr(prhs(XMAX_INDEX)) 
allocate (temp2(dims(2)), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating temp2.')
endif 
call mxCopyPtrToReal8 (xmax_pr, temp2, ny) 
domain%xmax = int(temp2) 
deallocate(temp2, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating temp2.')
endif 

a_pr = mxGetPr(prhs(A_INDEX)) 
call mxCopyPtrToReal8 (a_pr, a, MW_ONE) 

b_pr = mxGetPr(prhs(B_INDEX)) 
call mxCopyPtrToReal8 (b_pr, b, MW_ONE) 

d_pr = mxGetPr(prhs(D_INDEX)) 
call mxCopyPtrToReal8 (d_pr, d, dsize)

gamma_pr = mxGetPr(prhs(GAMMA_INDEX)) 
call mxCopyPtrToReal8 (gamma_pr, gamma, MW_ONE) 

theta_pr  = mxGetPr(prhs(THETA_INDEX))
ysplit = min(domain%ymax(1), domain%ymax(size(domain%ymax))) 
allocate (domain%theta(2, 2:ysplit-1), stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error allocating domain%theta.')
endif 
thetasize = 2 * (ysplit - 2)
call mxCopyPtrToReal8 (theta_pr, domain%theta, thetasize) 


! Call the computational routine 
call update (C, t0, nt, dt, h, domain, a, b, d, gamma, info) 

! Load the data into C_pr, which is the output to MATLAB. 
call mxCopyReal8ToPtr (C, C_pr_lhs, csize) 

! Deallocate arrays 
deallocate (C, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating C.')
endif 
deallocate (d, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating d.')
endif 
deallocate (domain%ymax, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating domain%ymax.')
endif 
deallocate (domain%xmin, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating domain%xmin.')
endif 
deallocate (domain%xmax, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating domain%xmax.')
endif 
deallocate (domain%theta, stat=stat) 
if (stat /= 0) then 
    call mexErrMsgIdAndTxt ('MATLAB:FDupdate', 'Error deallocating domain%theta.')
endif 

! Check for errors or warnings from "update" 
if (info < 0) then 
    write (msg, *) 'Input argument number ', -info, ' had an illegal value.' ! Note 1
    call mexWarnMsgIdAndTxt ('MATLAB:FDupdate', msg) 
elseif (info > 0) then 
    call mexWarnMsgIdAndTxt ('MATLAB:FDupdate', "Error returned by LAPACK routine DGTSV: There " // &
       "was a singular matrix in one or more of the tridiagonal systems solved in the ADI method.") 
endif 

return 
end subroutine mexFunction  


! NOTES 
! 1. If you increase the length of this string, be sure to increase the "len" attribute of msg 
!        when it is declared at the top! 












