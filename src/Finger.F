#include "fintrf.h" 
!===================================================================================================
#if 0
! FDupdate.F 
! .F file needs to be preprocessed to generate .f90 equivalent
! 
#endif
! Finger.f90 
! Kevin Coltin 

module Finger 

public 

! Structure containing parameters defining the finger domain. These are used by most functions 
! in FDupdate and Laplacian. 
! theta is a 2 x n-1 array, where n = min(ymax(1),ymax(end))-1. The second dimension is indexed
! from 2:n. Then, theta(1,iy) and theta(2,iy) are the values of theta at points (1,iy) and 
! (end,iy), respectively. 
type domain_params 
    integer, allocatable, dimension(:) :: ymax, xmin, xmax 
    double precision, allocatable, dimension(:,:) :: theta 
end type domain_params  
 
! Indices of the input args, as they appear in "FDupdate.mex*". 
integer, parameter :: C_INDEX = 1, T0_INDEX = 2, NT_INDEX = 3, DT_INDEX = 4, H_INDEX = 5, &
                        YMAX_INDEX = 6, XMIN_INDEX = 7, XMAX_INDEX = 8, A_INDEX = 9, B_INDEX = 10, &
                        D_INDEX = 11, GAMMA_INDEX = 12, THETA_INDEX = 13
integer, parameter :: NARGIN = 13  

! Constants defining the boundary conditions at different parts of the domain. Used by the  
! tridiagonal solvers and by Laplacian.F. 
integer, parameter :: BC_NEUMANN = 0, BC_DIRICHLET = 1, BC_OBLIQUE = 2, BC_MIXED = 3

end module Finger 








