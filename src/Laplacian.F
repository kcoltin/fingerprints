#include "fintrf.h" 
!===================================================================================================
#if 0
! Laplacian.F 
! .F file needs to be preprocessed to generate .f90 equivalent
! 
#endif
! Laplacian.f90 
! Kevin Coltin 

module Laplacian 

use Finger 

private 
public divdotdel, del2, xfirstdiff, yfirstdiff, xseconddiff, yseconddiff 

contains 

! -------------------------------------------------------------------------------------------------
! Returns del * (d del w), i.e. the divergence times d * del w. 
! -------------------------------------------------------------------------------------------------
function divdotdel (W, domain, h, d) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W  
    double precision, dimension(size(W,1),size(W,2)) :: divdotdel 
    double precision, dimension(size(W,1),size(W,2)), intent(in) :: d 
    double precision, intent(in) :: h
    type(domain_params), intent(in) :: domain 
    
    divdotdel = d * del2(W, domain, h) & 
                    + xfirstdiff(d, domain) / (2. * h) * xfirstdiff(W, domain) / (2. * h) & 
                    + yfirstdiff(d, domain) / (2. * h) * yfirstdiff(W, domain) / (2. * h)
    return 
end function divdotdel 

! -------------------------------------------------------------------------------------------------
! Computes the Laplacian del^2 w of a quantity w on a finger domain, 
! while maintaining mixed Dirichlet and no-flux boundary conditions: 
!    (n * del) w = 0 on the sides of the finger domain, and 
! w = 0 on the top and bottom.  
! 
! Arguments: 
! W should be a morphogen concentration, either u or v. It should be of
!    dimension NSTEPSX+1 x NSTEPSY+1. 
! domain is a structure containing YMAX, XMIN, and XMAX, the vectors of values indicating the 
! limits of the domain. 
! h is the step size in each spatial dimention. 
! -------------------------------------------------------------------------------------------------
function del2 (W, domain, h) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: del2 
    double precision, intent(in) :: h 
    type(domain_params), intent(in) :: domain 
    
    del2 = xseconddiff(W, domain) / h**2 + yseconddiff(W, domain) / h**2
    return 
end function del2 


! -------------------------------------------------------------------------------------------------
! Computes the central first difference in the x direction of a quantity w on a finger domain.  
! -------------------------------------------------------------------------------------------------
function xfirstdiff (W, domain, bc) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: xfirstdiff 
    type(domain_params), intent(in) :: domain 
    integer, intent(in), optional :: bc 
    double precision, dimension(size(W,1),size(W,2)) :: xp1, xm1 
    integer :: ysplit 
    
    ysplit = min(domain%ymax(1), domain%ymax(size(W,1))) 
    
    ! "x plus 1" and "x minus 1" 
    if (.not. present(bc) .or. (present(bc) .and. bc == BC_MIXED)) then 
        xp1 = val_right (W, domain%xmin, domain%xmax, domain%theta(2,:), ysplit)  
        xm1 = val_left (W, domain%xmin, domain%xmax, domain%theta(1,:), ysplit) 
    else ! Neumann b.c.'s 
        xp1 = val_right_nbc (W, domain%xmin, domain%xmax)  
        xm1 = val_left_nbc (W, domain%xmin, domain%xmax) 
    endif 
    
    xfirstdiff = xp1 - xm1 
    return 
end function xfirstdiff 


! -------------------------------------------------------------------------------------------------
! Computes the central first difference in the y direction of a quantity w on a finger domain. 
! -------------------------------------------------------------------------------------------------
function yfirstdiff (W, domain, bc) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: yfirstdiff 
    type(domain_params), intent(in) :: domain 
    integer, intent(in), optional :: bc 
    double precision, dimension(size(W,1),size(W,2)) :: yp1, ym1 
    
    if (.not. present(bc) .or. (present(bc) .and. bc == BC_MIXED)) then 
        yp1= val_top (W, domain%ymax) ! "y plus 1"
        ym1 = val_bottom (W, domain%ymax) ! "y minus 1" 
    else ! Neumann b.c.'s 
        yp1= val_top_nbc (W, domain%ymax) ! "y plus 1"
        ym1 = val_bottom_nbc (W, domain%ymax) ! "y minus 1" 
    endif 
    
    yfirstdiff = yp1 - ym1 
    return 
end function yfirstdiff 


! -------------------------------------------------------------------------------------------------
! Computes the central second difference in the x direction of a quantity w on a finger domain.  
! -------------------------------------------------------------------------------------------------
function xseconddiff (W, domain) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: xseconddiff 
    type(domain_params), intent(in) :: domain 
    double precision, dimension(size(W,1),size(W,2)) :: xp1, xm1 
    integer :: ysplit 
    
    ysplit = min(domain%ymax(1), domain%ymax(size(W,1))) 
    
    ! "x plus 1" and "x minus 1" 
    xp1 = val_right (W, domain%xmin, domain%xmax, domain%theta(2,:), ysplit)  
    xm1 = val_left (W, domain%xmin, domain%xmax, domain%theta(1,:), ysplit) 
    
    xseconddiff = xp1 - 2.*W + xm1 
    return 
end function xseconddiff 


! -------------------------------------------------------------------------------------------------
! Computes the central second difference in the y direction of a quantity w on a finger domain. 
! -------------------------------------------------------------------------------------------------
function yseconddiff (W, domain) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: yseconddiff 
    type(domain_params), intent(in) :: domain 
    double precision, dimension(size(W,1),size(W,2)) :: yp1, ym1 
    
    yp1= val_top (W, domain%ymax) ! "y plus 1"
    ym1 = val_bottom (W, domain%ymax) ! "y minus 1" 

    yseconddiff = yp1 - 2.*W + ym1 
    return 
end function yseconddiff 



! -------------------------------------------------------------------------------------------------
! Returns the value that is directly to the right of each given value -
! i.e., it returns the w_{x+1} term in the second order finite difference
! formula, while maintaining mixed boundary conditions: 
!    (n * del) w = 0 on the sides of the finger, and 
! w = 0 on the top and bottom of the finger.  
! 
! Arguments: 
! W should be a morphogen concentration, either u or v. It should be of
!    dimension NSTEPSX+1 x NSTEPSY+1 and should equal zero at all points
!    outside the boundary of the domain. 
! YMAX is the vector of values indicating the limits of the domain. 
! dx and dy are the step sizes in the x- and y-dimensions, respectively. 
! 
! On return, xp1(x,y) ("x plus 1") is the value (either an actual value or
! a value interpolated by the boundary conditions) at point (x+1,y). 
! -------------------------------------------------------------------------------------------------
function val_right (W, xmin, xmax, theta, ysplit) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_right  
    double precision, dimension(size(W,1),size(W,2)) :: xp1 
    integer, dimension(size(W,2)), intent(in) :: xmin, xmax 
    integer, intent(in) :: ysplit 
    double precision, dimension(2:(ysplit-1)), intent(in) :: theta  
    integer :: iy 
    
    ! Solve first row, with Neumann BC's  
    xp1(xmin(1):(xmax(1)-1),1) = W((xmin(1)+1):xmax(1),1) 
    xp1(xmax(1),1) = W(xmax(1)-1,1) 

    ! Solve for each row in rectangle domain, with oblique BC's 
    do iy = 2, ysplit - 1  
    
        xp1(xmin(iy):(xmax(iy)-1),iy) = W((xmin(iy)+1):xmax(iy),iy) 
        xp1(xmax(iy),iy) = W(xmax(iy)-1,iy) + (W(xmax(iy),iy+1) - W(xmax(iy),iy-1)) * tan(theta(iy))
    
    enddo 
    
    ! Solve for each row in semiellipse domain 
    do iy = ysplit, size(W, 2) 
    
        xp1(xmin(iy):(xmax(iy)-1),iy) = W((xmin(iy)+1):xmax(iy),iy) 

        ! For those points where w(n+1) is outside the domain, do nothing: we are "adding" the value 0
        ! on points outside the domain in conjunction with the homogeneous Dirichlet B.C. 

    enddo 
    
    val_right = xp1 
    return 
end function val_right 


! -------------------------------------------------------------------------------------------------
! Returns the value that is directly to the left of each given value -
! i.e., it returns the w_{x-1} term in the second order finite difference
! formula, while maintaining mixed boundary conditions: 
!    (n * del) w = 0 on the sides of the finger, and 
! w = 0 on the top and bottom of the finger.  
! 
! Arguments: 
! W should be a morphogen concentration, either u or v. It should be of
!    dimension NSTEPSX+1 x NSTEPSY+1 and should equal zero at all points
!    outside the boundary of the domain. 
! YMAX is the vector of values indicating the limits of the domain. 
! dx and dy are the step sizes in the x- and y-dimensions, respectively. 
! 
! On return, xm1(x,y) ("x minus 1") is the value (either an actual value or
! a value interpolated by the boundary conditions) at point (x-1,y). 
! -------------------------------------------------------------------------------------------------
function val_left (W, xmin, xmax, theta, ysplit) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_left  
    double precision, dimension(size(W,1),size(W,2)) :: xm1 
    integer, dimension(size(W,2)), intent(in) :: xmin, xmax 
    integer, intent(in) :: ysplit 
    double precision, dimension(2:(ysplit-1)), intent(in) :: theta 
    integer :: iy 
    
    ! Solve first row, with Neumann BC's  
    xm1((xmin(1)+1):xmax(1),1) = W(xmin(1):(xmax(1)-1),1) 
    xm1(xmin(1),1) = W(xmin(1)+1,1) 
    
    ! Solve for each row in rectangle domain, with oblique BC's 
    do iy = 2, ysplit - 1 
        
        xm1((xmin(iy)+1):xmax(iy),iy) = W(xmin(iy):(xmax(iy)-1),iy) 
        xm1(xmin(iy),iy) = W(xmin(iy)+1,iy) + (W(xmin(iy),iy+1) - W(xmin(iy),iy-1)) * tan(theta(iy)) 

    enddo 
    
    ! Solve for each row in semiellipse domain 
    do iy = ysplit, size(W, 2) 
    
        xm1((xmin(iy)+1):xmax(iy),iy) = W(xmin(iy):(xmax(iy)-1),iy) 

        ! For those points where w(n+1) is outside the domain, do nothing: we are "adding" the value 0
        ! on points outside the domain in conjunction with the homogeneous Dirichlet B.C. 

    enddo 
    
    val_left = xm1 
    return 
end function val_left 


! -------------------------------------------------------------------------------------------------
! Returns the value that is directly above each given value -
! i.e., it returns the w_{y+1} term in the second order finite difference
! formula, while maintaining Dirichlet boundary conditions.
! 
! Arguments: 
! W should be a morphogen concentration, either u or v. It should be of
!    dimension NSTEPSX+1 x NSTEPSY+1 and should equal zero at all points
!    outside the boundary of the domain. 
! YMAX is the vector of values indicating the limits of the domain. 
! dx and dy are the step sizes in the x- and y-dimensions, respectively. 
! 
! On return, yp1(x,y) ("y plus 1") is the value (either an actual value or
! a value interpolated by the boundary conditions) at point (x,y+1). 
! -------------------------------------------------------------------------------------------------
function val_top (W, ymax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_top 
    double precision, dimension(size(W,1),size(W,2)) :: yp1 
    integer, dimension(size(W,1)), intent(in) :: ymax 
    integer :: ix 
    
    ! Solve for each column 
    do ix = 1, size(W, 1) 
        
        yp1(ix,1:(ymax(ix)-1)) = W(ix,2:ymax(ix)) 

        ! No need to add top value: "adding" zero for homogeneous Dirichlet B.C. 
        
    enddo 
    
    val_top = yp1 
    return 
end function val_top 


! -------------------------------------------------------------------------------------------------
! Returns the value that is directly below each given value -
! i.e., it returns the w_{y-1} term in the second order finite difference
! formula, while maintaining Dirichlet boundary conditions.
! 
! Arguments: 
! W should be a morphogen concentration, either u or v. It should be of
!    dimension NSTEPSX+1 x NSTEPSY+1 and should equal zero at all points
!    outside the boundary of the domain. 
! YMAX is the vector of values indicating the limits of the domain. 
! dx and dy are the step sizes in the x- and y-dimensions, respectively. 
! 
! On return, ym1(x,y) ("y minus 1") is the value (either an actual value or
! a value interpolated by the boundary conditions) at point (x,y-1). 
! -------------------------------------------------------------------------------------------------
function val_bottom (W, ymax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_bottom 
    double precision, dimension(size(W,1),size(W,2)) :: ym1 
    integer, dimension(size(W,1)), intent(in) :: ymax 
    integer :: ix 
    
    ! Solve for each column  
    do ix = 1, size(W, 1)  
        
        ym1(ix,2:ymax(ix)) = W(ix,1:(ymax(ix)-1)) 
        
        ! No need to add bottom row; "adding" zero for homogeneous Dirichlet B.C. 
        
    enddo 
    
    val_bottom = ym1
    return 
end function val_bottom 


 !****************************************************************************** 
 ! The functions below this point are all for "nbc": Neumann Bounday Conditions. 
 ! These are used just for computing the derivatives of d, since it makes more 
 ! sense for d to have Neumann B.C.'s than the mixed B.C.'s used by the rest of
 ! the domain. 
 
 
! -------------------------------------------------------------------------------------------------
function val_right_nbc (W, xmin, xmax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_right_nbc  
    double precision, dimension(size(W,1),size(W,2)) :: xp1 
    integer, dimension(size(W,2)), intent(in) :: xmin, xmax 
    integer :: iy 
    
    do iy = 1, size(W, 2) 
        
        xp1(xmin(iy):(xmax(iy)-1),iy) = W((xmin(iy)+1):xmax(iy),iy) 
        xp1(xmax(iy),iy) = W(xmax(iy)-1,iy) 
        
    enddo 
    
    val_right_nbc = xp1 
    return 
end function val_right_nbc  


! -------------------------------------------------------------------------------------------------
function val_left_nbc (W, xmin, xmax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_left_nbc  
    double precision, dimension(size(W,1),size(W,2)) :: xm1 
    integer, dimension(size(W,2)), intent(in) :: xmin, xmax 
    integer :: iy 
    
    do iy = 1, size(W, 2) 
        
        xm1((xmin(iy)+1):xmax(iy),iy) = W(xmin(iy):(xmax(iy)-1),iy) 
        xm1(xmin(iy),iy) = W(xmin(iy)+1,iy) 

    enddo 
    
    val_left_nbc = xm1 
    return 
end function val_left_nbc  


! -------------------------------------------------------------------------------------------------
function val_top_nbc (W, ymax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_top_nbc 
    double precision, dimension(size(W,1),size(W,2)) :: yp1 
    integer, dimension(size(W,1)), intent(in) :: ymax 
    integer :: ix 
    
    do ix = 1, size(W, 1) 
        
        yp1(ix,1:(ymax(ix)-1)) = W(ix,2:ymax(ix)) 
        yp1(ix,ymax(ix)) = W(ix,ymax(ix)-1) 
        
    enddo 
    
    val_top_nbc = yp1 
    return 
end function val_top_nbc 

! -------------------------------------------------------------------------------------------------
function val_bottom_nbc (W, ymax) 
    implicit none 
    double precision, dimension(:,:), intent(in) :: W 
    double precision, dimension(size(W,1),size(W,2)) :: val_bottom_nbc 
    double precision, dimension(size(W,1),size(W,2)) :: ym1 
    integer, dimension(size(W,1)), intent(in) :: ymax 
    integer :: ix 
    
    do ix = 1, size(W, 1)  
        
        ym1(ix,2:ymax(ix)) = W(ix,1:(ymax(ix)-1)) 
        ym1(ix,1) = W(ix,2) 
        
    enddo 
    
    val_bottom_nbc = ym1
    return 
end function val_bottom_nbc 



end module Laplacian 





