#include "fintrf.h" 
!===================================================================================================
#if 0
! FDupdate.F 
! .F file needs to be preprocessed to generate .f90 equivalent
! 
#endif
! FDupdate.f90 
! Kevin Coltin 

module FDupdate ! "finite difference update" 

use Finger 
use Laplacian 

private 
public update 

contains 

! -------------------------------------------------------------------------------------------------
! Updates the finite difference scheme for the Gierer-Meinhardt system for nt timesteps. 
! The output is the updated version of C: C(:,:,1) = u, and C(:,:,2) = v. 
! Implements the Peaceman-Rachford ADI method. 
! -------------------------------------------------------------------------------------------------
subroutine update (C, t0, nt, dt, h, domain, a, b, d, gamma, info)
    implicit none 
    
    double precision, dimension(:,:,:), intent(inout) :: C 
    double precision, intent(in) :: t0, dt, h, a, b, gamma 
    double precision, dimension(size(C,1), size(C,2)), intent(in) :: d 
    integer, intent(in) :: nt 
    type(domain_params), intent(in) :: domain 
    integer, intent(out) :: info 
   
    integer :: i 
    integer, parameter :: U = 1, V = 2 ! indices in the third dimension of C 
    double precision, parameter :: MIN_VAL = .001 ! minimum allowable value of U, V 
    double precision, dimension(size(C,1),size(C,2)) :: fval, gval 
    double precision, dimension(size(C,1),size(C,2),2) :: dprime 
    double precision :: t 
   
    info = checkArgs (C, nt, dt, h, domain, d, gamma) 
    if (info /= 0) then 
       return
    endif 
   
   ! Compute partial derivatives of d in the x- and y-directions. 
   ! dprime(:,:,1) is the x derivative; dprime(:,:,2) is y. 
    dprime(:,:,1) = xfirstdiff(d, domain, BC_NEUMANN) / (2. * h) 
    dprime(:,:,2) = yfirstdiff(d, domain, BC_NEUMANN) / (2. * h)
   
    do i = 1, nt 
       
       t = t0 + (i - 1.) * dt 
       
       ! First half-step of the ADI method 
       fval = f (C(:,:,U), C(:,:,V), a, b) 
       gval = g (C(:,:,U), C(:,:,V)) ! Note 1 
       
       call adiStep1 (C(:,:,U), gamma*dt/2.*fval, dt, h, domain, info) 
       call adiStep1 (C(:,:,V), gamma*dt/2.*gval, dt, h, domain, info, d, dprime) 
       
       C = pmax (C, MIN_VAL) 
       
       ! Second half-step of ADI 
       fval = f (C(:,:,U), C(:,:,V), a, b) 
       gval = g (C(:,:,U), C(:,:,V)) 
       
       call adiStep2 (C(:,:,U), gamma*dt/2.*fval, dt, h, domain, info) 
       call adiStep2 (C(:,:,V), gamma*dt/2.*gval, dt, h, domain, info, d, dprime) 
       
       C = pmax (C, MIN_VAL) 
       
    enddo 

end subroutine update


! ------------------------------------------------------------------------------------------------- 
! The function f in the reaction-diffusion system. 
! -------------------------------------------------------------------------------------------------
elemental function f (u, v, a, b) 
    implicit none 
    double precision :: f 
    double precision, intent(in) :: u, v, a, b 
   
    f = a - b * u + u**2 / v 
    return 
end function f 


! -------------------------------------------------------------------------------------------------
! The function g in the reaction diffusion system. 
! -------------------------------------------------------------------------------------------------
elemental function g (u, v) 
    implicit none 
    double precision :: g 
    double precision, intent(in) :: u, v 
   
    g = u**2 - v 
    return 
end function g 



! -------------------------------------------------------------------------------------------------
! First half-step of the ADI method. 
! W is the quantity, either u or v, we are interested in. 
! f is a "function" or other quantity that essentially defines the inhomogeneous part of the PDE. 
!        E.g., it is the function f in the pde w_t - del^2 w = f. 
!        (Technically, f should be dt/2*f; i.e. it should have the timestep built in.) 
! d is an optional parameter, if included, the PDE becomes w_t - d*del^2 w = f. 
! -------------------------------------------------------------------------------------------------
subroutine adiStep1 (W, f, dt, h, domain, info, d, dprime) 
    implicit none 
    double precision, dimension(:,:), intent(inout) :: W 
    double precision, dimension(size(W,1), size(W,2)), intent(in) :: f
    double precision, intent(in) :: dt, h 
    type(domain_params), intent(in) :: domain 
    integer, intent(inout) :: info 
    double precision, dimension(size(W,1), size(W,2)), optional, intent(in) :: d 
    double precision, dimension(size(W,1), size(W,2), 2), optional, intent(in) :: dprime 
   
    double precision, dimension(size(W,1), size(W,2)) :: delta2x, deltax, ff
    integer :: ix, info_tmp
   
   ! Compute first and second derivative approximation in x-direction 
    delta2x = xseconddiff(W, domain) / h**2
    deltax = xfirstdiff(W, domain) / (2.*h) 
   
    if (.not. present(d)) then 
       ff = dt/2.*delta2x + f 
    else
       ff = dt/2.*d*delta2x + dt/2.*dprime(:,:,1)*deltax + f 
    endif 
   
   ! Iterate over each column of the domain 
    do ix = 1, size(W, 1) 
       if (.not. present(d)) then 
           call tridiagSolve (W(ix,1:domain%ymax(ix)), ff(ix,1:domain%ymax(ix)), dt, h, & 
                               BC_DIRICHLET, info_tmp) 
       else 
           call tridiagSolve (W(ix,1:domain%ymax(ix)), ff(ix,1:domain%ymax(ix)), dt, h, &
                               BC_DIRICHLET, info_tmp, d(ix,1:domain%ymax(ix)), &
                               dprime(ix,1:domain%ymax(ix),2)) 
       endif 
       info = max(info, info_tmp) ! Note 4  
    enddo 

end subroutine adiStep1 


! -------------------------------------------------------------------------------------------------
! Second half-step of the ADI method. 
! W is the quantity, either u or v, we are interested in. 
! f is a "function" or other quantity that essentially defines the inhomogeneous part of the PDE. 
!        E.g., it is the function f in the pde w_t - del^2 w = f. 
!        (Technically, f should be dt*f; i.e. it should have the timestep built in.) 
! d is an optional parameter, if included, the PDE becomes w_t - d*del^2 w = f. 
! -------------------------------------------------------------------------------------------------
subroutine adiStep2 (W, f, dt, h, domain, info, d, dprime) 
    implicit none 
    double precision, dimension(:,:), intent(inout) :: W 
    double precision, dimension(size(W,1), size(W,2)), intent(in) :: f
    double precision, intent(in) :: dt, h
    type(domain_params), intent(in) :: domain 
    integer, intent(inout) :: info 
    double precision, dimension(size(W,1), size(W,2)), optional, intent(in) :: d 
    double precision, dimension(size(W,1), size(W,2), 2), optional, intent(in) :: dprime 
   
    double precision, dimension(size(W,1), size(W,2)) :: delta2y, deltay, ff
    integer :: iy, ysplit, inf 
    double precision, dimension(2,2:(min(domain%ymax(1),domain%ymax(size(W,1)))-1)) :: dwdy ! Note 3 
   
   ! Compute first and second derivative approximation in y-direction 
    delta2y = yseconddiff(W, domain) / h**2
    deltay = yfirstdiff(W, domain) / (2.*h)
   
    if (.not. present(d)) then 
       ff = dt/2.*delta2y + f 
    else
       ff = dt/2.*d*delta2y + dt/2.*dprime(:,:,2)*deltay + f 
    endif 
   
   ! Split y into rectangle part and semicircle part 
    ysplit = min(domain%ymax(1), domain%ymax(size(W,1))) ! Note 2
   
   ! dwdy is the vector of values (u_{x,y+1} - u_{x,y-1}) on either side of the finger, at the 
   ! points where oblique boundary conditions apply. dwdy(1,iy) is the first difference at point
   ! (1,iy), and dwdy(2,iy) is the difference at point (NX,iy). 
    dwdy(1,:) = W(1,3:ysplit) - W(1,1:(ysplit-2)) 
    dwdy(2,:) = W(size(W,1),3:ysplit) - W(size(W,1),1:(ysplit-2)) 
   
   ! Solve for bottom row of rectangle 
    if (.not. present(d)) then 
       call tridiagSolve (W(:,1), ff(:,1), dt, h, BC_NEUMANN, inf) 
    else 
       call tridiagSolve (W(:,1), ff(:,1), dt, h, BC_NEUMANN, inf, d(:,1), dprime(:,1,1)) 
    endif 
    info = max(info, inf) ! Note 4  
   
   ! Iterate over rectangle, with oblique BCs 
    do iy = 2, ysplit - 1
       if (.not. present(d)) then 
           call tridiagSolve (W(:,iy), ff(:,iy), dt, h, BC_OBLIQUE, inf, dwdy=dwdy(:,iy), &
                                   theta=domain%theta(:,iy)) 
       else 
           call tridiagSolve (W(:,iy), ff(:,iy), dt, h, BC_OBLIQUE, inf, d(:,iy), dprime(:,iy,1), &
                                   dwdy(:,iy), domain%theta(:,iy)) 
       endif 
       info = max(info, inf) 
    enddo 
   
   ! Iterate over semicircle 
    do iy = ysplit, size(W, 2) 
       if (.not. present(d)) then 
           call tridiagSolve (W(domain%xmin(iy):domain%xmax(iy),iy), & 
                               ff(domain%xmin(iy):domain%xmax(iy),iy), dt, h, BC_DIRICHLET, inf) 
       else 
           call tridiagSolve (W(domain%xmin(iy):domain%xmax(iy),iy), & 
                               ff(domain%xmin(iy):domain%xmax(iy),iy), dt, h, BC_DIRICHLET, inf, & 
                               d(domain%xmin(iy):domain%xmax(iy),iy), &
                               dprime(domain%xmin(iy):domain%xmax(iy),iy,1)) 
       endif 
       info = max(info, inf) 
    enddo 
   
end subroutine adiStep2 



! -------------------------------------------------------------------------------------------------
! Solves the tridiagonal system
!   A w^{n+1} = w^n + rhs, 
! for an update of the ADI method along one line. 
! -------------------------------------------------------------------------------------------------
subroutine tridiagSolve (w, rhs, dt, h, bc, info, d, dprime, dwdy, theta) 
    implicit none 
    double precision, dimension(:), intent(inout) :: w 
    double precision, dimension(size(w)), intent(in) :: rhs 
    double precision, intent(in) :: dt, h
    integer, intent(in) :: bc 
    integer, intent(out) :: info
    double precision, dimension(size(w)), optional, intent(in) :: d, dprime 
    double precision, dimension(2), optional, intent(in) :: dwdy, theta 
   
    if (BC /= BC_OBLIQUE) then 
       if (.not. present(d)) then 
           call tridiagSolveNoD (w, rhs, dt, h, bc, info) 
       else 
           call tridiagSolveD (w, rhs, dt, h, bc, info, d, dprime) 
       endif 
    else 
       if (.not. present(d)) then 
           call tridiagSolveNoD (w, rhs, dt, h, bc, info, dwdy, theta) 
       else 
           call tridiagSolveD (w, rhs, dt, h, bc, info, d, dprime, dwdy, theta) 
       endif 
    endif 
   
end subroutine tridiagSolve 


! -------------------------------------------------------------------------------------------------
! Implements tridiagSolve when the optional argument d is absent. 
! -------------------------------------------------------------------------------------------------
subroutine tridiagSolveNoD (w, rhs, dt, h, bc, info, dwdy, theta) 
    implicit none 
    double precision, dimension(:), intent(inout) :: w 
    double precision, dimension(size(w)), intent(in) :: rhs 
    double precision, intent(in) :: dt, h
    integer, intent(in) :: bc 
    integer, intent(out) :: info 
    double precision, dimension(2), optional, intent(in) :: dwdy, theta 
   
    double precision, dimension(size(w)) :: dd ! Diagonal of tridiagonal matrix 
    double precision, dimension(size(w) - 1) :: dl, du ! Two off-diagonal vectors of tridiag matrix 
    ! Note: these variables must be integer*8 in order for them to be passed to LAPACK functions
    ! in mex-files for 64-bit MATLAB. If they are declared as just "integer", a segmentation fault
    ! may result. 
    integer*8 :: n, info8
    integer*8, parameter :: NRHS = 1 

   ! Add w to rhs 
    w = rhs + w 
    n = size(w) 
   
   ! Create tridiagonal matrix A 
    du = -dt/(2.*h**2)  
    dd = dt/h**2 + 1. 
    dl = -dt/(2.*h**2)   
   
    if (bc == BC_NEUMANN .or. bc == BC_OBLIQUE) then 
       du(1) = -dt / h**2  
       dl(n-1) = -dt / h**2 
    endif 
    if (bc == BC_OBLIQUE) then 
       w(1) = w(1) + dt/(2.*h**2) * dwdy(1) * tan(theta(1))   
       w(n) = w(n) + dt/(2.*h**2) * dwdy(2) * tan(theta(2)) 
    endif 
  
   ! Solve the tridiagonal system  
    call dgtsv (n, NRHS, dl, dd, du, w, n, info8) 

    ! Convert info8 to info
    info8 = min(info8, huge(info))
    info = max(info8, -huge(info))
   
end subroutine tridiagSolveNoD  


! -------------------------------------------------------------------------------------------------
! Implements tridiagSolve when the optional argument d is present. 
! -------------------------------------------------------------------------------------------------
subroutine tridiagSolveD (w, rhs, dt, h, bc, info, d, dprime, dwdy, theta) 
    implicit none 
    double precision, dimension(:), intent(inout) :: w 
    double precision, dimension(size(w)), intent(in) :: rhs 
    double precision, intent(in) :: dt, h
    integer, intent(in) :: bc 
    integer, intent(out) :: info 
    double precision, dimension(size(w)), intent(in) :: d, dprime 
    double precision, dimension(2), optional, intent(in) :: dwdy, theta 
   
    double precision, dimension(size(w)) :: dd ! Diagonal of tridiagonal matrix 
    double precision, dimension(size(w) - 1) :: dl, du ! Two off-diagonal vectors of tridiag matrix 
    ! Note: these variables must be integer*8 in order for them to be passed to LAPACK functions
    ! in mex-files for 64-bit MATLAB. If they are declared as just "integer", a segmentation fault
    ! may result. 
    integer*8 :: n, info8
    integer*8, parameter :: NRHS = 1 
   
   ! Add w to rhs 
    w = rhs + w 
    n = size(w) 
   
   ! Create tridiagonal matrix A 
    du = -dt/(2.*h**2) * d(1:(n-1)) - dt/(4.*h) * dprime(1:(n-1))
    dd = dt/h**2 * d + 1. 
    dl = -dt/(2.*h**2) * d(2:n) + dt/(4.*h) * dprime(2:n) 
   
    if (bc == BC_NEUMANN .or. bc == BC_OBLIQUE) then 
       du(1) = -dt/h**2 * d(1) 
       dl(n-1) = -dt/h**2 * d(n) 
    endif 
    if (bc == BC_OBLIQUE) then 
       w(1) = w(1) + dt/(2.*h**2) * d(1) * dwdy(1) * tan(theta(1)) 
       w(n) = w(n) + dt/(2.*h**2) * d(n) * dwdy(2) * tan(theta(2)) 
    endif 
   
   ! Solve the tridiagonal system  
    call dgtsv (n, NRHS, dl, dd, du, w, n, info8)

    ! Convert info8 to info
    info8 = min(info8, huge(info))
    info = max(info8, -huge(info))
   
end subroutine tridiagSolveD  



! -------------------------------------------------------------------------------------------------
! Checks input arguments. Returns zero if all are correct, or -i if the ith argument is invalid. 
! The order "i" refers to the order args are input into "update" and "FDupdate.mex*". 
! Note that checks on the *size* of parameters are done in MATLAB, in order to avoid memory errors
! that may occur while they are being passed in through mxFDupdate. 
! -------------------------------------------------------------------------------------------------
function checkArgs (C, nt, dt, h, domain, d, gamma)
    implicit none 
    
    integer :: checkArgs, info 
    double precision, dimension(:,:,:), intent(in) :: C 
    double precision, intent(in) :: dt, h, gamma 
    double precision, dimension(size(C,1), size(C,2)), intent(in) :: d 
    integer, intent(in) :: nt
    type(domain_params), intent(in) :: domain 
   
    double precision, parameter :: PI = 4. * atan(1.) 
   
    if (nt < 1) then 
       info = -NT_INDEX
    elseif (dt <= 0.) then 
       info = -DT_INDEX
    elseif (h <= 0.) then 
       info = -H_INDEX
    elseif (any(d < 0.)) then 
       info = -D_INDEX
    elseif (gamma < 0.) then 
       info = -GAMMA_INDEX
    elseif (any(domain%theta < 0.) .or. any(domain%theta >= PI / 2.)) then 
       info = -THETA_INDEX 
    else 
       info = 0
    endif 
   
    checkArgs = info 
    return    
end function checkArgs 


! -------------------------------------------------------------------------------------------------
! "Parallel max." Vectorized version of max function. 
! -------------------------------------------------------------------------------------------------
elemental function pmax (x, y) 
    implicit none 
    double precision :: pmax 
    double precision, intent(in) :: x, y 
    
    pmax = max(x, y)
    return 
end function pmax 

   
end module FDupdate 
      
   
! NOTES 
! 1. It is necessary to record gval as a variable; otherwise, we would be computing it based on the 
!        newly updated value of U but the old value of V.  
! 2. In the event that ymax(1) /= ymax(end), we need to apply Robin B.C.'s only up to the minimum, 
!        to avoid risking calling crobin(ymax(end)) if length(crobin) = ymax(1) < ymax(end). This also
!        has the advantage that we don't have to worry about the x-bounds of points up to ysplit. 
! 3. Note the irregular (two-based, rather than one-based) indices: these make indexing more
!        intuitive, since the entry "iy" corresponds to row iy.
! 4. "max" is used here for a couple reasons. First, if there have been no errors so far, info = 0. 
!        (If there are errors in the arguments to "update", then the routine will already have 
!        exited.) Second, it means we only catch errors from dgtsv that result in info > 0. Errors with
!        info < 0 mean there was an illegal argument to dgtsv, which shouldn't happen anyway in the 
!        code - but just to be safe, this avoids causing confusion between errors from dgtsv and errors
!        with illegal arguments to "update." 





