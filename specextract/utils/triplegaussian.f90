function triplegaussian(nfit,sol,x)
!Jason Rowe 2015 - jasonfrowe@gmail.com
!returns the sum of three Gaussians and an offset
use precision
implicit none
integer :: nfit,i,np
real(double) :: triplegaussian,x
real(double), dimension(nfit) :: sol

triplegaussian=sol(1)
!first Gaussian
triplegaussian=triplegaussian+                                          &
               sol(2)*sol(8)*exp(-(x-sol(3)-sol(9))**2.0d0/             &
               (2.0d0*sol(4)*sol(4)))
!second Gaussian
triplegaussian=triplegaussian+                                          &
               sol(5)*sol(8)*exp(-(x-sol(6)-sol(9))**2.0d0/             &
               (2.0d0*sol(7)*sol(7)))
!third Gaussian
triplegaussian=triplegaussian+                                          &
               sol(8)*exp(-(x-sol(9))**2.0d0/             &
               (2.0d0*sol(10)*sol(10)))


return
end function triplegaussian
