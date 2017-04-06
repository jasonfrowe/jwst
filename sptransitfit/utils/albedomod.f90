!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real(double) function albedomod(Pi,ag,phi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) Pi,phi,alpha,phase,ag

phi=phi+Pi
if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


alpha=abs(phi)
alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere

albedomod=ag*phase

return
end
