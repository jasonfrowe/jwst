CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION erf(x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      REAL*8 x
C USES gammp
C Returns the error function erf(x).
      REAL*8 gammp
      if(x.lt.0.)then
         erf=-gammp(.5,x**2)
      else
         erf=gammp(.5,x**2)
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION gammp(a,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 a,x
C USES gcf,gser
C Returns the incomplete gamma function P(a, x).
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause "bad arguments in gammp"
      if(x.lt.a+1.)then !Use the series representation.
         call gser(gamser,a,x,gln)
         gammp=gamser
      else !Use the continued fraction representation
         call gcf(gammcf,a,x,gln)
         gammp=1.-gammcf !and take its complement.
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gser(gamser,a,x,gln)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
C USES gammln
C Returns the incomplete gamma function P(a, x) evaluated by its series representation as
C gamser. Also returns lnΓ(a) as gln.
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
         if(x.lt.0.)pause "x < 0 in gser"
         gamser=0.
         return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*EPS)goto 1
 11   continue
      pause "a too large, ITMAX too small in gser"
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gcf(gammcf,a,x,gln)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
C USES gammln
C Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
C as gammcf. Also returns ln Γ(a) as gln.
C Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accuracy;
C FPMIN is a number near the smallest representable floating-point number.
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a !Set up for evaluating continued fraction by modified
C               Lentz’s method (§5.2) with b0 = 0.
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX !Iterate to convergence.
         an=-i*(i-a)
         b=b+2.
         d=an*d+b
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1./d
         del=d*c
         h=h*del
         if(abs(del-1.).lt.EPS)goto 1
 11   continue
      pause "a too large, ITMAX too small in gcf"
 1    gammcf=exp(-x+a*log(x)-gln)*h !Put factors in front.
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION gammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     From Numerical Recipes
      double precision xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END

