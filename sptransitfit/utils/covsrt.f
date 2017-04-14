CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
         do 11 j=1,i
            covar(i,j)=0.
            covar(j,i)=0.
 11      enddo
 12   enddo
      k=mfit
      do 15 j=ma,1,-1
         if(ia(j).ne.0)then
            do 13 i=1,ma
               swap=covar(i,k)
               covar(i,k)=covar(i,j)
               covar(i,j)=swap
 13         enddo
            do 14 i=1,ma
               swap=covar(k,i)
               covar(k,i)=covar(j,i)
               covar(j,i)=swap
 14         enddo
            k=k-1
         endif
 15   enddo
      return
      END
