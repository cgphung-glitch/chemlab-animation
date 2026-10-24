c
c Arrange all the valence electrons to the necessary molecular orbitals 
c so that the potential energy of the system is minimized.
c
      subroutine filloc(eva,etb,occ,ntotal,ndim)
c
      include "maxdim.h"
c
      implicit none
      integer ndim,ntotal,i,j,nbottom,ntop
      double precision eva(max9)
      double precision occ(max9),coef(max9),total,u,etb
c
      u=0.d0
c
c Initialize occupation number set
c
      do 10 i=1,ndim
        occ(i)=0.d0
        coef(i)=0.d0
   10 continue
c
      etb=0.d0
c
      if (ntotal .gt. 2*ndim) stop "There is not enough MOs"
c
c Arrange electrons simply from lowest orbitals to highest orbitals
c
      nbottom=ntotal/2
      if (ntotal .eq. 1) then
        occ(1)=1.d0
        etb=eva(1)
        return
      else
        do 20 i=1, ntotal/2
          occ(i)=2.d0
   20   continue
        if (2*(ntotal/2) .ne. ntotal) then
          occ(i)=1.d0
          ntop=2
        else
          ntop=1
        endif
      endif
c
c Rearrange the occupation according to penalty energy.
c
c216  if (eva(nbottom+ntop) .le. (eva(nbottom)+u+0.00001)) then
c       if (nbottom .eq. 0) goto 116
c       if ((nbottom+ntop) .gt. ndim) goto 116
c       occ(nbottom)=occ(nbottom)-1.d0
c       occ(nbottom+ntop)=occ(nbottom+ntop)+1.d0
c       nbottom=nbottom-1
c       ntop=ntop+2
c       goto 216
c     endif
c
c Calculate the total valent energy.
c
 116  total=0.d0
      do 30 i=1,ndim
        total=total+occ(i)
        if (abs(occ(i)-2.d0) .lt. 1.e-05) then
          etb=etb+occ(i)*eva(i)+u
        else
          etb=etb+occ(i)*eva(i)
        endif
 30   continue
      if (abs(ntotal-total).gt.1.e-5) stop "Occupatations is wrong."

      return
      end
