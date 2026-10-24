      program main

      implicit none

      include "maxdim.h"

      integer i,j,N3TM
      character*2 symb(maxatm)
      double precision X(max3),V,DX(max3)

      N3TM = max3

      do i = 1, maxatm
        read (98,*,end=21) symb(i), (X(3*(i-1)+j),j=1,3)
      enddo

 21   continue

      call surf(V,X,DX,N3TM)

      end
