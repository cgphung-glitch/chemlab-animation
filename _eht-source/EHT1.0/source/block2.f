      BLOCK DATA set2

      implicit none

      integer maxatm,i
      parameter (maxatm=60)
      logical first
      double precision fact(0:12)
      integer idtype(maxatm),natom,ndim,nElect
      
      common /factorial/ fact
      common /firsRD/ natom,ndim,nElect,idtype,first

      save /firsRD/

      data fact/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
     *          362880.0,3628800.0,39916800.0,479001600.0/

      data natom,ndim,nElect /0,0,0/
      
      data first /.true./

      data (idtype(i),i=1,maxatm) /maxatm * 0/

      end
