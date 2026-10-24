c This fragment calculate directional cosines between atom i and atom j in 
c the Hamiltonian matrix.

      subroutine lmndlm(i,j,ratom)

      implicit none
      integer i,j,n3tm
      parameter (n3tm=180)
      double precision ratom(3*n3tm)
      double precision dx,dy,dz,r
      double precision l,m,n,lladmm,llmimm,nn
      double precision dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz

      common /lmn/ l,m,n,lladmm,llmimm,nn
      common /dlmn/ dlx,dly,dlz,dmx,dmy,dmz,dnx,dny,dnz
      common /distan/ r

      dx=ratom(3*(i-1)+1)-ratom(3*(j-1)+1)
      dy=ratom(3*(i-1)+2)-ratom(3*(j-1)+2)
      dz=ratom(3*(i-1)+3)-ratom(3*(j-1)+3)
      
      r=dsqrt(dx*dx+dy*dy+dz*dz)

      l=dx/r
      m=dy/r
      n=dz/r

      lladmm=l*l+m*m
      llmimm=l*l-m*m
      nn=n*n

      dlx=(1.0-l*l)/r
      dly=-l*m/r
      dlz=-l*n/r
      dmx=dly
      dmy=(1.0-m*m)/r
      dmz=-m*n/r
      dnx=dlz
      dny=dmz
      dnz=(1.0-n*n)/r

      return
      end
