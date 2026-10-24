C This subroutine is an interface to call two subroutines to compute the 
C overlap integrals (one has the same zeta values, the other has different
C zeta values.
      subroutine ovlcon(Np,Lp,zetap,N,L,M,zeta,intgr,dintgr)
c
      implicit none
c
      integer Np,Lp,N,L,M
c
      double precision zetap,zeta,intgr,dintgr
c
      if (Np.le.Lp .or. N.le.L) then
         intgr=0.0 
        dintgr=0.0 
        return
      endif
c
      if (abs(zetap-zeta) .lt. 0.001) then
        zeta=(zeta+zetap)/2.0
c
        call ovlapd(Np,Lp,zeta,N,L,M,intgr,dintgr)
      else
c
        call dovlap(Np,Lp,zetap,N,L,M,zeta,intgr,dintgr)
      endif
c
      return
      end
