c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between Ni and C.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine exhukl(i,j,hamil,dhamil,S)
c
      implicit none
c
c i,j are atomic numbers of the two atoms.
c sNp,pNp,dNp,sN,pN,dN are principal quantum numbers of atom i (s, p, and d)
c and atom j (s, p, and d).
c
      integer i,j
      integer sNp,pNp,dNp,sN,pN,dN                     !quantumn numbers
c
c Provides an submatrices for interactions between the two atoms.
c  hamil --- Hamiltonian matrix between atom i and atom j.
c dhamil --- First derivatives of hmail.
c      S --- Overlap integral matrix between atom i and atom j.
c The basis functions are in order:
c s, px, py, pz, dxy, dyz, dzx, dxx-yy, and dzz.
c
      double precision hamil(9,9),S(9,9),dhamil(9,9,3),t,dt(3)
c
c The interactions between s of atom i and orbital s, p, and d of atom j.
c
      double precision  sss, sps, sds
      double precision dsss,dsps,dsds
c
c The interactions between p of atom i and orbital s, p, and d of atom j.
c
      double precision  pss, pps, ppp, pds, pdp
      double precision dpss,dpps,dppp,dpds,dpdp
c
c The interactions between d of atom i and orbital s, p, and d of atom j.
c
      double precision   dss, dps, dpp, dds, ddp, ddd
      double precision  ddss,ddps,ddpp,ddds,dddp,dddd
c
c Temporary variables.
c
      double precision zetap,zeta,intgr,dintgr    
c
c STO exponents: single for s and p, double for d.
c
      double precision expss(103),expps(103),expdd1(103),expdd2(103)
c
c Temporary variables.
c
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
c
c nsd --- The principal quantum number of s orbital
c npd --- The principal quantum number of p orbital
c ndd --- The principal quantum number of d orbital
c
      integer          nsd(103),npd(103),ndd(103)
      double precision vsip(103,3),cdd1(103),cdd2(103)
c
c Common block "nspdd" will share with subroutine "ehuckl,surf,data block"
c Common block "coulom" will share with subroutine "ehuckl,data block"
c Common block "dcoeff" will share with subroutine "ehuckl,data block"
c Common block "dtzeta" will share with subroutine "ehuckl,data block"
c
      common /nspdd/ nsd,npd,ndd
      common /coulom/ vsip
      common /dcoeff/ cdd1,cdd2
      common /dtzeta/ expss,expps,expdd1,expdd2
c
c STO parameters for the atom with atomic number i
c
      sNp=nsd(i)
      pNp=npd(i)
      dNp=ndd(i)
c
c STO parameters for the atom with atomic number j
c
       sN=nsd(j)
       pN=npd(j)
       dN=ndd(j)
c
c Wolfberg-Holmholz constants
c
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
c s s sigma bonding overlap integral in z-axis
      zetap=expss(i)
       zeta=expss(j)
      call ovlcon(sNp,0,zetap,sN,0,0,zeta,intgr,dintgr)
       sss=intgr
      dsss=dintgr

c s p sigma bonding overlap integral in z-axis
      if (pN.eq.0) then
         sps=0.0
        dsps=0.0
      else
        zetap=expss(i)
         zeta=expps(j)
        call ovlcon(sNp,0,zetap,pN,1,0,zeta,intgr,dintgr)
         sps=intgr
        dsps=dintgr
      endif

c s d sigma bonding overlap integral in z-axis
      if (dN.eq.0) then
         sds=0.0
        dsds=0.0
      else
        zetap=expss(i)
         zeta=expdd1(j)
        call ovlcon(sNp,0,zetap,dN,2,0,zeta,intgr,dintgr)
         sds= intgr*cdd1(j)
        dsds=dintgr*cdd1(j)

        zetap=expss(i)
         zeta=expdd2(j)
        call ovlcon(sNp,0,zetap,dN,2,0,zeta,intgr,dintgr)
         sds= intgr*cdd2(j)+ sds
        dsds=dintgr*cdd2(j)+dsds
      endif

c p s sigma bonding overlap integral in z-axis
      if (pNp.eq.0) then
         pss=0.0
        dpss=0.0
      else
        zetap=expps(i)
         zeta=expss(j)
        call ovlcon(pNp,1,zetap,sN,0,0,zeta,intgr,dintgr)
         pss=intgr
        dpss=dintgr
      endif

c p p sigma bonding overlap integral in z-axis
      if (pNp.eq.0 .or. pN.eq.0) then
         pps=0.0
        dpps=0.0
      else
        zetap=expps(i)
         zeta=expps(j)
        call ovlcon(pNp,1,zetap,pN,1,0,zeta,intgr,dintgr)
         pps=intgr
        dpps=dintgr
      endif

c p p pi bonding overlap integral in z-axis
      if (pNp.eq.0 .or. pN.eq.0) then
         ppp=0.0
        dppp=0.0
      else
        zetap=expps(i)
         zeta=expps(j)
        call ovlcon(pNp,1,zetap,pN,1,1,zeta,intgr,dintgr)
         ppp=intgr
        dppp=dintgr
      endif

c p d sigma bonding overlap integral in z-axis
      if (pNp.eq.0 .or. dN.eq.0) then
         pds=0.0
        dpds=0.0
      else
        zetap=expps(i)
         zeta=expdd1(j)
        call ovlcon(pNp,1,zetap,dN,2,0,zeta,intgr,dintgr)
         pds= intgr*cdd1(j)
        dpds=dintgr*cdd1(j)

        zetap=expps(i)
         zeta=expdd2(j)
        call ovlcon(pNp,1,zetap,dN,2,0,zeta,intgr,dintgr)
         pds= intgr*cdd2(j)+ pds
        dpds=dintgr*cdd2(j)+dpds
      endif

c p d pi bonding overlap integral in z-axis
      if (pNp.eq.0 .or. dN.eq.0) then
         pdp=0.0
        dpdp=0.0
      else
        zetap=expps(i)
         zeta=expdd1(j)
        call ovlcon(pNp,1,zetap,dN,2,1,zeta,intgr,dintgr)
         pdp= intgr*cdd1(j)
        dpdp=dintgr*cdd1(j)

        zetap=expps(i)
         zeta=expdd2(j)
        call ovlcon(pNp,1,zetap,dN,2,1,zeta,intgr,dintgr)
         pdp= intgr*cdd2(j)+ pdp
        dpdp=dintgr*cdd2(j)+dpdp
      endif

c d s sigma bonding overlap integral in z-axis
      if (dNp.eq.0) then
         dss=0.0
        ddss=0.0
      else
        zetap=expdd1(i)
         zeta=expss(j)
        call ovlcon(dNp,2,zetap,sN,0,0,zeta,intgr,dintgr)
         dss= intgr*cdd1(i)
        ddss=dintgr*cdd1(i)

        zetap=expdd2(i)
         zeta=expss(j)
        call ovlcon(dNp,2,zetap,sN,0,0,zeta,intgr,dintgr)
         dss= intgr*cdd2(i)+ dss
        ddss=dintgr*cdd2(i)+ddss
      endif

c d p sigma bonding overlap integral in z-axis
      if (dNp.eq.0 .or. pN.eq.0) then
         dps=0.0
        ddps=0.0
      else
        zetap=expdd1(i)
         zeta=expps(j)
        call ovlcon(dNp,2,zetap,pN,1,0,zeta,intgr,dintgr)
         dps= intgr*cdd1(i)
        ddps=dintgr*cdd1(i)

        zetap=expdd2(i)
         zeta=expps(j)
        call ovlcon(dNp,2,zetap,pN,1,0,zeta,intgr,dintgr)
         dps= intgr*cdd2(i)+ dps
        ddps=dintgr*cdd2(i)+ddps
      endif

c d p pi bonding overlap integral in z-axis
      if (dNp.eq.0 .or. pN.eq.0) then
         dps=0.0
        ddps=0.0
      else
        zetap=expdd1(i)
         zeta=expps(j)
        call ovlcon(dNp,2,zetap,pN,1,1,zeta,intgr,dintgr)
         dpp= intgr*cdd1(i)
        ddpp=dintgr*cdd1(i)

        zetap=expdd2(i)
         zeta=expps(j)
        call ovlcon(dNp,2,zetap,pN,1,1,zeta,intgr,dintgr)
         dpp= intgr*cdd2(i)+ dpp
        ddpp=dintgr*cdd2(i)+ddpp
      endif

c d d sigma bonding overlap integral in z-axis
      if (dNp.eq.0 .or. dN.eq.0) then
         dds=0.0
        ddds=0.0
      else
        zetap=expdd1(i)
       zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,0,zeta,intgr,dintgr)
         dds= intgr*cdd1(i)*cdd1(j)
        ddds=dintgr*cdd1(i)*cdd1(j)

        zetap=expdd1(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,0,zeta,intgr,dintgr)
         dds= intgr*cdd1(i)*cdd2(j)+ dds
        ddds=dintgr*cdd1(i)*cdd2(j)+ddds

        zetap=expdd2(i)
         zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,0,zeta,intgr,dintgr)
         dds= intgr*cdd2(i)*cdd1(j)+ dds
        ddds=dintgr*cdd2(i)*cdd1(j)+ddds

        zetap=expdd2(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,0,zeta,intgr,dintgr)
         dds= intgr*cdd2(i)*cdd2(j)+ dds
        ddds=dintgr*cdd2(i)*cdd2(j)+ddds
      endif

c d d pi bonding overlap integral in z-axis
      if (dNp.eq.0 .or. dN.eq.0) then
         ddp=0.0
        dddp=0.0
      else
        zetap=expdd1(i)
         zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,1,zeta,intgr,dintgr)
         ddp= intgr*cdd1(i)*cdd1(j)
        dddp=dintgr*cdd1(i)*cdd1(j)

        zetap=expdd1(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,1,zeta,intgr,dintgr)
         ddp= intgr*cdd1(i)*cdd2(j)+ ddp
        dddp=dintgr*cdd1(i)*cdd2(j)+dddp

        zetap=expdd2(i)
         zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,1,zeta,intgr,dintgr)
         ddp= intgr*cdd2(i)*cdd1(j)+ ddp
        dddp=dintgr*cdd2(i)*cdd1(j)+dddp

        zetap=expdd2(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,1,zeta,intgr,dintgr)
         ddp= intgr*cdd2(i)*cdd2(j)+ ddp
        dddp=dintgr*cdd2(i)*cdd2(j)+dddp
      endif

c d d delta bonding overlap integral in z-axis
      if (dNp.eq.0 .or. pN.eq.0) then
         ddd=0.0
        dddd=0.0
      else
        zetap=expdd1(i)
         zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,2,zeta,intgr,dintgr)
         ddd= intgr*cdd1(i)*cdd1(j)
        dddd=dintgr*cdd1(i)*cdd1(j)

        zetap=expdd1(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,2,zeta,intgr,dintgr)
         ddd= intgr*cdd1(i)*cdd2(j)+ ddd
        dddd=dintgr*cdd1(i)*cdd2(j)+dddd

        zetap=expdd2(i)
         zeta=expdd1(j)
        call ovlcon(dNp,2,zetap,dN,2,2,zeta,intgr,dintgr)
         ddd= intgr*cdd2(i)*cdd1(j)+ ddd
        dddd=dintgr*cdd2(i)*cdd1(j)+dddd

        zetap=expdd2(i)
         zeta=expdd2(j)
        call ovlcon(dNp,2,zetap,dN,2,2,zeta,intgr,dintgr)
         ddd= intgr*cdd2(i)*cdd2(j)+ ddd
        dddd=dintgr*cdd2(i)*cdd2(j)+dddd
      endif
c
c Multiply kmn with (Hm(Ni)+Hn(C))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=1.75*(vsip(i,1)+vsip(j,1))/2.0
      tmp2=1.75*(vsip(i,1)+vsip(j,2))/2.0
      tmp3=1.75*(vsip(i,1)+vsip(j,3))/2.0
      tmp4=1.75*(vsip(i,2)+vsip(j,1))/2.0
      tmp5=1.75*(vsip(i,2)+vsip(j,2))/2.0
      tmp6=1.75*(vsip(i,2)+vsip(j,3))/2.0
      tmp7=1.75*(vsip(i,3)+vsip(j,1))/2.0
      tmp8=1.75*(vsip(i,3)+vsip(j,2))/2.0
      tmp9=1.75*(vsip(i,3)+vsip(j,3))/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      S(1,1)=t
      hamil(1,1)=t*tmp1
      dhamil(1,1,1)=dt(1)*tmp1
      dhamil(1,1,2)=dt(2)*tmp1
      dhamil(1,1,3)=dt(3)*tmp1

      call sxyz('x',sps,dsps,t,dt)       !s-px
      S(1,2)=t
      hamil(1,2)=t*tmp2
      dhamil(1,2,1)=dt(1)*tmp2
      dhamil(1,2,2)=dt(2)*tmp2
      dhamil(1,2,3)=dt(3)*tmp2

      call sxyz('y',sps,dsps,t,dt)       !s-py
      S(1,3)=t
      hamil(1,3)=t*tmp2
      dhamil(1,3,1)=dt(1)*tmp2
      dhamil(1,3,2)=dt(2)*tmp2
      dhamil(1,3,3)=dt(3)*tmp2

      call sxyz('z',sps,dsps,t,dt)       !s-pz
      S(1,4)=t
      hamil(1,4)=t*tmp2
      dhamil(1,4,1)=dt(1)*tmp2
      dhamil(1,4,2)=dt(2)*tmp2
      dhamil(1,4,3)=dt(3)*tmp2

      call sxy('xy',sds,dsds,t,dt)       !s-dxy
      S(1,5)=t
      hamil(1,5)=t*tmp3
      dhamil(1,5,1)=dt(1)*tmp3
      dhamil(1,5,2)=dt(2)*tmp3
      dhamil(1,5,3)=dt(3)*tmp3

      call sxy('yz',sds,dsds,t,dt)       !s-dyz
      S(1,6)=t
      hamil(1,6)=t*tmp3
      dhamil(1,6,1)=dt(1)*tmp3
      dhamil(1,6,2)=dt(2)*tmp3
      dhamil(1,6,3)=dt(3)*tmp3

      call sxy('zx',sds,dsds,t,dt)       !s-dzx
      S(1,7)=t
      hamil(1,7)=t*tmp3
      dhamil(1,7,1)=dt(1)*tmp3
      dhamil(1,7,2)=dt(2)*tmp3
      dhamil(1,7,3)=dt(3)*tmp3

      call sxxyy(sds,dsds,t,dt)          !s-dx2_y2
      S(1,8)=t
      hamil(1,8)=t*tmp3
      dhamil(1,8,1)=dt(1)*tmp3
      dhamil(1,8,2)=dt(2)*tmp3
      dhamil(1,8,3)=dt(3)*tmp3

      call szz(sds,dsds,t,dt)            !s-dz2
      S(1,9)=t
      hamil(1,9)=t*tmp3
      dhamil(1,9,1)=dt(1)*tmp3
      dhamil(1,9,2)=dt(2)*tmp3
      dhamil(1,9,3)=dt(3)*tmp3

      call sxyz('x',pss,dpss,t,dt)       !px-s
      S(2,1)=t
      hamil(2,1)=t*tmp4
      dhamil(2,1,1)=dt(1)*tmp4
      dhamil(2,1,2)=dt(2)*tmp4
      dhamil(2,1,3)=dt(3)*tmp4

      call xx('x',pps,dpps,ppp,dppp,t,dt)   !px-px
      S(2,2)=t
      hamil(2,2)=t*tmp5
      dhamil(2,2,1)=dt(1)*tmp5
      dhamil(2,2,2)=dt(2)*tmp5
      dhamil(2,2,3)=dt(3)*tmp5

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !px-py
      S(2,3)=t
      hamil(2,3)=t*tmp5
      dhamil(2,3,1)=dt(1)*tmp5
      dhamil(2,3,2)=dt(2)*tmp5
      dhamil(2,3,3)=dt(3)*tmp5

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !px-pz
      S(2,4)=t
      hamil(2,4)=t*tmp5
      dhamil(2,4,1)=dt(1)*tmp5
      dhamil(2,4,2)=dt(2)*tmp5
      dhamil(2,4,3)=dt(3)*tmp5

      call xxy('xy',pds,dpds,pdp,dpdp,t,dt)         !px-dxy
      S(2,5)=t
      hamil(2,5)=t*tmp6
      dhamil(2,5,1)=dt(1)*tmp6
      dhamil(2,5,2)=dt(2)*tmp6
      dhamil(2,5,3)=dt(3)*tmp6

      call xyz(pds,dpds,pdp,dpdp,t,dt)              !px-dyz
      S(2,6)=t
      hamil(2,6)=t*tmp6
      dhamil(2,6,1)=dt(1)*tmp6
      dhamil(2,6,2)=dt(2)*tmp6
      dhamil(2,6,3)=dt(3)*tmp6

      call xxy('xz',pds,dpds,pdp,dpdp,t,dt)         !px-dzx
      S(2,7)=t
      hamil(2,7)=t*tmp6
      dhamil(2,7,1)=dt(1)*tmp6
      dhamil(2,7,2)=dt(2)*tmp6
      dhamil(2,7,3)=dt(3)*tmp6

      call xxxyy(pds,dpds,pdp,dpdp,t,dt)            !px-dx2_y2
      S(2,8)=t
      hamil(2,8)=t*tmp6
      dhamil(2,8,1)=dt(1)*tmp6
      dhamil(2,8,2)=dt(2)*tmp6
      dhamil(2,8,3)=dt(3)*tmp6

      call xzz('x',pds,dpds,pdp,dpdp,t,dt)              !px-dz2
      S(2,9)=t
      hamil(2,9)=t*tmp6
      dhamil(2,9,1)=dt(1)*tmp6
      dhamil(2,9,2)=dt(2)*tmp6
      dhamil(2,9,3)=dt(3)*tmp6

      call sxyz('y',pss,dpss,t,dt)          !py-s
      S(3,1)=t
      hamil(3,1)=t*tmp4
      dhamil(3,1,1)=dt(1)*tmp4
      dhamil(3,1,2)=dt(2)*tmp4
      dhamil(3,1,3)=dt(3)*tmp4

      call xy('xy',pps,dpps,ppp,dppp,t,dt)  !py-px
      S(3,2)=t
      hamil(3,2)=t*tmp5
      dhamil(3,2,1)=dt(1)*tmp5
      dhamil(3,2,2)=dt(2)*tmp5
      dhamil(3,2,3)=dt(3)*tmp5

      call xx('y',pps,dpps,ppp,dppp,t,dt)   !py-py
      S(3,3)=t
      hamil(3,3)=t*tmp5
      dhamil(3,3,1)=dt(1)*tmp5
      dhamil(3,3,2)=dt(2)*tmp5
      dhamil(3,3,3)=dt(3)*tmp5

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !py-pz
      S(3,4)=t
      hamil(3,4)=t*tmp5
      dhamil(3,4,1)=dt(1)*tmp5
      dhamil(3,4,2)=dt(2)*tmp5
      dhamil(3,4,3)=dt(3)*tmp5

      call xxy('yx',pds,dpds,pdp,dpdp,t,dt)         !py-dxy
      S(3,5)=t
      hamil(3,5)=t*tmp6
      dhamil(3,5,1)=dt(1)*tmp6
      dhamil(3,5,2)=dt(2)*tmp6
      dhamil(3,5,3)=dt(3)*tmp6

      call xxy('yz',pds,dpds,pdp,dpdp,t,dt)         !py-dyz
      S(3,6)=t
      hamil(3,6)=t*tmp6
      dhamil(3,6,1)=dt(1)*tmp6
      dhamil(3,6,2)=dt(2)*tmp6
      dhamil(3,6,3)=dt(3)*tmp6

      call xyz(pds,dpds,pdp,dpdp,t,dt)              !py-dzx
      S(3,7)=t
      hamil(3,7)=t*tmp6
      dhamil(3,7,1)=dt(1)*tmp6
      dhamil(3,7,2)=dt(2)*tmp6
      dhamil(3,7,3)=dt(3)*tmp6

      call yxxyy(pds,dpds,pdp,dpdp,t,dt)            !py-dx2_y2
      S(3,8)=t
      hamil(3,8)=t*tmp6
      dhamil(3,8,1)=dt(1)*tmp6
      dhamil(3,8,2)=dt(2)*tmp6
      dhamil(3,8,3)=dt(3)*tmp6

      call xzz('y',pds,dpds,pdp,dpdp,t,dt)          !py-dz2
      S(3,9)=t
      hamil(3,9)=t*tmp6
      dhamil(3,9,1)=dt(1)*tmp6
      dhamil(3,9,2)=dt(2)*tmp6
      dhamil(3,9,3)=dt(3)*tmp6

      call sxyz('z',pss,dpss,t,dt)          !pz-s
      S(4,1)=t
      hamil(4,1)=t*tmp4
      dhamil(4,1,1)=dt(1)*tmp4
      dhamil(4,1,2)=dt(2)*tmp4
      dhamil(4,1,3)=dt(3)*tmp4

      call xy('xz',pps,dpps,ppp,dppp,t,dt)  !pz-px
      S(4,2)=t
      hamil(4,2)=t*tmp5
      dhamil(4,2,1)=dt(1)*tmp5
      dhamil(4,2,2)=dt(2)*tmp5
      dhamil(4,2,3)=dt(3)*tmp5

      call xy('yz',pps,dpps,ppp,dppp,t,dt)  !pz-py
      S(4,3)=t
      hamil(4,3)=t*tmp5
      dhamil(4,3,1)=dt(1)*tmp5
      dhamil(4,3,2)=dt(2)*tmp5
      dhamil(4,3,3)=dt(3)*tmp5

      call xx('z',pps,dpps,ppp,dppp,t,dt)   !pz-pz
      S(4,4)=t
      hamil(4,4)=t*tmp5
      dhamil(4,4,1)=dt(1)*tmp5
      dhamil(4,4,2)=dt(2)*tmp5
      dhamil(4,4,3)=dt(3)*tmp5

      call xyz(pds,dpds,pdp,dpdp,t,dt)               !pz-dxy
      S(4,5)=t
      hamil(4,5)=t*tmp6
      dhamil(4,5,1)=dt(1)*tmp6
      dhamil(4,5,2)=dt(2)*tmp6
      dhamil(4,5,3)=dt(3)*tmp6

      call xxy('zy',pds,dpds,pdp,dpdp,t,dt)          !pz-dyz
      S(4,6)=t
      hamil(4,6)=t*tmp6
      dhamil(4,6,1)=dt(1)*tmp6
      dhamil(4,6,2)=dt(2)*tmp6
      dhamil(4,6,3)=dt(3)*tmp6

      call xxy('zx',pds,dpds,pdp,dpdp,t,dt)          !pz-dzx
      S(4,7)=t
      hamil(4,7)=t*tmp6
      dhamil(4,7,1)=dt(1)*tmp6
      dhamil(4,7,2)=dt(2)*tmp6
      dhamil(4,7,3)=dt(3)*tmp6

      call zxxyy(pds,dpds,pdp,dpdp,t,dt)             !pz-dx2_y2
      S(4,8)=t
      hamil(4,8)=t*tmp6
      dhamil(4,8,1)=dt(1)*tmp6
      dhamil(4,8,2)=dt(2)*tmp6
      dhamil(4,8,3)=dt(3)*tmp6

      call zzz(pds,dpds,pdp,dpdp,t,dt)               !pz-dz2
      S(4,9)=t
      hamil(4,9)=t*tmp6
      dhamil(4,9,1)=dt(1)*tmp6
      dhamil(4,9,2)=dt(2)*tmp6
      dhamil(4,9,3)=dt(3)*tmp6

      call sxy('xy',dss,ddss,t,dt)                  !dxy-s
      S(5,1)=t
      hamil(5,1)=t*tmp7
      dhamil(5,1,1)=dt(1)*tmp7
      dhamil(5,1,2)=dt(2)*tmp7
      dhamil(5,1,3)=dt(3)*tmp7

      call xxy('xy',dps,ddps,dpp,ddpp,t,dt)         !dxy-px
      S(5,2)=t
      hamil(5,2)=t*tmp8
      dhamil(5,2,1)=dt(1)*tmp8
      dhamil(5,2,2)=dt(2)*tmp8
      dhamil(5,2,3)=dt(3)*tmp8

      call xxy('yx',dps,ddps,dpp,ddpp,t,dt)         !dxy-py
      S(5,3)=t
      hamil(5,3)=t*tmp8
      dhamil(5,3,1)=dt(1)*tmp8
      dhamil(5,3,2)=dt(2)*tmp8
      dhamil(5,3,3)=dt(3)*tmp8

      call xyz(dps,ddps,dpp,ddpp,t,dt)              !dxy-pz
      S(5,4)=t
      hamil(5,4)=t*tmp8
      dhamil(5,4,1)=dt(1)*tmp8
      dhamil(5,4,2)=dt(2)*tmp8
      dhamil(5,4,3)=dt(3)*tmp8

      call xyxy('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dxy
      S(5,5)=t
      hamil(5,5)=t*tmp9
      dhamil(5,5,1)=dt(1)*tmp9
      dhamil(5,5,2)=dt(2)*tmp9
      dhamil(5,5,3)=dt(3)*tmp9

      call xyyz('xz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dyz
      S(5,6)=t
      hamil(5,6)=t*tmp9
      dhamil(5,6,1)=dt(1)*tmp9
      dhamil(5,6,2)=dt(2)*tmp9
      dhamil(5,6,3)=dt(3)*tmp9

      call xyyz('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dzx
      S(5,7)=t
      hamil(5,7)=t*tmp9
      dhamil(5,7,1)=dt(1)*tmp9
      dhamil(5,7,2)=dt(2)*tmp9
      dhamil(5,7,3)=dt(3)*tmp9

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dxy-dx2_y2
      S(5,8)=t
      hamil(5,8)=t*tmp9
      dhamil(5,8,1)=dt(1)*tmp9
      dhamil(5,8,2)=dt(2)*tmp9
      dhamil(5,8,3)=dt(3)*tmp9

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      S(5,9)=t
      hamil(5,9)=t*tmp9
      dhamil(5,9,1)=dt(1)*tmp9
      dhamil(5,9,2)=dt(2)*tmp9
      dhamil(5,9,3)=dt(3)*tmp9

      call sxy('yz',dss,ddss,t,dt)                  !dyz-s
      S(6,1)=t
      hamil(6,1)=t*tmp7
      dhamil(6,1,1)=dt(1)*tmp7
      dhamil(6,1,2)=dt(2)*tmp7
      dhamil(6,1,3)=dt(3)*tmp7

      call xyz(dps,ddps,dpp,ddpp,t,dt)              !dyz-px
      S(6,2)=t
      hamil(6,2)=t*tmp8
      dhamil(6,2,1)=dt(1)*tmp8
      dhamil(6,2,2)=dt(2)*tmp8
      dhamil(6,2,3)=dt(3)*tmp8

      call xxy('yz',dps,ddps,dpp,ddpp,t,dt)         !dyz-py
      S(6,3)=t
      hamil(6,3)=t*tmp8
      dhamil(6,3,1)=dt(1)*tmp8
      dhamil(6,3,2)=dt(2)*tmp8
      dhamil(6,3,3)=dt(3)*tmp8

      call xxy('zy',dps,ddps,dpp,ddpp,t,dt)          !dyz-pz
      S(6,4)=t
      hamil(6,4)=t*tmp8
      dhamil(6,4,1)=dt(1)*tmp8
      dhamil(6,4,2)=dt(2)*tmp8
      dhamil(6,4,3)=dt(3)*tmp8

      call xyyz('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dxy
      S(6,5)=t
      hamil(6,5)=t*tmp9
      dhamil(6,5,1)=dt(1)*tmp9
      dhamil(6,5,2)=dt(2)*tmp9
      dhamil(6,5,3)=dt(3)*tmp9

      call xyxy('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dyz
      S(6,6)=t
      hamil(6,6)=t*tmp9
      dhamil(6,6,1)=dt(1)*tmp9
      dhamil(6,6,2)=dt(2)*tmp9
      dhamil(6,6,3)=dt(3)*tmp9

      call xyyz('yx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dzx
      S(6,7)=t
      hamil(6,7)=t*tmp9
      dhamil(6,7,1)=dt(1)*tmp9
      dhamil(6,7,2)=dt(2)*tmp9
      dhamil(6,7,3)=dt(3)*tmp9

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dyz-dx2_y2
      S(6,8)=t
      hamil(6,8)=t*tmp9
      dhamil(6,8,1)=dt(1)*tmp9
      dhamil(6,8,2)=dt(2)*tmp9
      dhamil(6,8,3)=dt(3)*tmp9

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      S(6,9)=t
      hamil(6,9)=t*tmp9
      dhamil(6,9,1)=dt(1)*tmp9
      dhamil(6,9,2)=dt(2)*tmp9
      dhamil(6,9,3)=dt(3)*tmp9

      call sxy('zx',dss,ddss,t,dt)                  !dzx-s
      S(7,1)=t
      hamil(7,1)=t*tmp7
      dhamil(7,1,1)=dt(1)*tmp7
      dhamil(7,1,2)=dt(2)*tmp7
      dhamil(7,1,3)=dt(3)*tmp7

      call xxy('xz',dps,ddps,dpp,ddpp,t,dt)         !dzx-px
      S(7,2)=t
      hamil(7,2)=t*tmp8
      dhamil(7,2,1)=dt(1)*tmp8
      dhamil(7,2,2)=dt(2)*tmp8
      dhamil(7,2,3)=dt(3)*tmp8

      call xyz(dps,ddps,dpp,ddpp,t,dt)              !dzx-py
      S(7,3)=t
      hamil(7,3)=t*tmp8
      dhamil(7,3,1)=dt(1)*tmp8
      dhamil(7,3,2)=dt(2)*tmp8
      dhamil(7,3,3)=dt(3)*tmp8

      call xxy('zx',dps,ddps,dpp,ddpp,t,dt)          !dzx-pz
      S(7,4)=t
      hamil(7,4)=t*tmp8
      dhamil(7,4,1)=dt(1)*tmp8
      dhamil(7,4,2)=dt(2)*tmp8
      dhamil(7,4,3)=dt(3)*tmp8

      call xyyz('zy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dxy
      S(7,5)=t
      hamil(7,5)=t*tmp9
      dhamil(7,5,1)=dt(1)*tmp9
      dhamil(7,5,2)=dt(2)*tmp9
      dhamil(7,5,3)=dt(3)*tmp9

      call xyyz('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dyz
      S(7,6)=t
      hamil(7,6)=t*tmp9
      dhamil(7,6,1)=dt(1)*tmp9
      dhamil(7,6,2)=dt(2)*tmp9
      dhamil(7,6,3)=dt(3)*tmp9

      call xyxy('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dzx
      S(7,7)=t
      hamil(7,7)=t*tmp9
      dhamil(7,7,1)=dt(1)*tmp9
      dhamil(7,7,2)=dt(2)*tmp9
      dhamil(7,7,3)=dt(3)*tmp9

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      S(7,8)=t
      hamil(7,8)=t*tmp9
      dhamil(7,8,1)=dt(1)*tmp9
      dhamil(7,8,2)=dt(2)*tmp9
      dhamil(7,8,3)=dt(3)*tmp9

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dzx-dz2
      S(7,9)=t
      hamil(7,9)=t*tmp9
      dhamil(7,9,1)=dt(1)*tmp9
      dhamil(7,9,2)=dt(2)*tmp9
      dhamil(7,9,3)=dt(3)*tmp9

      call sxxyy(dss,ddss,t,dt)                     !dx2_y2-s
      S(8,1)=t
      hamil(8,1)=t*tmp7
      dhamil(8,1,1)=dt(1)*tmp7
      dhamil(8,1,2)=dt(2)*tmp7
      dhamil(8,1,3)=dt(3)*tmp7

      call xxxyy(dps,ddps,dpp,ddpp,t,dt)            !dx2_y2-px
      S(8,2)=t
      hamil(8,2)=t*tmp8
      dhamil(8,2,1)=dt(1)*tmp8
      dhamil(8,2,2)=dt(2)*tmp8
      dhamil(8,2,3)=dt(3)*tmp8

      call yxxyy(dps,ddps,dpp,ddpp,t,dt)            !dx2_y2-py
      S(8,3)=t
      hamil(8,3)=t*tmp8
      dhamil(8,3,1)=dt(1)*tmp8
      dhamil(8,3,2)=dt(2)*tmp8
      dhamil(8,3,3)=dt(3)*tmp8

      call zxxyy(dps,ddps,dpp,ddpp,t,dt)            !dx2_y2-pz
      S(8,4)=t
      hamil(8,4)=t*tmp8
      dhamil(8,4,1)=dt(1)*tmp8
      dhamil(8,4,2)=dt(2)*tmp8
      dhamil(8,4,3)=dt(3)*tmp8

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      S(8,5)=t
      hamil(8,5)=t*tmp9
      dhamil(8,5,1)=dt(1)*tmp9
      dhamil(8,5,2)=dt(2)*tmp9
      dhamil(8,5,3)=dt(3)*tmp9

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      S(8,6)=t
      hamil(8,6)=t*tmp9
      dhamil(8,6,1)=dt(1)*tmp9
      dhamil(8,6,2)=dt(2)*tmp9
      dhamil(8,6,3)=dt(3)*tmp9

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dzx
      S(8,7)=t
      hamil(8,7)=t*tmp9
      dhamil(8,7,1)=dt(1)*tmp9
      dhamil(8,7,2)=dt(2)*tmp9
      dhamil(8,7,3)=dt(3)*tmp9

      call xxyy_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dx2_y2
      S(8,8)=t
      hamil(8,8)=t*tmp9
      dhamil(8,8,1)=dt(1)*tmp9
      dhamil(8,8,2)=dt(2)*tmp9
      dhamil(8,8,3)=dt(3)*tmp9

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      S(8,9)=t
      hamil(8,9)=t*tmp9
      dhamil(8,9,1)=dt(1)*tmp9
      dhamil(8,9,2)=dt(2)*tmp9
      dhamil(8,9,3)=dt(3)*tmp9

      call szz(dss,ddss,t,dt)            !dz2-s
      S(9,1)=t
      hamil(9,1)=t*tmp7
      dhamil(9,1,1)=dt(1)*tmp7
      dhamil(9,1,2)=dt(2)*tmp7
      dhamil(9,1,3)=dt(3)*tmp7

      call xzz('x',dps,ddps,dpp,ddpp,t,dt)              !dz2-px
      S(9,2)=t
      hamil(9,2)=t*tmp8
      dhamil(9,2,1)=dt(1)*tmp8
      dhamil(9,2,2)=dt(2)*tmp8
      dhamil(9,2,3)=dt(3)*tmp8

      call xzz('y',dps,ddps,dpp,ddpp,t,dt)              !dz2-py
      S(9,3)=t
      hamil(9,3)=t*tmp8
      dhamil(9,3,1)=dt(1)*tmp8
      dhamil(9,3,2)=dt(2)*tmp8
      dhamil(9,3,3)=dt(3)*tmp8

      call zzz(dps,ddps,dpp,ddpp,t,dt)               !dz2-pz
      S(9,4)=t
      hamil(9,4)=t*tmp8
      dhamil(9,4,1)=dt(1)*tmp8
      dhamil(9,4,2)=dt(2)*tmp8
      dhamil(9,4,3)=dt(3)*tmp8

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      S(9,5)=t
      hamil(9,5)=t*tmp9
      dhamil(9,5,1)=dt(1)*tmp9
      dhamil(9,5,2)=dt(2)*tmp9
      dhamil(9,5,3)=dt(3)*tmp9

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      S(9,6)=t
      hamil(9,6)=t*tmp9
      dhamil(9,6,1)=dt(1)*tmp9
      dhamil(9,6,2)=dt(2)*tmp9
      dhamil(9,6,3)=dt(3)*tmp9

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dzx
      S(9,7)=t
      hamil(9,7)=t*tmp9
      dhamil(9,7,1)=dt(1)*tmp9
      dhamil(9,7,2)=dt(2)*tmp9
      dhamil(9,7,3)=dt(3)*tmp9

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dz2-dx2_y2
      S(9,8)=t
      hamil(9,8)=t*tmp9
      dhamil(9,8,1)=dt(1)*tmp9
      dhamil(9,8,2)=dt(2)*tmp9
      dhamil(9,8,3)=dt(3)*tmp9

      call zz_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dz2
      S(9,9)=t
      hamil(9,9)=t*tmp9
      dhamil(9,9,1)=dt(1)*tmp9
      dhamil(9,9,2)=dt(2)*tmp9
      dhamil(9,9,3)=dt(3)*tmp9

      return
      end
