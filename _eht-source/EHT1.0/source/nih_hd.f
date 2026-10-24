c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between C and H.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine nih_hd (hamil,dhamil)
c
      implicit none
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      integer sNp,sLp,dNp,dLp,sN,sL,ssM,dsM
      double precision szeta,szetapm,sss,dsss,dss,ddss,vsip(3,6)
      double precision coef1,coef2,dzeta1,dzeta2,dss1,ddss1,dss2,ddss2
      double precision tmp1,tmp2,kss,kds
      common /coulom/ vsip
c
c STO parameters for nickel element
c
      DATA sNp,sLp,dNp,dLp,ssM,dsM /4,0,3,2,0,0/
      DATA szetapm,coef1,coef2,dzeta1,dzeta2 /1.825,0.5683,0.6292,
     *         5.75,2.0/
c
c STO parameters for hydrogen element
c
      DATA sN,sL /1,0/
      DATA szeta /1.3d0/
c
c Wolfberg-Hol relation k constant.
c
      DATA kss,kds /1.2458,1.0782/
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      call dovlap(sNp,sLp,szetapm,sN,sL,ssM,szeta,sss ,dsss)
      call dovlap(dNp,dLp,dzeta1, sN,sL,dsM,szeta,dss1,ddss1)
      call dovlap(dNp,dLp,dzeta2, sN,sL,dsM,szeta,dss2,ddss2)
c 
      dss=coef1*dss1+coef2*dss2
      ddss=coef1*ddss1+coef2*ddss2
c
c Multiply kmn with (Hm(Ni)+Hn(H))/2 to make a Wolfsberg-Helmholtz relation.
c
      tmp1=(vsip(3,1)+vsip(1,1))/2.0
      tmp2=(vsip(3,2)+vsip(1,1))/2.0
c
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)                          !s-s
      hamil(1,1)=t*kss*tmp1
      dhamil(1,1,1)=dt(1)*kss*tmp1
      dhamil(1,1,2)=dt(2)*kss*tmp1
      dhamil(1,1,3)=dt(3)*kss*tmp1

      call sxy('xy',dss,ddss,t,dt)                    !dxy-s
      hamil(2,1)=t*kds*tmp2
      dhamil(2,1,1)=dt(1)*kds*tmp2
      dhamil(2,1,2)=dt(2)*kds*tmp2
      dhamil(2,1,3)=dt(3)*kds*tmp2

      call sxy('yz',dss,ddss,t,dt)                    !dyz-s
      hamil(3,1)=t*kds*tmp2
      dhamil(3,1,1)=dt(1)*kds*tmp2
      dhamil(3,1,2)=dt(2)*kds*tmp2
      dhamil(3,1,3)=dt(3)*kds*tmp2

      call sxy('zx',dss,ddss,t,dt)                    !dzx-s
      hamil(4,1)=t*kds*tmp2
      dhamil(4,1,1)=dt(1)*kds*tmp2
      dhamil(4,1,2)=dt(2)*kds*tmp2
      dhamil(4,1,3)=dt(3)*kds*tmp2

      call sxxyy(dss,ddss,t,dt)                       !dx2_y2-s
      hamil(5,1)=t*kds*tmp2
      dhamil(5,1,1)=dt(1)*kds*tmp2
      dhamil(5,1,2)=dt(2)*kds*tmp2
      dhamil(5,1,3)=dt(3)*kds*tmp2

      call szz(dss,ddss,t,dt)                         !dz2-s
      hamil(6,1)=t*kds*tmp2
      dhamil(6,1,1)=dt(1)*kds*tmp2
      dhamil(6,1,2)=dt(2)*kds*tmp2
      dhamil(6,1,3)=dt(3)*kds*tmp2

      return
      end
