c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the Hamiltonian elements and the first    
c derivatives of these Hamiltonian elements between two same transition
c metal atoms.
c If you want to some calculation between two different metal atoms, find
c the parameters first.
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trtrhd (ip,hamil,dhamil)
c
      implicit none
      integer i,ip,jp
      double precision hamil(6,6),dhamil(6,6,3),t,dt(3)
      double precision ssseta0(6),rd(6),d(6),ssstau,sssalpha(6),sss,dsss
      double precision sdseta0,             sdstau,sdsalpha(6),sds,dsds
      double precision dsseta0,             dsstau,dssalpha(6),dss,ddss
      double precision ddseta0,             ddstau,ddsalpha(6),dds,ddds
      double precision ddpeta0,             ddptau,ddpalpha(6),ddp,dddp
      double precision dddeta0,             dddtau,dddalpha(6),ddd,dddd
c
c Hopping parameters for metal and metal
c
      DATA ( d(i),i=1,6) /2.4920,2.3978,3.0196,2.8890,2.6673,2.5292/
      DATA (rd(i),i=1,6) /0.71,0.67,0.94,0.89,1.04,1.01/
c s s sigma
      DATA (ssseta0(i),i=1,6) /-0.4238,-0.7083,-0.1992,-0.7475,
     &                         -0.5319,-0.4356/
      DATA ssstau /0.d0/
      DATA (sssalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c s d sigma
      DATA sdseta0 /-3.16d0/
      DATA sdstau /1.5d0/
      DATA (sdsalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c d s sigma
      DATA dsseta0 /-3.16d0/
      DATA dsstau /1.5d0/
      DATA (dssalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c d d sigma
      DATA ddseta0 /-16.2d0/
      DATA ddstau /3.d0/
      DATA (ddsalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c d d pi
      DATA ddpeta0 /8.75d0/
      DATA ddptau /3.d0/
      DATA (ddpalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c d d delta
      DATA dddeta0 /0.0d0/
      DATA dddtau /3.d0/
      DATA (dddalpha(i),i=1,6) /0.46463,0.43344,0.55001,0.32780,
     &                          0.66301,0.55180/
c
c Call subroutine to compute the two center integral of the two atoms when
c one atom is at the origin and the other atom is at z-axis.
c
      jp=ip-2
c
      call lathio(ssseta0(jp),rd(jp),d(jp),ssstau,sssalpha(jp),sss,dsss)
      call lathio(sdseta0,    rd(jp),d(jp),sdstau,sdsalpha(jp),sds,dsds)
      call lathio(dsseta0,    rd(jp),d(jp),dsstau,dssalpha(jp),dss,ddss)
      call lathio(ddseta0,    rd(jp),d(jp),ddstau,ddsalpha(jp),dds,ddds)
      call lathio(ddpeta0,    rd(jp),d(jp),ddptau,ddpalpha(jp),ddp,dddp)
      call lathio(dddeta0,    rd(jp),d(jp),dddtau,dddalpha(jp),ddd,dddd)
c 
c Project the two-center integral along the basis set.
c
      call ss(sss,dsss,t,dt)             !s-s
      hamil(1,1)=t
      dhamil(1,1,1)=dt(1)
      dhamil(1,1,2)=dt(2)
      dhamil(1,1,3)=dt(3)

      call sxy('xy',sds,dsds,t,dt)       !s-dxy
      hamil(1,2)=t
      dhamil(1,2,1)=dt(1)
      dhamil(1,2,2)=dt(2)
      dhamil(1,2,3)=dt(3)

      call sxy('yz',sds,dsds,t,dt)       !s-dyz
      hamil(1,3)=t
      dhamil(1,3,1)=dt(1)
      dhamil(1,3,2)=dt(2)
      dhamil(1,3,3)=dt(3)

      call sxy('zx',sds,dsds,t,dt)       !s-dzx
      hamil(1,4)=t
      dhamil(1,4,1)=dt(1)
      dhamil(1,4,2)=dt(2)
      dhamil(1,4,3)=dt(3)

      call sxxyy(sds,dsds,t,dt)          !s-dx2_y2
      hamil(1,5)=t
      dhamil(1,5,1)=dt(1)
      dhamil(1,5,2)=dt(2)
      dhamil(1,5,3)=dt(3)

      call szz(sds,dsds,t,dt)            !s-dz2
      hamil(1,6)=t
      dhamil(1,6,1)=dt(1)
      dhamil(1,6,2)=dt(2)
      dhamil(1,6,3)=dt(3)

      call sxy('xy',dss,ddss,t,dt)             !dxy-s
      hamil(2,1)=t
      dhamil(2,1,1)=dt(1)
      dhamil(2,1,2)=dt(2)
      dhamil(2,1,3)=dt(3)

      call xyxy('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dxy
      hamil(2,2)=t
      dhamil(2,2,1)=dt(1)
      dhamil(2,2,2)=dt(2)
      dhamil(2,2,3)=dt(3)

      call xyyz('xz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dyz
      hamil(2,3)=t
      dhamil(2,3,1)=dt(1)
      dhamil(2,3,2)=dt(2)
      dhamil(2,3,3)=dt(3)

      call xyyz('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dxy-dzx
      hamil(2,4)=t
      dhamil(2,4,1)=dt(1)
      dhamil(2,4,2)=dt(2)
      dhamil(2,4,3)=dt(3)

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dxy-dx2_y2
      hamil(2,5)=t
      dhamil(2,5,1)=dt(1)
      dhamil(2,5,2)=dt(2)
      dhamil(2,5,3)=dt(3)

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      hamil(2,6)=t
      dhamil(2,6,1)=dt(1)
      dhamil(2,6,2)=dt(2)
      dhamil(2,6,3)=dt(3)

      call sxy('yz',dss,ddss,t,dt)                     !dyz-s
      hamil(3,1)=t
      dhamil(3,1,1)=dt(1)
      dhamil(3,1,2)=dt(2)
      dhamil(3,1,3)=dt(3)

      call xyyz('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dxy
      hamil(3,2)=t
      dhamil(3,2,1)=dt(1)
      dhamil(3,2,2)=dt(2)
      dhamil(3,2,3)=dt(3)

      call xyxy('yz',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dyz
      hamil(3,3)=t
      dhamil(3,3,1)=dt(1)
      dhamil(3,3,2)=dt(2)
      dhamil(3,3,3)=dt(3)

      call xyyz('yx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dyz-dzx
      hamil(3,4)=t
      dhamil(3,4,1)=dt(1)
      dhamil(3,4,2)=dt(2)
      dhamil(3,4,3)=dt(3)

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dxy-dx2_y2
      hamil(3,5)=t
      dhamil(3,5,1)=dt(1)
      dhamil(3,5,2)=dt(2)
      dhamil(3,5,3)=dt(3)

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dxy-dz2
      hamil(3,6)=t
      dhamil(3,6,1)=dt(1)
      dhamil(3,6,2)=dt(2)
      dhamil(3,6,3)=dt(3)

      call sxy('zx',dss,ddss,t,dt)                     !dzx-s
      hamil(4,1)=t
      dhamil(4,1,1)=dt(1)
      dhamil(4,1,2)=dt(2)
      dhamil(4,1,3)=dt(3)

      call xyyz('zy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dxy
      hamil(4,2)=t
      dhamil(4,2,1)=dt(1)
      dhamil(4,2,2)=dt(2)
      dhamil(4,2,3)=dt(3)

      call xyyz('xy',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dyz
      hamil(4,3)=t
      dhamil(4,3,1)=dt(1)
      dhamil(4,3,2)=dt(2)
      dhamil(4,3,3)=dt(3)

      call xyxy('zx',dds,ddds,ddp,dddp,ddd,dddd,t,dt)  !dzx-dzx
      hamil(4,4)=t
      dhamil(4,4,1)=dt(1)
      dhamil(4,4,2)=dt(2)
      dhamil(4,4,3)=dt(3)

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      hamil(4,5)=t
      dhamil(4,5,1)=dt(1)
      dhamil(4,5,2)=dt(2)
      dhamil(4,5,3)=dt(3)

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dzx-dz2
      hamil(4,6)=t
      dhamil(4,6,1)=dt(1)
      dhamil(4,6,2)=dt(2)
      dhamil(4,6,3)=dt(3)

      call sxxyy(sds,dsds,t,dt)                       !dx2_y2-s
      hamil(5,1)=t
      dhamil(5,1,1)=dt(1)
      dhamil(5,1,2)=dt(2)
      dhamil(5,1,3)=dt(3)

      call xyxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      hamil(5,2)=t
      dhamil(5,2,1)=dt(1)
      dhamil(5,2,2)=dt(2)
      dhamil(5,2,3)=dt(3)

      call yzxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dxy
      hamil(5,3)=t
      dhamil(5,3,1)=dt(1)
      dhamil(5,3,2)=dt(2)
      dhamil(5,3,3)=dt(3)

      call zxxxyy(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dzx
      hamil(5,4)=t
      dhamil(5,4,1)=dt(1)
      dhamil(5,4,2)=dt(2)
      dhamil(5,4,3)=dt(3)

      call xxyy_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dx2_y2-dx2_y2
      hamil(5,5)=t
      dhamil(5,5,1)=dt(1)
      dhamil(5,5,2)=dt(2)
      dhamil(5,5,3)=dt(3)

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dzx-dx2_y2
      hamil(5,6)=t
      dhamil(5,6,1)=dt(1)
      dhamil(5,6,2)=dt(2)
      dhamil(5,6,3)=dt(3)

      call szz(sds,dsds,t,dt)                          !dz2-s
      hamil(6,1)=t
      dhamil(6,1,1)=dt(1)
      dhamil(6,1,2)=dt(2)
      dhamil(6,1,3)=dt(3)

      call xyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      hamil(6,2)=t
      dhamil(6,2,1)=dt(1)
      dhamil(6,2,2)=dt(2)
      dhamil(6,2,3)=dt(3)

      call yzzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dxy
      hamil(6,3)=t
      dhamil(6,3,1)=dt(1)
      dhamil(6,3,2)=dt(2)
      dhamil(6,3,3)=dt(3)

      call zxzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dzx
      hamil(6,4)=t
      dhamil(6,4,1)=dt(1)
      dhamil(6,4,2)=dt(2)
      dhamil(6,4,3)=dt(3)

      call xxyyzz(dds,ddds,ddp,dddp,ddd,dddd,t,dt)     !dz2-dx2_y2
      hamil(6,5)=t
      dhamil(6,5,1)=dt(1)
      dhamil(6,5,2)=dt(2)
      dhamil(6,5,3)=dt(3)

      call zz_2(dds,ddds,ddp,dddp,ddd,dddd,t,dt)       !dz2-dz2
      hamil(6,6)=t
      dhamil(6,6,1)=dt(1)
      dhamil(6,6,2)=dt(2)
      dhamil(6,6,3)=dt(3)

      return
      end
