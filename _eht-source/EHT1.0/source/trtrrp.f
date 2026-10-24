c 88888888888888888888888888888888888888888888888888888888888888888888888
c This subroutine calculates the pairwise repulsion and the first    
c derivatives of the repulsive energy between Ni and Ni
c 88888888888888888888888888888888888888888888888888888888888888888888888
c
      subroutine trtrrp (ip,repul,drepul)
c
      implicit none
      integer i,ip,jp
      double precision repul,dt,drepul(3)
      double precision l,m,n,lladmm,llmimm,nn
      double precision alpha(6),gamma(6),beta(6)
c
      common /lmn/ l,m,n,lladmm,llmimm,nn
c
c Repulsive parameters for two same transition metal atom.
c
      DATA (alpha(i),i=1,6) /5018.6 ,88898.0,40079.0,104820,
     &                       25219.0,57529.0/
      DATA (gamma(i),i=1,6) /0.80114,0.46178,0.35787,2.29571,
     &                       1.82081,0.74050/
      DATA (beta(i), i=1,6) /4.3383 ,6.0109,4.0670 ,5.6531 ,
     &                       4.5748,4.6922/
c
      jp=ip-2
c
      call grepfm (alpha(jp),gamma(jp),beta(jp),repul,dt)
      
      drepul(1)=dt*l
      drepul(2)=dt*m
      drepul(3)=dt*n
 
      return
      end
