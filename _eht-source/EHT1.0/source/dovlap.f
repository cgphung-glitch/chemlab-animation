C The program is coded for Equation (8) in reference:
C   H.W. Jones, Int. J. Quantum Chem. 18, 709 (1980).
C The program basically can work for any combinations of the quantum numbers,
C but the C matrix coefficientd provided in the data block preclude that kind
C of job with the principal quantum number bigger than 4.

      subroutine dovlap(Nprime,Lprime,zetapm,N,L,M,zeta,S,dS)

      implicit none
      integer Nprime, Lprime,N,L,M,ns,k,i,j,tmp              !ns lower n
      double precision fact(0:12)
      double precision a,S,Keba,cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)
      double precision zeta,zetapm,tmpa,tmpb,ratioz,azeta,azetap
      double precision expzet,expzpm,azeadd,azminu,Kbrati
      double precision dS,dtempa,dtempb,r

      common /factorial/ fact
      common /cmatrix/ cmatrx
      common /distan/ r

      ratioz=(zetapm/zeta)**(Nprime+0.50)

      a=r/0.52917706

      if (r.gt.10.0) then
        S=0.0
        dS=0.0
        return
      endif

      azetap=a*zetapm
      azeta=a*zeta

      expzet=dexp(-azeta)
      expzpm=dexp(-azetap)

      azeadd=azetap+azeta
      azminu=azetap-azeta

      Keba=2**(Nprime+N)*(-1)**M

      tmpa=(2*L+1)*(2*Lprime+1)
      tmp=L+M
      tmpa=tmpa*fact(tmp)
      tmp=Lprime-M
      tmpa=tmpa*fact(tmp)
      tmp=2*Nprime
      tmpa=tmpa/fact(tmp)
      tmp=2*N
      tmpa=tmpa/fact(tmp)
      tmp=Lprime+M
      tmpa=tmpa/fact(tmp)
      tmp=L-M
      tmpa=tmpa/fact(tmp)

      keba=Keba*dsqrt(tmpa)

      S=0.0
      dS=0.0

      do 1 i=0,(N+L+Lprime)
        do 2 j=0,(N+Lprime)
          ns=Nprime-Lprime+j
          tmpa=0.0
          dtempa=0.0
          do 3 k=0,ns
            tmpb=(-1)**j/azminu**(k+1)
            tmpb=(-1)**i/azeadd**(k+1)-tmpb 
            tmp=ns-k
            tmpb=tmpb/fact(tmp)
C Prepare to compute the first derivative
            dtempb=-(k+1)*tmpb/a
            tmpa=tmpa+tmpb
            dtempa=dtempa+dtempb
 3        continue   
          tmpa=tmpa*expzpm
          dtempa=dtempa*expzpm-zetapm*tmpa
C
          tmpb=1.0/azeadd**(ns+1)
          tmpb=(-1)**j/azminu**(ns+1)-tmpb
          dtempb=-(ns+1)*tmpb/a
          tmpb=tmpb*expzet
          tmpa=tmpa+tmpb
          dtempb=dtempb*expzet-zeta*tmpb
          dtempa=dtempa+dtempb

          tmp=Nprime-2*Lprime-L+i+j
          tmpb=fact(ns)*cmatrx(Lprime,N,L,M,i,j)*azeta**tmp
          dtempb=tmp*tmpb/a

          S=S+tmpa*tmpb
          dS=dS+tmpa*dtempb+dtempa*tmpb

 2      continue
 1    continue
 
      Kbrati=Keba*ratioz

      S=S*Keba*ratioz*(-1)**(M+Lprime+L)
      dS=dS*Keba*ratioz*(-1)**(M+Lprime+L)
c
c Convert from bohr -1 to A -1
c
      dS=dS/0.52917706
c     if (M.eq.2) then
c     endif

      return
      end
