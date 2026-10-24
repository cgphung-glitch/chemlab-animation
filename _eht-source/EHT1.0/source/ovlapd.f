C The program is coded for Equation (8) in reference:
C   H.W. Jones, Int. J. Quantum Chem. 19, 567 (1981).
C The program basically can work for any combinations of the quantum numbers,
C but the C matrix coefficients provided in the data block preclude that kind
C of job with the principal quantum number bigger than 4.
C For the same zeta.

      subroutine ovlapd(CN,CL,zeta,n,l,M,S,dS)
  
      implicit none
      integer CN,CL,n,l,M,i,j,t,u,v,k,w,tmp
      double precision fact(0:12)
      double precision zeta,S,dS,r,a,p,Keba
      double precision cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)
      double precision tmpa,tmpb,tmpc

      common /factorial/ fact
      common /cmatrix/ cmatrx
      common /distan/ r

      a=r/0.52917706
      p=zeta*a
 
      if (r.gt.10.0) then
        S=0.0
        dS=0.0
        return
      endif

      Keba=2**(CN+n)*(-1)**M 

      tmpa=(2*CL+1)*(2*l+1)
      tmp=l+M
      tmpa=tmpa*fact(tmp)
      tmp=CL-M
      tmpa=tmpa*fact(tmp)
      tmp=2*CN
      tmpa=tmpa/fact(tmp)
      tmp=2*n
      tmpa=tmpa/fact(tmp)
      tmp=CL+M
      tmpa=tmpa/fact(tmp)
      tmp=l-M
      tmpa=tmpa/fact(tmp)

      Keba=Keba*sqrt(tmpa)

      S=0.0
      dS=0.0

      do 1 i=0,n+CL+l
        do 2 j=0,n+CL
          t=j+CN-CL
          u=j+CN-2*CL-l+i
          v=i-CL-l-1
          S=S+cmatrx(CL,n,l,M,i,j)*
     &       ((-1)**j*p**u/(t+1)-fact(t)*p**v/2**(t+1))
          dS=dS+cmatrx(CL,n,l,M,i,j)*zeta*
     &       (u*(-1)**j*p**(u-1)/(t+1)-v*fact(t)*p**(v-1)/2**(t+1))
          tmpa=0.0
          tmpb=0.0
          do 3 k=0,t
            tmp=t-k
            w=t-k+i-l-CL-1
c           print 99,t,u,v,w
c99    format("t,u,v,and w ",4(I2,1x))
            tmpa=tmpa+fact(t)/fact(tmp)*p**w/2**(k+1)
            tmpb=tmpb+w*fact(t)*zeta/fact(tmp)*p**(w-1)/2**(k+1)
 3        continue
c         print *,"tmpa tmpb ",tmpa,tmpb
          S=S+cmatrx(CL,n,l,M,i,j)*(-1)**i*tmpa
          dS=dS+cmatrx(CL,n,l,M,i,j)*(-1)**i*tmpb
 2      continue
 1    continue

      dS=Keba*(exp(-p)*dS-zeta*exp(-p)*S)*(-1)**(M+CL+l)
      S=Keba*exp(-p)*S*(-1)**(M+CL+l)

c
c Convert from bohr -1 to A -1
c
      dS=dS/0.52917706

      return
      end
