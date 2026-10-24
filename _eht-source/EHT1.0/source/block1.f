CC Note: The following section gives a fortran program to compute
CC coefficients by removing "C "(C and a space).

CC This program is made to obtain the C matrix  for a given l,N,L,M.
CC The formula for C matrix is given in reference:
CC H. W. Jones and C.A. Weatherford, Int. J. Quantum Chem. Symp. 12, 483 (1978).
CC In Page 486 they wrote: We build up a C matrix Cl(NA,NR) using eight nested 
CC DO loops to sort and collect the sums of all coefficients of all like powers
CC of r and a.
CC
CC This program only gives all C matrix coefficients, the modulus of 
CC N-L+2p+2q+2qprime-k-kprime to 2, the exponential index of a, and that of r.
CC After this program is compiled with f77, you can obtain a series of data
CC for C matrix.
C
C      program gen
C
C      implicit none
C      integer sl,N,L,M,p,q,v,pp,qp,vp,t,k,kp
C      integer result2,result3,result4,tmp
C      double precision result,alpha,beta,efunc
C       
C      print *, "Please input l, N, L, M"
C      read *,sl,N,L,M
C      if (sl .lt. M) stop "Irrational pair: l and M"
C      print 98,sl,N,L,M
C      do 1 p=0,(L+M)/2
C        do 2 q=0,(L+M-2*p)
C          do 3 v=0,(L+M-2*p-q)
C            do 4 pp=0,((sl-M)/2)
C              do 5 qp=0,sl-M-2*pp
C                do 6 vp=0,sl-M-2*pp-qp
C                  t=N-L+2*p+2*q+2*qp
C                  do 7 k=0,t
C                    do 8 kp=0,t-k
C                      result=(-1)**(L+v+qp)
C                      result=result/4.0**(L+sl-p-pp)
C                      tmp=L-M
C                      result=result*alpha(L,tmp,p)
C                      tmp=L+M-2*p
C                      result=result*beta(tmp,q)
C                      tmp=L+M-2*p-q
C                      result=result*beta(tmp,v)
C                      tmp=sl+M
C                      result=result*alpha(sl,tmp,pp)
C                      tmp=sl-M-2*pp
C                      result=result*beta(tmp,qp)
C                      tmp=sl-M-2*pp-qp
C                      result=result*beta(tmp,vp)
C                      tmp=N-L+2*p+2*q+2*qp
C                      result=result*efunc(tmp,k)
C                      tmp=N-L+2*p+2*q+2*qp-k
C                      result=result*beta(tmp,kp)
C                 
C                      result2=N-L+2*p+2*q+2*qp-k-kp
C                      result2=MOD(result2,2)
C
C                      result3=N+L+2*sl-2*pp-2*vp-2*v-k-kp
C
C                      result4=2*pp+2*v+2*vp+kp
C                     
C                      print 99, result,result2,result3,result4
C
C 8                  continue
C 7                continue
C 6              continue
C 5            continue
C 4          continue
C 3        continue
C 2      continue
C 1    continue
C
C 98   format(5x,"tttt ",4(I4,2x))
C 99   format(5x,f20.8,3(I4,2x))
C
C      end
C       
C
C      double precision function alpha(sl,n,p)
C      implicit none
C      integer sl,n,p,tmpa,tmpb,tmpc
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=2*sl-2*p
C      tmpb=2*sl-2*p-n
C      tmpc=sl-p
C
C      alpha=(-1)**p*fact(tmpa)/(fact(tmpb)*fact(tmpc)*fact(p))
C
C      return
C      end
C
C      double precision function beta(n,k)
C      implicit none
C      integer n,k,tmpa
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=n-k
C
C      beta=fact(n)/(fact(tmpa)*fact(k))
C
C      return
C      end
C
C      double precision function efunc(t,k)
C      implicit none
C      integer t,k,tmpa
C      double precision fact(0:14)
C      common /factorial/ fact
C
C      tmpa=t-k
C
C      efunc=fact(t)/fact(tmpa)
C
C      return
C      end
C
C      BLOCK DATA
C      double precision fact(0:14)
C      common /factorial/ fact
C      data fact/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
C     &          362880.0,3628800.0,39916800.0,479001600.0,622020800.0,
C     &          87178291200.0/
C      end
C
CC Suppose the executable file of this program is a.out and you are in a
CC UNIX system, the following commands let you sort and collect the sums of 
CC all coefficients of all like powers C of r and a.
CC     a.out | sort +2 -4 | awk -f awk.prog 
CC awk.prog is a file to do some special job for unix command awk.
CC File awk.prog begins.
CC /ttt/ { sl=$2;n=$3;l=$4;m=$5}
CC {pop[$3,$4] += $1}
CC END { printf("\nC Those are data for Lprime,N,L,M %d%d%d%d.\n\n",sl,n,l,m);
CC   for (i=0;i<=(sl+n+l);i++) 
CC     for (j=0;j<=(sl+n);j++)  
CC       printf("      DATA cmatrx(%2d,%2d,%2d,%2d,%2d,%2d)/%20.8f/\n",sl,n,l,
CC         m,i,j,pop[i,j]);
CC }
CC File awk.prog ends.
CC This is how the all above data are produced.

C C This is how the all following data are produced.

C Main data begins here.

       BLOCK DATA set1

      implicit none

      double precision cmatrx(0:2,1:6,0:2,0:2,0:10,0:8)

      common /cmatrix/ cmatrx

C This is all the data block for the C matrix elements
C   l (lprime): 0-2
C   N         : 1-6
C   L         : 0-2
C   M         : 0-2
C If you want to calculate the overlap integral involving principal quantum
C number 7 or lagger, go ahead made your extra C matrix.
C I have some little program and notes followed those data.
C Advice: Do not trust me too much.

C Those are data for Lprime,N,L,M 0100.

      DATA cmatrx( 0, 1, 0, 0, 0, 0)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 0, 1)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 1, 0)/          1.00000000/
      DATA cmatrx( 0, 1, 0, 0, 1, 1)/          0.00000000/

C Those are data for Lprime,N,L,M 0200.

      DATA cmatrx( 0, 2, 0, 0, 0, 0)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 0, 1)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 0, 2)/          1.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 0)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 1)/          2.00000000/
      DATA cmatrx( 0, 2, 0, 0, 1, 2)/          0.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 0)/          1.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 1)/          0.00000000/
      DATA cmatrx( 0, 2, 0, 0, 2, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 0210.

      DATA cmatrx( 0, 2, 1, 0, 0, 0)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 0, 1)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 0, 2)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 0)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 1)/         -3.00000000/
      DATA cmatrx( 0, 2, 1, 0, 1, 2)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 0)/         -2.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 1)/         -2.00000000/
      DATA cmatrx( 0, 2, 1, 0, 2, 2)/          0.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 0)/         -1.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 1)/          0.00000000/
      DATA cmatrx( 0, 2, 1, 0, 3, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 0300.

      DATA cmatrx( 0, 3, 0, 0, 0, 0)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 1)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 0, 3)/          1.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 0)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 1)/          6.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 1, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 0)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 1)/          3.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 0)/          1.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0310.

      DATA cmatrx( 0, 3, 1, 0, 0, 0)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 1)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 2)/         -5.00000000/
      DATA cmatrx( 0, 3, 1, 0, 0, 3)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 0)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 1)/        -12.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 2)/         -5.00000000/
      DATA cmatrx( 0, 3, 1, 0, 1, 3)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 0)/         -7.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 1)/         -7.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 2)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 0)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 1)/         -3.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 0)/         -1.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 1, 0, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0320.

      DATA cmatrx( 0, 3, 2, 0, 0, 0)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 1)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 2)/         18.00000000/
      DATA cmatrx( 0, 3, 2, 0, 0, 3)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 0)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 1)/         45.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 2)/         18.00000000/
      DATA cmatrx( 0, 3, 2, 0, 1, 3)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 0)/         24.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 1)/         24.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 2)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 0)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 1)/          9.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 2)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 0)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 1)/          3.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 3, 2, 0, 5, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 0400.

      DATA cmatrx( 0, 4, 0, 0, 0, 0)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 1)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 2)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 0, 4)/          1.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 0)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 1)/         24.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 2)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 1, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 0)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 1)/         12.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 2)/          6.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 0)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 1)/          4.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 0)/          1.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0410.

      DATA cmatrx( 0, 4, 1, 0, 0, 0)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 1)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 2)/        -27.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 3)/         -7.00000000/
      DATA cmatrx( 0, 4, 1, 0, 0, 4)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 0)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 1)/        -60.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 2)/        -27.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 3)/         -7.00000000/
      DATA cmatrx( 0, 4, 1, 0, 1, 4)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 0)/        -33.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 1)/        -33.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 2)/        -15.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 3)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 0)/        -13.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 1)/        -13.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 2)/         -6.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 0)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 1)/         -4.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 0)/         -1.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0420.

      DATA cmatrx( 0, 4, 2, 0, 0, 0)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 1)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 2)/        117.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 3)/         27.00000000/
      DATA cmatrx( 0, 4, 2, 0, 0, 4)/          3.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 0)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 1)/        270.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 2)/        117.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 3)/         27.00000000/
      DATA cmatrx( 0, 4, 2, 0, 1, 4)/          3.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 0)/        141.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 1)/        141.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 2)/         60.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 3)/         13.00000000/
      DATA cmatrx( 0, 4, 2, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 0)/         51.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 1)/         51.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 2)/         21.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 3)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 0)/         15.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 1)/         15.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 2)/          6.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 0)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 1)/          4.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 4, 2, 0, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 0500.

      DATA cmatrx( 0, 5, 0, 0, 0, 0)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 1)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 2)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 3)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 0, 5)/          1.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 0)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 1)/        120.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 2)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 3)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 1, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 0)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 1)/         60.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 2)/         30.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 3)/         10.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 0)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 1)/         20.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 2)/         10.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 0)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 1)/          5.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0510.

      DATA cmatrx( 0, 5, 1, 0, 0, 0)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 1)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 2)/       -168.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 3)/        -48.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 0, 5, 1, 0, 0, 5)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 0)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 1)/       -360.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 2)/       -168.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 3)/        -48.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 0, 5, 1, 0, 1, 5)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 0)/       -192.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 1)/       -192.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 2)/        -90.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 3)/        -26.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 4)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 0)/        -72.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 1)/        -72.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 2)/        -34.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 3)/        -10.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 0)/        -21.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 1)/        -21.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 2)/        -10.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 0)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 1)/         -5.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0520.

      DATA cmatrx( 0, 5, 2, 0, 0, 0)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 1)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 2)/        855.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 3)/        225.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 4)/         36.00000000/
      DATA cmatrx( 0, 5, 2, 0, 0, 5)/          3.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 0)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 1)/       1890.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 2)/        855.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 3)/        225.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 4)/         36.00000000/
      DATA cmatrx( 0, 5, 2, 0, 1, 5)/          3.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 0)/        975.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 1)/        975.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 2)/        438.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 3)/        113.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 4)/         17.00000000/
      DATA cmatrx( 0, 5, 2, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 0)/        345.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 1)/        345.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 2)/        153.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 3)/         38.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 4)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 0)/         96.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 1)/         96.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 2)/         42.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 3)/         10.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 0)/         23.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 1)/         23.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 2)/         10.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 0)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 1)/          5.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 0)/          1.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 1)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 5, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 0600.

      DATA cmatrx( 0, 6, 0, 0, 0, 0)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 1)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 2)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 3)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 4)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 0, 6)/          1.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 0)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 1)/        720.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 2)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 3)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 4)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 1, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 0)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 1)/        360.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 2)/        180.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 3)/         60.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 4)/         15.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 0)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 1)/        120.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 2)/         60.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 3)/         20.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 0)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 1)/         30.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 2)/         15.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 0)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 1)/          6.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 0610.

      DATA cmatrx( 0, 6, 1, 0, 0, 0)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 1)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 2)/      -1200.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 3)/       -360.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 4)/        -75.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 5)/        -11.00000000/
      DATA cmatrx( 0, 6, 1, 0, 0, 6)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 0)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 1)/      -2520.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 2)/      -1200.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 3)/       -360.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 4)/        -75.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 5)/        -11.00000000/
      DATA cmatrx( 0, 6, 1, 0, 1, 6)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 0)/      -1320.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 1)/      -1320.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 2)/       -630.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 3)/       -190.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 4)/        -40.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 5)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 0)/       -480.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 1)/       -480.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 2)/       -230.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 3)/        -70.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 4)/        -15.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 0)/       -135.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 1)/       -135.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 2)/        -65.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 3)/        -20.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 0)/        -31.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 1)/        -31.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 2)/        -15.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 0)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 1)/         -6.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 0620.

      DATA cmatrx( 0, 6, 2, 0, 0, 0)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 1)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 2)/       7020.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 3)/       1980.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 4)/        369.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 5)/         45.00000000/
      DATA cmatrx( 0, 6, 2, 0, 0, 6)/          3.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 0)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 1)/      15120.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 2)/       7020.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 3)/       1980.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 4)/        369.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 5)/         45.00000000/
      DATA cmatrx( 0, 6, 2, 0, 1, 6)/          3.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 0)/       7740.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 1)/       7740.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 2)/       3582.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 3)/       1002.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 4)/        183.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 5)/         21.00000000/
      DATA cmatrx( 0, 6, 2, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 0)/       2700.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 1)/       2700.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 2)/       1242.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 3)/        342.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 4)/         60.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 5)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 0)/        729.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 1)/        729.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 2)/        333.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 3)/         90.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 4)/         15.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 0)/        165.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 1)/        165.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 2)/         75.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 3)/         20.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 0)/         33.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 1)/         33.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 2)/         15.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 0)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 1)/          6.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 0)/          1.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 1)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 0, 6, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1100.

      DATA cmatrx( 1, 1, 0, 0, 0, 0)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 0, 1)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 0, 2)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 0)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 1)/         -3.00000000/
      DATA cmatrx( 1, 1, 0, 0, 1, 2)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 0)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 1)/         -1.00000000/
      DATA cmatrx( 1, 1, 0, 0, 2, 2)/          0.00000000/

C Those are data for Lprime,N,L,M 1200.

      DATA cmatrx( 1, 2, 0, 0, 0, 0)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 1)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 2)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 0, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 0)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 1)/        -12.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 2)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 1, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 0)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 1)/         -5.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 2)/         -2.00000000/
      DATA cmatrx( 1, 2, 0, 0, 2, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 0)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 1)/         -1.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1210.

      DATA cmatrx( 1, 2, 1, 0, 0, 0)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 1)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 2)/         12.00000000/
      DATA cmatrx( 1, 2, 1, 0, 0, 3)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 0)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 1)/         30.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 2)/         12.00000000/
      DATA cmatrx( 1, 2, 1, 0, 1, 3)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 0)/         15.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 1)/         15.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 2)/          6.00000000/
      DATA cmatrx( 1, 2, 1, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 0)/          5.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 1)/          5.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 2)/          2.00000000/
      DATA cmatrx( 1, 2, 1, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 0)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 1)/          1.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 0, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1211.

      DATA cmatrx( 1, 2, 1, 1, 0, 0)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 1)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 2)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 0, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 0)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 1)/        -15.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 2)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 1, 3)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 0)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 1)/         -6.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 2)/         -2.00000000/
      DATA cmatrx( 1, 2, 1, 1, 2, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 0)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 1)/         -1.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 0)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 1)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 2, 1, 1, 4, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 1300.

      DATA cmatrx( 1, 3, 0, 0, 0, 0)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 1)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 2)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 3)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 0, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 0)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 1)/        -60.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 2)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 3)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 1, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 0)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 1)/        -27.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 2)/        -12.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 0, 0, 2, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 0)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 1)/         -7.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1310.

      DATA cmatrx( 1, 3, 1, 0, 0, 0)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 1)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 2)/         78.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 3)/         18.00000000/
      DATA cmatrx( 1, 3, 1, 0, 0, 4)/          2.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 0)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 1)/        180.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 2)/         78.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 3)/         18.00000000/
      DATA cmatrx( 1, 3, 1, 0, 1, 4)/          2.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 0)/         90.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 1)/         90.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 2)/         39.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 3)/          9.00000000/
      DATA cmatrx( 1, 3, 1, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 0)/         30.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 1)/         30.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 2)/         13.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 3)/          3.00000000/
      DATA cmatrx( 1, 3, 1, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 0)/          7.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 1)/          7.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 2)/          3.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 0)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 1)/          1.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1311.

      DATA cmatrx( 1, 3, 1, 1, 0, 0)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 1)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 2)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 3)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 0, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 0)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 1)/        -90.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 2)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 3)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 1, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 0)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 1)/        -39.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 2)/        -16.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 1, 1, 2, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 0)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 1)/         -9.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 0)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 1)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 1, 1, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1320.

      DATA cmatrx( 1, 3, 2, 0, 0, 0)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 1)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 2)/       -405.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 3)/        -90.00000000/
      DATA cmatrx( 1, 3, 2, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 0)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 1)/       -945.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 2)/       -405.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 3)/        -90.00000000/
      DATA cmatrx( 1, 3, 2, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 0)/       -465.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 1)/       -465.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 2)/       -198.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 3)/        -43.00000000/
      DATA cmatrx( 1, 3, 2, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 0)/       -150.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 1)/       -150.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 2)/        -63.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 3)/        -13.00000000/
      DATA cmatrx( 1, 3, 2, 0, 3, 4)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 0)/        -36.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 1)/        -36.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 2)/        -15.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 3)/         -3.00000000/
      DATA cmatrx( 1, 3, 2, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 0)/         -7.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 1)/         -7.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 2)/         -3.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 0, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1321.

      DATA cmatrx( 1, 3, 2, 1, 0, 0)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 1)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 2)/        135.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 3)/         30.00000000/
      DATA cmatrx( 1, 3, 2, 1, 0, 4)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 0)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 1)/        315.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 2)/        135.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 3)/         30.00000000/
      DATA cmatrx( 1, 3, 2, 1, 1, 4)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 0)/        150.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 1)/        150.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 2)/         63.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 3)/         13.00000000/
      DATA cmatrx( 1, 3, 2, 1, 2, 4)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 0)/         45.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 1)/         45.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 2)/         18.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 3)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 0)/          9.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 1)/          9.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 2)/          3.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 0)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 1)/          1.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 3, 2, 1, 6, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 1400.

      DATA cmatrx( 1, 4, 0, 0, 0, 0)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 1)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 2)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 3)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 0, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 0)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 1)/       -360.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 2)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 3)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 1, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 0)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 1)/       -168.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 2)/        -78.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 3)/        -22.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 0, 0, 2, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 0)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 1)/        -48.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 2)/        -22.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 0)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 1)/         -9.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1410.

      DATA cmatrx( 1, 4, 1, 0, 0, 0)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 1)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 2)/        570.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 3)/        150.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 4)/         24.00000000/
      DATA cmatrx( 1, 4, 1, 0, 0, 5)/          2.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 0)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 1)/       1260.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 2)/        570.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 3)/        150.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 4)/         24.00000000/
      DATA cmatrx( 1, 4, 1, 0, 1, 5)/          2.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 0)/        630.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 1)/        630.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 2)/        285.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 3)/         75.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 4)/         12.00000000/
      DATA cmatrx( 1, 4, 1, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 0)/        210.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 1)/        210.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 2)/         95.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 3)/         25.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 4)/          4.00000000/
      DATA cmatrx( 1, 4, 1, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 0)/         51.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 1)/         51.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 2)/         23.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 3)/          6.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 0)/          9.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 1)/          9.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 2)/          4.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 0)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 1)/          1.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1411.

      DATA cmatrx( 1, 4, 1, 1, 0, 0)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 1)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 2)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 3)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 4)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 0, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 0)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 1)/       -630.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 2)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 3)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 4)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 1, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 0)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 1)/       -285.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 2)/       -126.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 3)/        -31.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 1, 1, 2, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 0)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 1)/        -75.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 2)/        -31.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 0)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 1)/        -12.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 1, 1, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1420.

      DATA cmatrx( 1, 4, 2, 0, 0, 0)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 1)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 2)/      -3375.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 3)/       -855.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 4)/       -126.00000000/
      DATA cmatrx( 1, 4, 2, 0, 0, 5)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 0)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 1)/      -7560.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 2)/      -3375.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 3)/       -855.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 4)/       -126.00000000/
      DATA cmatrx( 1, 4, 2, 0, 1, 5)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 0)/      -3735.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 1)/      -3735.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 2)/      -1662.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 3)/       -417.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 4)/        -60.00000000/
      DATA cmatrx( 1, 4, 2, 0, 2, 5)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 0)/      -1215.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 1)/      -1215.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 2)/       -537.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 3)/       -132.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 4)/        -18.00000000/
      DATA cmatrx( 1, 4, 2, 0, 3, 5)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 0)/       -294.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 1)/       -294.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 2)/       -129.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 3)/        -31.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 4)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 0)/        -57.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 1)/        -57.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 2)/        -25.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 3)/         -6.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 0)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 1)/         -9.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 2)/         -4.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1421.

      DATA cmatrx( 1, 4, 2, 1, 0, 0)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 1)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 2)/       1125.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 3)/        285.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 4)/         42.00000000/
      DATA cmatrx( 1, 4, 2, 1, 0, 5)/          3.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 0)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 1)/       2520.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 2)/       1125.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 3)/        285.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 4)/         42.00000000/
      DATA cmatrx( 1, 4, 2, 1, 1, 5)/          3.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 0)/       1215.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 1)/       1215.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 2)/        537.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 3)/        132.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 4)/         18.00000000/
      DATA cmatrx( 1, 4, 2, 1, 2, 5)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 0)/        375.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 1)/        375.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 2)/        162.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 3)/         37.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 4)/          4.00000000/
      DATA cmatrx( 1, 4, 2, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 0)/         81.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 1)/         81.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 2)/         33.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 3)/          6.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 0)/         12.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 1)/         12.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 2)/          4.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 0)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 1)/          1.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 4, 2, 1, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 1500.

      DATA cmatrx( 1, 5, 0, 0, 0, 0)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 1)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 2)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 3)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 4)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 5)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 0, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 0)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 1)/      -2520.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 2)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 3)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 4)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 5)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 1, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 0)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 1)/      -1200.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 2)/       -570.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 3)/       -170.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 4)/        -35.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 0, 0, 2, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 0)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 1)/       -360.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 2)/       -170.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 3)/        -50.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 0)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 1)/        -75.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 2)/        -35.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 0)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 1)/        -11.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1510.

      DATA cmatrx( 1, 5, 1, 0, 0, 0)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 1)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 2)/       4680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 3)/       1320.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 4)/        246.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 5)/         30.00000000/
      DATA cmatrx( 1, 5, 1, 0, 0, 6)/          2.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 0)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 1)/      10080.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 2)/       4680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 3)/       1320.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 4)/        246.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 5)/         30.00000000/
      DATA cmatrx( 1, 5, 1, 0, 1, 6)/          2.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 0)/       5040.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 1)/       5040.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 2)/       2340.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 3)/        660.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 4)/        123.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 5)/         15.00000000/
      DATA cmatrx( 1, 5, 1, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 0)/       1680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 1)/       1680.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 2)/        780.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 3)/        220.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 4)/         41.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 5)/          5.00000000/
      DATA cmatrx( 1, 5, 1, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 0)/        414.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 1)/        414.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 2)/        192.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 3)/         54.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 4)/         10.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 0)/         78.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 1)/         78.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 2)/         36.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 3)/         10.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 0)/         11.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 1)/         11.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 2)/          5.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 0)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 1)/          1.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1511.

      DATA cmatrx( 1, 5, 1, 1, 0, 0)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 1)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 2)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 3)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 4)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 5)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 0, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 0)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 1)/      -5040.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 2)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 3)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 4)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 5)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 1, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 0)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 1)/      -2340.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 2)/      -1074.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 3)/       -294.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 4)/        -51.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 1, 1, 2, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 0)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 1)/       -660.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 2)/       -294.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 3)/        -74.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 0)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 1)/       -123.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 2)/        -51.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 0)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 1)/        -15.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 1, 1, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1520.

      DATA cmatrx( 1, 5, 2, 0, 0, 0)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 1)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 2)/     -31185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 3)/      -8505.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 4)/      -1485.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 5)/       -162.00000000/
      DATA cmatrx( 1, 5, 2, 0, 0, 6)/         -9.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 0)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 1)/     -68040.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 2)/     -31185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 3)/      -8505.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 4)/      -1485.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 5)/       -162.00000000/
      DATA cmatrx( 1, 5, 2, 0, 1, 6)/         -9.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 0)/     -33705.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 1)/     -33705.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 2)/     -15420.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 3)/      -4185.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 4)/       -723.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 5)/        -77.00000000/
      DATA cmatrx( 1, 5, 2, 0, 2, 6)/         -4.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 0)/     -11025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 1)/     -11025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 2)/      -5025.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 3)/      -1350.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 4)/       -228.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 5)/        -23.00000000/
      DATA cmatrx( 1, 5, 2, 0, 3, 6)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 0)/      -2685.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 1)/      -2685.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 2)/      -1218.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 3)/       -323.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 4)/        -53.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 5)/         -5.00000000/
      DATA cmatrx( 1, 5, 2, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 0)/       -522.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 1)/       -522.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 2)/       -236.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 3)/        -62.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 4)/        -10.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 0)/        -84.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 1)/        -84.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 2)/        -38.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 3)/        -10.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 0)/        -11.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 1)/        -11.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 2)/         -5.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 0)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 1)/         -1.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1521.

      DATA cmatrx( 1, 5, 2, 1, 0, 0)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 1)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 2)/      10395.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 3)/       2835.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 4)/        495.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 5)/         54.00000000/
      DATA cmatrx( 1, 5, 2, 1, 0, 6)/          3.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 0)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 1)/      22680.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 2)/      10395.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 3)/       2835.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 4)/        495.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 5)/         54.00000000/
      DATA cmatrx( 1, 5, 2, 1, 1, 6)/          3.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 0)/      11025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 1)/      11025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 2)/       5025.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 3)/       1350.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 4)/        228.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 5)/         23.00000000/
      DATA cmatrx( 1, 5, 2, 1, 2, 6)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 0)/       3465.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 1)/       3465.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 2)/       1560.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 3)/        405.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 4)/         63.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 5)/          5.00000000/
      DATA cmatrx( 1, 5, 2, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 0)/        780.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 1)/        780.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 2)/        342.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 3)/         82.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 4)/         10.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 0)/        129.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 1)/        129.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 2)/         53.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 3)/         10.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 0)/         15.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 1)/         15.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 2)/          5.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 0)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 1)/          1.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 5, 2, 1, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 1600.

      DATA cmatrx( 1, 6, 0, 0, 0, 0)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 1)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 2)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 3)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 4)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 5)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 6)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 0, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 0)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 1)/     -20160.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 2)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 3)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 4)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 5)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 6)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 1, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 0)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 1)/      -9720.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 2)/      -4680.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 3)/      -1440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 4)/       -315.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 5)/        -51.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 0, 0, 2, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 0)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 1)/      -3000.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 2)/      -1440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 3)/       -440.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 4)/        -95.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 0)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 1)/       -660.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 2)/       -315.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 3)/        -95.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 0)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 1)/       -108.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 2)/        -51.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 0)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 1)/        -13.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 0, 0, 7, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1610.

      DATA cmatrx( 1, 6, 1, 0, 0, 0)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 1)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 2)/      42840.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 3)/      12600.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 4)/       2550.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 5)/        366.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 6)/         36.00000000/
      DATA cmatrx( 1, 6, 1, 0, 0, 7)/          2.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 0)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 1)/      90720.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 2)/      42840.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 3)/      12600.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 4)/       2550.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 5)/        366.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 6)/         36.00000000/
      DATA cmatrx( 1, 6, 1, 0, 1, 7)/          2.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 0)/      45360.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 1)/      45360.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 2)/      21420.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 3)/       6300.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 4)/       1275.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 5)/        183.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 6)/         18.00000000/
      DATA cmatrx( 1, 6, 1, 0, 2, 7)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 0)/      15120.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 1)/      15120.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 2)/       7140.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 3)/       2100.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 4)/        425.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 5)/         61.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 6)/          6.00000000/
      DATA cmatrx( 1, 6, 1, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 0)/       3750.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 1)/       3750.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 2)/       1770.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 3)/        520.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 4)/        105.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 5)/         15.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 0)/        726.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 1)/        726.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 2)/        342.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 3)/        100.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 4)/         20.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 0)/        111.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 1)/        111.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 2)/         52.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 3)/         15.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 0)/         13.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 1)/         13.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 2)/          6.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 0)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 1)/          1.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 0, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1611.

      DATA cmatrx( 1, 6, 1, 1, 0, 0)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 1)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 2)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 3)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 4)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 5)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 6)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 0, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 0)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 1)/     -45360.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 2)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 3)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 4)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 5)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 6)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 1, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 0)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 1)/     -21420.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 2)/     -10050.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 3)/      -2910.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 4)/       -570.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 5)/        -76.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 1, 1, 2, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 0)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 1)/      -6300.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 2)/      -2910.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 3)/       -810.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 4)/       -145.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 0)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 1)/      -1275.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 2)/       -570.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 3)/       -145.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 0)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 1)/       -183.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 2)/        -76.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 0)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 1)/        -18.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 1, 1, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1620.

      DATA cmatrx( 1, 6, 2, 0, 0, 0)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 1)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 2)/    -317520.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 3)/     -90720.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 4)/     -17415.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 5)/      -2295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 6)/       -198.00000000/
      DATA cmatrx( 1, 6, 2, 0, 0, 7)/         -9.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 0)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 1)/    -680400.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 2)/    -317520.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 3)/     -90720.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 4)/     -17415.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 5)/      -2295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 6)/       -198.00000000/
      DATA cmatrx( 1, 6, 2, 0, 1, 7)/         -9.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 0)/    -337680.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 1)/    -337680.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 2)/    -157410.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 3)/     -44850.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 4)/      -8562.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 5)/      -1116.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 6)/        -94.00000000/
      DATA cmatrx( 1, 6, 2, 0, 2, 7)/         -4.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 0)/    -110880.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 1)/    -110880.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 2)/     -51570.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 3)/     -14610.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 4)/      -2757.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 5)/       -351.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 6)/        -28.00000000/
      DATA cmatrx( 1, 6, 2, 0, 3, 7)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 0)/     -27135.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 1)/     -27135.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 2)/     -12582.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 3)/      -3537.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 4)/       -657.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 5)/        -81.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 6)/         -6.00000000/
      DATA cmatrx( 1, 6, 2, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 0)/      -5295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 1)/      -5295.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 2)/      -2448.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 3)/       -683.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 4)/       -125.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 5)/        -15.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 0)/       -858.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 1)/       -858.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 2)/       -396.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 3)/       -110.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 4)/        -20.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 0)/       -117.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 1)/       -117.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 2)/        -54.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 3)/        -15.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 0)/        -13.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 1)/        -13.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 2)/         -6.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 0)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 1)/         -1.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 0, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 1621.

      DATA cmatrx( 1, 6, 2, 1, 0, 0)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 1)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 2)/     105840.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 3)/      30240.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 4)/       5805.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 5)/        765.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 6)/         66.00000000/
      DATA cmatrx( 1, 6, 2, 1, 0, 7)/          3.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 0)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 1)/     226800.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 2)/     105840.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 3)/      30240.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 4)/       5805.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 5)/        765.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 6)/         66.00000000/
      DATA cmatrx( 1, 6, 2, 1, 1, 7)/          3.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 0)/     110880.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 1)/     110880.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 2)/      51570.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 3)/      14610.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 4)/       2757.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 5)/        351.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 6)/         28.00000000/
      DATA cmatrx( 1, 6, 2, 1, 2, 7)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 0)/      35280.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 1)/      35280.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 2)/      16290.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 3)/       4530.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 4)/        822.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 5)/         96.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 6)/          6.00000000/
      DATA cmatrx( 1, 6, 2, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 0)/       8145.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 1)/       8145.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 2)/       3708.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 3)/        993.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 4)/        165.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 5)/         15.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 0)/       1425.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 1)/       1425.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 2)/        630.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 3)/        155.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 4)/         20.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 0)/        189.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 1)/        189.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 2)/         78.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 3)/         15.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 0)/         18.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 1)/         18.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 2)/          6.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 0)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 1)/          1.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 1, 6, 2, 1, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2100.

      DATA cmatrx( 2, 1, 0, 0, 0, 0)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 1)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 2)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 0, 3)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 0)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 1)/         45.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 2)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 1, 3)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 0)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 1)/         18.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 2)/          7.00000000/
      DATA cmatrx( 2, 1, 0, 0, 2, 3)/          1.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 0)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 1)/          3.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 2)/          1.00000000/
      DATA cmatrx( 2, 1, 0, 0, 3, 3)/          0.00000000/

C Those are data for Lprime,N,L,M 2200.

      DATA cmatrx( 2, 2, 0, 0, 0, 0)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 1)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 2)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 3)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 0, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 0)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 1)/        270.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 2)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 3)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 1, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 0)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 1)/        117.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 2)/         50.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 3)/         11.00000000/
      DATA cmatrx( 2, 2, 0, 0, 2, 4)/          1.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 0)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 1)/         27.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 2)/         11.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 3)/          2.00000000/
      DATA cmatrx( 2, 2, 0, 0, 3, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 0)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 1)/          3.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 2)/          1.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 0, 0, 4, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2210.

      DATA cmatrx( 2, 2, 1, 0, 0, 0)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 1)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 2)/       -405.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 3)/        -90.00000000/
      DATA cmatrx( 2, 2, 1, 0, 0, 4)/         -9.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 0)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 1)/       -945.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 2)/       -405.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 3)/        -90.00000000/
      DATA cmatrx( 2, 2, 1, 0, 1, 4)/         -9.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 0)/       -450.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 1)/       -450.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 2)/       -192.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 3)/        -42.00000000/
      DATA cmatrx( 2, 2, 1, 0, 2, 4)/         -4.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 0)/       -135.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 1)/       -135.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 2)/        -57.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 3)/        -12.00000000/
      DATA cmatrx( 2, 2, 1, 0, 3, 4)/         -1.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 0)/        -27.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 1)/        -27.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 2)/        -11.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 3)/         -2.00000000/
      DATA cmatrx( 2, 2, 1, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 0)/         -3.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 1)/         -3.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 2)/         -1.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 0, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2211.

      DATA cmatrx( 2, 2, 1, 1, 0, 0)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 1)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 2)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 3)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 0, 4)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 0)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 1)/        945.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 2)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 3)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 1, 4)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 0)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 1)/        405.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 2)/        171.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 3)/         36.00000000/
      DATA cmatrx( 2, 2, 1, 1, 2, 4)/          3.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 0)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 1)/         90.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 2)/         36.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 3)/          6.00000000/
      DATA cmatrx( 2, 2, 1, 1, 3, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 0)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 1)/          9.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 2)/          3.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 0)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 1)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 2)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 2, 1, 1, 5, 4)/          0.00000000/

C Those are data for Lprime,N,L,M 2300.

      DATA cmatrx( 2, 3, 0, 0, 0, 0)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 1)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 2)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 3)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 4)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 0, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 0)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 1)/       1890.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 2)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 3)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 4)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 1, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 0)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 1)/        855.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 2)/        384.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 3)/         99.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 4)/         15.00000000/
      DATA cmatrx( 2, 3, 0, 0, 2, 5)/          1.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 0)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 1)/        225.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 2)/         99.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 3)/         24.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 4)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 0)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 1)/         36.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 2)/         15.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 3)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 0)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 1)/          3.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 2)/          1.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 0, 0, 5, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2310.

      DATA cmatrx( 2, 3, 1, 0, 0, 0)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 1)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 2)/      -3375.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 3)/       -855.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 4)/       -126.00000000/
      DATA cmatrx( 2, 3, 1, 0, 0, 5)/         -9.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 0)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 1)/      -7560.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 2)/      -3375.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 3)/       -855.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 4)/       -126.00000000/
      DATA cmatrx( 2, 3, 1, 0, 1, 5)/         -9.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 0)/      -3645.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 1)/      -3645.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 2)/      -1623.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 3)/       -408.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 4)/        -59.00000000/
      DATA cmatrx( 2, 3, 1, 0, 2, 5)/         -4.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 0)/      -1125.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 1)/      -1125.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 2)/       -498.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 3)/       -123.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 4)/        -17.00000000/
      DATA cmatrx( 2, 3, 1, 0, 3, 5)/         -1.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 0)/       -243.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 1)/       -243.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 2)/       -106.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 3)/        -25.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 4)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 0)/        -36.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 1)/        -36.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 2)/        -15.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 3)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 0)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 1)/         -3.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 2)/         -1.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 0, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2311.

      DATA cmatrx( 2, 3, 1, 1, 0, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 2)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 3)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 4)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 0, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 2)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 3)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 4)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 1, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 0)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 1)/       3375.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 2)/       1494.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 3)/        369.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 4)/         51.00000000/
      DATA cmatrx( 2, 3, 1, 1, 2, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 0)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 1)/        855.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 2)/        369.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 3)/         84.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 4)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 0)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 1)/        126.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 2)/         51.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 3)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 0)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 1)/          9.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 2)/          3.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 1, 1, 6, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2320.

      DATA cmatrx( 2, 3, 2, 0, 0, 0)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 1)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 2)/      22680.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 3)/       5670.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 4)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 0, 5)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 0)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 1)/      51030.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 2)/      22680.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 3)/       5670.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 4)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 1, 5)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 0)/      24570.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 1)/      24570.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 2)/      10890.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 3)/       2700.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 4)/        378.00000000/
      DATA cmatrx( 2, 3, 2, 0, 2, 5)/         24.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 0)/       7560.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 1)/       7560.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 2)/       3330.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 3)/        810.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 4)/        108.00000000/
      DATA cmatrx( 2, 3, 2, 0, 3, 5)/          6.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 0)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 1)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 2)/        726.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 3)/        171.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 4)/         21.00000000/
      DATA cmatrx( 2, 3, 2, 0, 4, 5)/          1.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 0)/        279.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 1)/        279.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 2)/        120.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 3)/         27.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 4)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 0)/         36.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 1)/         36.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 2)/         15.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 3)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 0)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 1)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 2)/          1.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 0, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2321.

      DATA cmatrx( 2, 3, 2, 1, 0, 0)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 1)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 2)/     -15120.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 3)/      -3780.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 4)/       -540.00000000/
      DATA cmatrx( 2, 3, 2, 1, 0, 5)/        -36.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 0)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 1)/     -34020.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 2)/     -15120.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 3)/      -3780.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 4)/       -540.00000000/
      DATA cmatrx( 2, 3, 2, 1, 1, 5)/        -36.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 0)/     -16065.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 1)/     -16065.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 2)/      -7110.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 3)/      -1755.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 4)/       -243.00000000/
      DATA cmatrx( 2, 3, 2, 1, 2, 5)/        -15.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 0)/      -4725.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 1)/      -4725.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 2)/      -2070.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 3)/       -495.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 4)/        -63.00000000/
      DATA cmatrx( 2, 3, 2, 1, 3, 5)/         -3.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 0)/       -945.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 1)/       -945.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 2)/       -405.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 3)/        -90.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 4)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 0)/       -126.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 1)/       -126.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 2)/        -51.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 3)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 0)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 1)/         -9.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 2)/         -3.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 1, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2322.

      DATA cmatrx( 2, 3, 2, 2, 0, 0)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 1)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 2)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 3)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 4)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 0, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 0)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 1)/       8505.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 2)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 3)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 4)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 1, 5)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 0)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 1)/       3780.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 2)/       1665.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 3)/        405.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 4)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 2, 2, 5)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 0)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 1)/        945.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 2)/        405.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 3)/         90.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 4)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 3, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 0)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 1)/        135.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 2)/         54.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 3)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 0)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 1)/          9.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 2)/          3.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 3, 2, 2, 7, 5)/          0.00000000/

C Those are data for Lprime,N,L,M 2400.

      DATA cmatrx( 2, 4, 0, 0, 0, 0)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 1)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 2)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 3)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 4)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 5)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 0, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 0)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 1)/      15120.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 2)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 3)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 4)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 5)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 1, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 0)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 1)/       7020.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 2)/       3246.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 3)/        906.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 4)/        165.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 5)/         19.00000000/
      DATA cmatrx( 2, 4, 0, 0, 2, 6)/          1.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 0)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 1)/       1980.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 2)/        906.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 3)/        246.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 4)/         42.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 5)/          4.00000000/
      DATA cmatrx( 2, 4, 0, 0, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 0)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 1)/        369.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 2)/        165.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 3)/         42.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 4)/          6.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 0)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 1)/         45.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 2)/         19.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 3)/          4.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 0)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 1)/          3.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 2)/          1.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 0, 0, 6, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2410.

      DATA cmatrx( 2, 4, 1, 0, 0, 0)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 1)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 2)/     -31185.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 3)/      -8505.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 4)/      -1485.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 5)/       -162.00000000/
      DATA cmatrx( 2, 4, 1, 0, 0, 6)/         -9.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 0)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 1)/     -68040.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 2)/     -31185.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 3)/      -8505.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 4)/      -1485.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 5)/       -162.00000000/
      DATA cmatrx( 2, 4, 1, 0, 1, 6)/         -9.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 0)/     -33075.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 1)/     -33075.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 2)/     -15135.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 3)/      -4110.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 4)/       -711.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 5)/        -76.00000000/
      DATA cmatrx( 2, 4, 1, 0, 2, 6)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 0)/     -10395.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 1)/     -10395.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 2)/      -4740.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 3)/      -1275.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 4)/       -216.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 5)/        -22.00000000/
      DATA cmatrx( 2, 4, 1, 0, 3, 6)/         -1.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 0)/      -2340.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 1)/      -2340.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 2)/      -1059.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 3)/       -279.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 4)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 5)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 0)/       -387.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 1)/       -387.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 2)/       -172.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 3)/        -43.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 4)/         -6.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 0)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 1)/        -45.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 2)/        -19.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 3)/         -4.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 0)/         -3.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 1)/         -3.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 2)/         -1.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 0, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2411.

      DATA cmatrx( 2, 4, 1, 1, 0, 0)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 1)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 2)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 3)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 4)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 5)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 0, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 0)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 1)/      68040.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 2)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 3)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 4)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 5)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 1, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 0)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 1)/      31185.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 2)/      14220.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 3)/       3825.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 4)/        648.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 5)/         66.00000000/
      DATA cmatrx( 2, 4, 1, 1, 2, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 0)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 1)/       8505.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 2)/       3825.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 3)/        990.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 4)/        153.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 5)/         12.00000000/
      DATA cmatrx( 2, 4, 1, 1, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 0)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 1)/       1485.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 2)/        648.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 3)/        153.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 4)/         18.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 0)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 1)/        162.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 2)/         66.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 3)/         12.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 0)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 1)/          9.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 2)/          3.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 1, 1, 7, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2420.

      DATA cmatrx( 2, 4, 2, 0, 0, 0)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 1)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 2)/     232470.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 3)/      62370.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 4)/      10530.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 5)/       1080.00000000/
      DATA cmatrx( 2, 4, 2, 0, 0, 6)/         54.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 0)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 1)/     510300.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 2)/     232470.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 3)/      62370.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 4)/      10530.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 5)/       1080.00000000/
      DATA cmatrx( 2, 4, 2, 0, 1, 6)/         54.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 0)/     247590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 1)/     247590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 2)/     112590.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 3)/      30060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 4)/       5022.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 5)/        504.00000000/
      DATA cmatrx( 2, 4, 2, 0, 2, 6)/         24.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 0)/      77490.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 1)/      77490.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 2)/      35100.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 3)/       9270.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 4)/       1512.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 5)/        144.00000000/
      DATA cmatrx( 2, 4, 2, 0, 3, 6)/          6.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 0)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 1)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 2)/       7899.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 3)/       2049.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 4)/        321.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 5)/         28.00000000/
      DATA cmatrx( 2, 4, 2, 0, 4, 6)/          1.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 0)/       3060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 1)/       3060.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 2)/       1365.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 3)/        345.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 4)/         51.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 5)/          4.00000000/
      DATA cmatrx( 2, 4, 2, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 0)/        423.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 1)/        423.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 2)/        186.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 3)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 4)/          6.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 0)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 1)/         45.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 2)/         19.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 3)/          4.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 0)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 1)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 2)/          1.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 0, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2421.

      DATA cmatrx( 2, 4, 2, 1, 0, 0)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 1)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 2)/    -154980.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 3)/     -41580.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 4)/      -7020.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 5)/       -720.00000000/
      DATA cmatrx( 2, 4, 2, 1, 0, 6)/        -36.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 0)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 1)/    -340200.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 2)/    -154980.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 3)/     -41580.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 4)/      -7020.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 5)/       -720.00000000/
      DATA cmatrx( 2, 4, 2, 1, 1, 6)/        -36.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 0)/    -162540.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 1)/    -162540.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 2)/     -73845.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 3)/     -19665.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 4)/      -3267.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 5)/       -324.00000000/
      DATA cmatrx( 2, 4, 2, 1, 2, 6)/        -15.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 0)/     -49140.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 1)/     -49140.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 2)/     -22185.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 3)/      -5805.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 4)/       -927.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 5)/        -84.00000000/
      DATA cmatrx( 2, 4, 2, 1, 3, 6)/         -3.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 0)/     -10395.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 1)/     -10395.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 2)/      -4635.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 3)/      -1170.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 4)/       -171.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 5)/        -12.00000000/
      DATA cmatrx( 2, 4, 2, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 0)/      -1575.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 1)/      -1575.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 2)/       -684.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 3)/       -159.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 4)/        -18.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 0)/       -162.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 1)/       -162.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 2)/        -66.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 3)/        -12.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 0)/         -9.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 1)/         -9.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 2)/         -3.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 1, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2422.

      DATA cmatrx( 2, 4, 2, 2, 0, 0)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 1)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 2)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 3)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 4)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 5)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 0, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 0)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 1)/      85050.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 2)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 3)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 4)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 5)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 1, 6)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 0)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 1)/      38745.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 2)/      17550.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 3)/       4635.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 4)/        756.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 5)/         72.00000000/
      DATA cmatrx( 2, 4, 2, 2, 2, 6)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 0)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 1)/      10395.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 2)/       4635.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 3)/       1170.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 4)/        171.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 5)/         12.00000000/
      DATA cmatrx( 2, 4, 2, 2, 3, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 0)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 1)/       1755.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 2)/        756.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 3)/        171.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 4)/         18.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 0)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 1)/        180.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 2)/         72.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 3)/         12.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 0)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 1)/          9.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 2)/          3.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 4, 2, 2, 8, 6)/          0.00000000/

C Those are data for Lprime,N,L,M 2500.

      DATA cmatrx( 2, 5, 0, 0, 0, 0)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 1)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 2)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 3)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 4)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 5)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 6)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 0, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 0)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 1)/     136080.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 2)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 3)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 4)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 5)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 6)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 1, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 0)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 1)/      64260.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 2)/      30270.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 3)/       8850.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 4)/       1770.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 5)/        248.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 6)/         23.00000000/
      DATA cmatrx( 2, 5, 0, 0, 2, 7)/          1.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 0)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 1)/      18900.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 2)/       8850.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 3)/       2550.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 4)/        495.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 5)/         65.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 6)/          5.00000000/
      DATA cmatrx( 2, 5, 0, 0, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 0)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 1)/       3825.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 2)/       1770.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 3)/        495.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 4)/         90.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 5)/         10.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 0)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 1)/        549.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 2)/        248.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 3)/         65.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 4)/         10.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 0)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 1)/         54.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 2)/         23.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 3)/          5.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 0)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 1)/          3.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 2)/          1.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 0, 0, 7, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2510.

      DATA cmatrx( 2, 5, 1, 0, 0, 0)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 1)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 2)/    -317520.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 3)/     -90720.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 4)/     -17415.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 5)/      -2295.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 6)/       -198.00000000/
      DATA cmatrx( 2, 5, 1, 0, 0, 7)/         -9.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 0)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 1)/    -680400.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 2)/    -317520.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 3)/     -90720.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 4)/     -17415.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 5)/      -2295.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 6)/       -198.00000000/
      DATA cmatrx( 2, 5, 1, 0, 1, 7)/         -9.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 0)/    -332640.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 1)/    -332640.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 2)/    -155070.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 3)/     -44190.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 4)/      -8439.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 5)/      -1101.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 6)/        -93.00000000/
      DATA cmatrx( 2, 5, 1, 0, 2, 7)/         -4.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 0)/    -105840.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 1)/    -105840.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 2)/     -49230.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 3)/     -13950.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 4)/      -2634.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 5)/       -336.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 6)/        -27.00000000/
      DATA cmatrx( 2, 5, 1, 0, 3, 7)/         -1.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 0)/     -24435.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 1)/     -24435.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 2)/     -11316.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 3)/      -3171.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 4)/       -585.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 5)/        -71.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 6)/         -5.00000000/
      DATA cmatrx( 2, 5, 1, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 0)/      -4275.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 1)/      -4275.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 2)/      -1962.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 3)/       -537.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 4)/        -94.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 5)/        -10.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 0)/       -567.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 1)/       -567.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 2)/       -255.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 3)/        -66.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 4)/        -10.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 0)/        -54.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 1)/        -54.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 2)/        -23.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 3)/         -5.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 0)/         -3.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 1)/         -3.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 2)/         -1.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 0, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2511.

      DATA cmatrx( 2, 5, 1, 1, 0, 0)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 1)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 2)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 3)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 4)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 5)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 6)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 0, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 0)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 1)/     680400.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 2)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 3)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 4)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 5)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 6)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 1, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 0)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 1)/     317520.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 2)/     147690.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 3)/      41850.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 4)/       7902.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 5)/       1008.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 6)/         81.00000000/
      DATA cmatrx( 2, 5, 1, 1, 2, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 0)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 1)/      90720.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 2)/      41850.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 3)/      11610.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 4)/       2097.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 5)/        243.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 6)/         15.00000000/
      DATA cmatrx( 2, 5, 1, 1, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 0)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 1)/      17415.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 2)/       7902.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 3)/       2097.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 4)/        342.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 5)/         30.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 0)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 1)/       2295.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 2)/       1008.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 3)/        243.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 4)/         30.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 0)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 1)/        198.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 2)/         81.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 3)/         15.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 0)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 1)/          9.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 2)/          3.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 1, 1, 8, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2520.

      DATA cmatrx( 2, 5, 2, 0, 0, 0)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 1)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 2)/    2602530.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 3)/     731430.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 4)/     136080.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 5)/      17010.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 6)/       1350.00000000/
      DATA cmatrx( 2, 5, 2, 0, 0, 7)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 0)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 1)/    5613300.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 2)/    2602530.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 3)/     731430.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 4)/     136080.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 5)/      17010.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 6)/       1350.00000000/
      DATA cmatrx( 2, 5, 2, 0, 1, 7)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 0)/    2738610.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 1)/    2738610.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 2)/    1268190.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 3)/     355320.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 4)/      65700.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 5)/       8118.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 6)/        630.00000000/
      DATA cmatrx( 2, 5, 2, 0, 2, 7)/         24.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 0)/     867510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 1)/     867510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 2)/     400680.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 3)/     111510.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 4)/      20340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 5)/       2448.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 6)/        180.00000000/
      DATA cmatrx( 2, 5, 2, 0, 3, 7)/          6.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 0)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 1)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 2)/      92145.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 3)/      25365.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 4)/       4524.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 5)/        521.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 6)/         35.00000000/
      DATA cmatrx( 2, 5, 2, 0, 4, 7)/          1.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 0)/      35910.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 1)/      35910.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 2)/      16419.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 3)/       4449.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 4)/        768.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 5)/         83.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 6)/          5.00000000/
      DATA cmatrx( 2, 5, 2, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 0)/       5175.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 1)/       5175.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 2)/       2346.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 3)/        621.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 4)/        102.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 5)/         10.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 0)/        603.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 1)/        603.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 2)/        269.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 3)/         68.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 4)/         10.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 0)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 1)/         54.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 2)/         23.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 3)/          5.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 0)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 1)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 2)/          1.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 0, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2521.

      DATA cmatrx( 2, 5, 2, 1, 0, 0)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 0, 1)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 0, 2)/   -1735019.99999999/
      DATA cmatrx( 2, 5, 2, 1, 0, 3)/    -487619.99999999/
      DATA cmatrx( 2, 5, 2, 1, 0, 4)/     -90720.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 5)/     -11340.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 6)/       -900.00000000/
      DATA cmatrx( 2, 5, 2, 1, 0, 7)/        -36.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 0)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 1, 1)/   -3742199.99999997/
      DATA cmatrx( 2, 5, 2, 1, 1, 2)/   -1735019.99999998/
      DATA cmatrx( 2, 5, 2, 1, 1, 3)/    -487619.99999999/
      DATA cmatrx( 2, 5, 2, 1, 1, 4)/     -90720.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 5)/     -11340.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 6)/       -900.00000000/
      DATA cmatrx( 2, 5, 2, 1, 1, 7)/        -36.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 0)/   -1803059.99999999/
      DATA cmatrx( 2, 5, 2, 1, 2, 1)/   -1803059.99999998/
      DATA cmatrx( 2, 5, 2, 1, 2, 2)/    -834434.99999999/
      DATA cmatrx( 2, 5, 2, 1, 2, 3)/    -233415.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 4)/     -43020.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 5)/      -5283.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 6)/       -405.00000000/
      DATA cmatrx( 2, 5, 2, 1, 2, 7)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 0)/    -555659.99999999/
      DATA cmatrx( 2, 5, 2, 1, 3, 1)/    -555659.99999999/
      DATA cmatrx( 2, 5, 2, 1, 3, 2)/    -256095.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 3)/     -70875.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 4)/     -12780.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 5)/      -1503.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 6)/       -105.00000000/
      DATA cmatrx( 2, 5, 2, 1, 3, 7)/         -3.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 0)/    -121905.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 1)/    -121905.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 2)/     -55755.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 3)/     -15120.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 4)/      -2610.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 5)/       -279.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 6)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 0)/     -19845.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 1)/     -19845.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 2)/      -8946.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 3)/      -2331.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 4)/       -366.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 5)/        -30.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 0)/      -2385.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 1)/      -2385.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 2)/      -1044.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 3)/       -249.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 4)/        -30.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 0)/       -198.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 1)/       -198.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 2)/        -81.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 3)/        -15.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 0)/         -9.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 1)/         -9.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 2)/         -3.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 1, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2522.

      DATA cmatrx( 2, 5, 2, 2, 0, 0)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 0, 1)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 0, 2)/     433755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 3)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 4)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 5)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 6)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 0, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 0)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 1)/     935549.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 2)/     433754.99999999/
      DATA cmatrx( 2, 5, 2, 2, 1, 3)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 4)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 5)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 6)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 1, 7)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 0)/     433755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 1)/     433754.99999999/
      DATA cmatrx( 2, 5, 2, 2, 2, 2)/     200340.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 3)/      55755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 4)/      10170.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 5)/       1224.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 6)/         90.00000000/
      DATA cmatrx( 2, 5, 2, 2, 2, 7)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 0)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 1)/     121905.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 2)/      55755.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 3)/      15120.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 4)/       2610.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 5)/        279.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 6)/         15.00000000/
      DATA cmatrx( 2, 5, 2, 2, 3, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 0)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 1)/      22680.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 2)/      10170.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 3)/       2610.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 4)/        396.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 5)/         30.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 0)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 1)/       2835.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 2)/       1224.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 3)/        279.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 4)/         30.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 0)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 1)/        225.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 2)/         90.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 3)/         15.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 0)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 1)/          9.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 2)/          3.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 5, 2, 2, 9, 7)/          0.00000000/

C Those are data for Lprime,N,L,M 2600.

      DATA cmatrx( 2, 6, 0, 0, 0, 0)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 1)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 2)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 3)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 4)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 5)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 6)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 7)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 0, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 0)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 1)/    1360800.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 2)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 3)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 4)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 5)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 6)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 7)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 1, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 0)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 1)/     650160.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 2)/     310140.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 3)/      93420.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 4)/      19755.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 5)/       3063.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 6)/        348.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 7)/         27.00000000/
      DATA cmatrx( 2, 6, 0, 0, 2, 8)/          1.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 0)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 1)/     196560.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 2)/      93420.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 3)/      27900.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 4)/       5805.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 5)/        873.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 6)/         93.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 7)/          6.00000000/
      DATA cmatrx( 2, 6, 0, 0, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 0)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 1)/      41850.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 2)/      19755.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 3)/       5805.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 4)/       1170.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 5)/        165.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 6)/         15.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 0)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 1)/       6570.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 2)/       3063.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 3)/        873.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 4)/        165.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 5)/         20.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 0)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 1)/        765.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 2)/        348.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 3)/         93.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 4)/         15.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 0)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 1)/         63.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 2)/         27.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 3)/          6.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 0)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 1)/          3.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 2)/          1.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 0, 0, 8, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2610.

      DATA cmatrx( 2, 6, 1, 0, 0, 0)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 1)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 2)/   -3538080.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 3)/   -1043279.99999999/
      DATA cmatrx( 2, 6, 1, 0, 0, 4)/    -212625.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 5)/     -31185.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 6)/      -3285.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 7)/       -234.00000000/
      DATA cmatrx( 2, 6, 1, 0, 0, 8)/         -9.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 0)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 1)/   -7484400.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 2)/   -3538080.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 3)/   -1043279.99999999/
      DATA cmatrx( 2, 6, 1, 0, 1, 4)/    -212625.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 5)/     -31185.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 6)/      -3285.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 7)/       -234.00000000/
      DATA cmatrx( 2, 6, 1, 0, 1, 8)/         -9.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 0)/   -3674160.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 1)/   -3674160.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 2)/   -1735650.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 3)/    -510930.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 4)/    -103800.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 5)/     -15138.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 6)/      -1578.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 7)/       -110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 2, 8)/         -4.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 0)/   -1179359.99999999/
      DATA cmatrx( 2, 6, 1, 0, 3, 1)/   -1179359.99999999/
      DATA cmatrx( 2, 6, 1, 0, 3, 2)/    -556290.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 3)/    -163170.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 4)/     -32925.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 5)/      -4743.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 6)/       -483.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 7)/        -32.00000000/
      DATA cmatrx( 2, 6, 1, 0, 3, 8)/         -1.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 0)/    -276885.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 1)/    -276885.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 2)/    -130245.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 3)/     -37950.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 4)/      -7560.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 5)/      -1063.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 6)/       -103.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 7)/         -6.00000000/
      DATA cmatrx( 2, 6, 1, 0, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 0)/     -50085.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 1)/     -50085.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 2)/     -23439.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 3)/      -6744.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 4)/      -1310.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 5)/       -175.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 6)/        -15.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 0)/      -7110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 1)/      -7110.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 2)/      -3294.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 3)/       -924.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 4)/       -170.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 5)/        -20.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 0)/       -783.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 1)/       -783.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 2)/       -355.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 3)/        -94.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 4)/        -15.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 0)/        -63.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 1)/        -63.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 2)/        -27.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 3)/         -6.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 0)/         -3.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 1)/         -3.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 2)/         -1.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 0, 9, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2611.

      DATA cmatrx( 2, 6, 1, 1, 0, 0)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 1)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 2)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 3)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 0, 4)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 5)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 6)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 7)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 0, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 0)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 1)/    7484400.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 2)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 3)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 1, 4)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 5)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 6)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 7)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 1, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 0)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 1)/    3538080.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 2)/    1668870.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 3)/     489510.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 4)/      98775.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 5)/      14229.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 6)/       1449.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 7)/         96.00000000/
      DATA cmatrx( 2, 6, 1, 1, 2, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 0)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 3, 1)/    1043279.99999999/
      DATA cmatrx( 2, 6, 1, 1, 3, 2)/     489510.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 3)/     141750.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 4)/      27900.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 5)/       3834.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 6)/        354.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 7)/         18.00000000/
      DATA cmatrx( 2, 6, 1, 1, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 0)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 1)/     212625.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 2)/      98775.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 3)/      27900.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 4)/       5220.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 5)/        645.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 6)/         45.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 0)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 1)/      31185.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 2)/      14229.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 3)/       3834.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 4)/        645.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 5)/         60.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 0)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 1)/       3285.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 2)/       1449.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 3)/        354.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 4)/         45.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 0)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 1)/        234.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 2)/         96.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 3)/         18.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 0)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 1)/          9.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 2)/          3.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 1, 1, 9, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2620.

      DATA cmatrx( 2, 6, 2, 0, 0, 0)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 1)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 2)/   31638600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 3)/    9185399.99999988/
      DATA cmatrx( 2, 6, 2, 0, 0, 4)/    1820069.99999997/
      DATA cmatrx( 2, 6, 2, 0, 0, 5)/     255149.99999999/
      DATA cmatrx( 2, 6, 2, 0, 0, 6)/      25110.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 7)/       1620.00000000/
      DATA cmatrx( 2, 6, 2, 0, 0, 8)/         54.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 0)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 1)/   67359600.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 2)/   31638599.99999970/
      DATA cmatrx( 2, 6, 2, 0, 1, 3)/    9185399.99999988/
      DATA cmatrx( 2, 6, 2, 0, 1, 4)/    1820069.99999997/
      DATA cmatrx( 2, 6, 2, 0, 1, 5)/     255150.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 6)/      25110.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 7)/       1620.00000000/
      DATA cmatrx( 2, 6, 2, 0, 1, 8)/         54.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 0)/   32999400.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 1)/   32999399.99999970/
      DATA cmatrx( 2, 6, 2, 0, 2, 2)/   15486659.99999980/
      DATA cmatrx( 2, 6, 2, 0, 2, 3)/    4486859.99999994/
      DATA cmatrx( 2, 6, 2, 0, 2, 4)/     885599.99999999/
      DATA cmatrx( 2, 6, 2, 0, 2, 5)/     123300.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 6)/      11988.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 7)/        756.00000000/
      DATA cmatrx( 2, 6, 2, 0, 2, 8)/         24.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 0)/   10546199.99999980/
      DATA cmatrx( 2, 6, 2, 0, 3, 1)/   10546199.99999980/
      DATA cmatrx( 2, 6, 2, 0, 3, 2)/    4940459.99999994/
      DATA cmatrx( 2, 6, 2, 0, 3, 3)/    1425059.99999998/
      DATA cmatrx( 2, 6, 2, 0, 3, 4)/     278909.99999999/
      DATA cmatrx( 2, 6, 2, 0, 3, 5)/      38250.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 6)/       3618.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 7)/        216.00000000/
      DATA cmatrx( 2, 6, 2, 0, 3, 8)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 0)/    2470229.99999997/
      DATA cmatrx( 2, 6, 2, 0, 4, 1)/    2470229.99999997/
      DATA cmatrx( 2, 6, 2, 0, 4, 2)/    1153889.99999999/
      DATA cmatrx( 2, 6, 2, 0, 4, 3)/     330479.99999999/
      DATA cmatrx( 2, 6, 2, 0, 4, 4)/      63801.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 5)/       8535.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 6)/        771.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 7)/         42.00000000/
      DATA cmatrx( 2, 6, 2, 0, 4, 8)/          1.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 0)/     451709.99999999/
      DATA cmatrx( 2, 6, 2, 0, 5, 1)/     451710.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 2)/     210150.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 3)/      59580.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 4)/      11277.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 5)/       1455.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 6)/        123.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 7)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 0)/      66960.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 1)/      66960.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 2)/      30978.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 3)/       8658.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 4)/       1593.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 5)/        195.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 6)/         15.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 0)/       8190.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 1)/       8190.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 2)/       3756.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 3)/       1026.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 4)/        180.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 5)/         20.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 0)/        819.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 1)/        819.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 2)/        369.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 3)/         96.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 4)/         15.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 0)/         63.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 1)/         63.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 2)/         27.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 3)/          6.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 0)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 1)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 2)/          1.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 0,10, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2621.

      DATA cmatrx( 2, 6, 2, 1, 0, 0)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 0, 1)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 0, 2)/  -21092399.99999970/
      DATA cmatrx( 2, 6, 2, 1, 0, 3)/   -6123599.99999985/
      DATA cmatrx( 2, 6, 2, 1, 0, 4)/   -1213379.99999996/
      DATA cmatrx( 2, 6, 2, 1, 0, 5)/    -170099.99999999/
      DATA cmatrx( 2, 6, 2, 1, 0, 6)/     -16740.00000000/
      DATA cmatrx( 2, 6, 2, 1, 0, 7)/      -1080.00000000/
      DATA cmatrx( 2, 6, 2, 1, 0, 8)/        -36.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 0)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 1, 1)/  -44906399.99999950/
      DATA cmatrx( 2, 6, 2, 1, 1, 2)/  -21092399.99999960/
      DATA cmatrx( 2, 6, 2, 1, 1, 3)/   -6123599.99999985/
      DATA cmatrx( 2, 6, 2, 1, 1, 4)/   -1213379.99999998/
      DATA cmatrx( 2, 6, 2, 1, 1, 5)/    -170099.99999999/
      DATA cmatrx( 2, 6, 2, 1, 1, 6)/     -16740.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 7)/      -1080.00000000/
      DATA cmatrx( 2, 6, 2, 1, 1, 8)/        -36.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 0)/  -21772799.99999970/
      DATA cmatrx( 2, 6, 2, 1, 2, 1)/  -21772799.99999960/
      DATA cmatrx( 2, 6, 2, 1, 2, 2)/  -10213559.99999980/
      DATA cmatrx( 2, 6, 2, 1, 2, 3)/   -2955959.99999994/
      DATA cmatrx( 2, 6, 2, 1, 2, 4)/    -582254.99999998/
      DATA cmatrx( 2, 6, 2, 1, 2, 5)/     -80775.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 6)/      -7803.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 7)/       -486.00000000/
      DATA cmatrx( 2, 6, 2, 1, 2, 8)/        -15.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 0)/   -6803999.99999985/
      DATA cmatrx( 2, 6, 2, 1, 3, 1)/   -6803999.99999985/
      DATA cmatrx( 2, 6, 2, 1, 3, 2)/   -3182759.99999991/
      DATA cmatrx( 2, 6, 2, 1, 3, 3)/    -914759.99999997/
      DATA cmatrx( 2, 6, 2, 1, 3, 4)/    -177794.99999999/
      DATA cmatrx( 2, 6, 2, 1, 3, 5)/     -24075.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 6)/      -2223.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 7)/       -126.00000000/
      DATA cmatrx( 2, 6, 2, 1, 3, 8)/         -3.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 0)/   -1530899.99999996/
      DATA cmatrx( 2, 6, 2, 1, 4, 1)/   -1530899.99999996/
      DATA cmatrx( 2, 6, 2, 1, 4, 2)/    -712529.99999998/
      DATA cmatrx( 2, 6, 2, 1, 4, 3)/    -202229.99999999/
      DATA cmatrx( 2, 6, 2, 1, 4, 4)/     -38340.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 5)/      -4950.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 6)/       -414.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 7)/        -18.00000000/
      DATA cmatrx( 2, 6, 2, 1, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 0)/    -260819.99999999/
      DATA cmatrx( 2, 6, 2, 1, 5, 1)/    -260819.99999999/
      DATA cmatrx( 2, 6, 2, 1, 5, 2)/    -120330.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 3)/     -33390.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 4)/      -6039.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 5)/       -705.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 6)/        -45.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 0)/     -34155.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 1)/     -34155.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 2)/     -15507.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 3)/      -4122.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 4)/       -675.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 5)/        -60.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 0)/      -3375.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 1)/      -3375.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 2)/      -1485.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 3)/       -360.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 4)/        -45.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 0)/       -234.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 1)/       -234.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 2)/        -96.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 3)/        -18.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 0)/         -9.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 1)/         -9.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 2)/         -3.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 1,10, 8)/          0.00000000/

C Those are data for Lprime,N,L,M 2622.

      DATA cmatrx( 2, 6, 2, 2, 0, 0)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 0, 1)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 0, 2)/    5273099.99999994/
      DATA cmatrx( 2, 6, 2, 2, 0, 3)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 0, 4)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 0, 5)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 6)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 7)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 0, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 0)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 1, 1)/   11226599.99999980/
      DATA cmatrx( 2, 6, 2, 2, 1, 2)/    5273099.99999991/
      DATA cmatrx( 2, 6, 2, 2, 1, 3)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 1, 4)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 1, 5)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 6)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 7)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 1, 8)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 0)/    5273099.99999994/
      DATA cmatrx( 2, 6, 2, 2, 2, 1)/    5273099.99999991/
      DATA cmatrx( 2, 6, 2, 2, 2, 2)/    2470229.99999996/
      DATA cmatrx( 2, 6, 2, 2, 2, 3)/     712529.99999998/
      DATA cmatrx( 2, 6, 2, 2, 2, 4)/     139455.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 5)/      19125.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 6)/       1809.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 7)/        108.00000000/
      DATA cmatrx( 2, 6, 2, 2, 2, 8)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 0)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 3, 1)/    1530899.99999996/
      DATA cmatrx( 2, 6, 2, 2, 3, 2)/     712529.99999998/
      DATA cmatrx( 2, 6, 2, 2, 3, 3)/     202229.99999999/
      DATA cmatrx( 2, 6, 2, 2, 3, 4)/      38340.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 5)/       4950.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 6)/        414.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 7)/         18.00000000/
      DATA cmatrx( 2, 6, 2, 2, 3, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 0)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 4, 1)/     303344.99999999/
      DATA cmatrx( 2, 6, 2, 2, 4, 2)/     139455.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 3)/      38340.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 4)/       6804.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 5)/        765.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 6)/         45.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 4, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 0)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 1)/      42525.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 2)/      19125.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 3)/       4950.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 4)/        765.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 5)/         60.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 5, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 0)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 1)/       4185.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 2)/       1809.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 3)/        414.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 4)/         45.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 6, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 0)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 1)/        270.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 2)/        108.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 3)/         18.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 7, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 0)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 1)/          9.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 2)/          3.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 8, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2, 9, 8)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 0)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 1)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 2)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 3)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 4)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 5)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 6)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 7)/          0.00000000/
      DATA cmatrx( 2, 6, 2, 2,10, 8)/          0.00000000/

      END
