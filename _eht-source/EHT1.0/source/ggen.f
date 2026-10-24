C Note: The following section gives a fortran program to compute
C coefficients by removing "C "(C and a space).

C This program is made to obtain the C matrix  for a given l,N,L,M.
C The formula for C matrix is given in reference:
C H. W. Jones and C.A. Weatherford, Int. J. Quantum Chem. Symp. 12, 483 (1978).
C In Page 486 they wrote: We build up a C matrix Cl(NA,NR) using eight nested 
C DO loops to sort and collect the sums of all coefficients of all like powers
C of r and a.
C
C This program only gives all C matrix coefficients, the modulus of 
C N-L+2p+2q+2qprime-k-kprime to 2, the exponential index of a, and that of r.
C After this program is compiled with f77, you can obtain a series of data
C for C matrix.

      program gen

      implicit none
      integer sl,N,L,M,p,q,v,pp,qp,vp,t,k,kp
      integer result2,result3,result4,tmp
      double precision result,alpha,beta,efunc
       
      print *, "Please input l, N, L, M"
      read *,sl,N,L,M
      if (sl .lt. M) stop "Irrational pair: l and M"
      print 98,sl,N,L,M
      do 1 p=0,(L+M)/2
        do 2 q=0,(L+M-2*p)
          do 3 v=0,(L+M-2*p-q)
            do 4 pp=0,((sl-M)/2)
              do 5 qp=0,sl-M-2*pp
                do 6 vp=0,sl-M-2*pp-qp
                  t=N-L+2*p+2*q+2*qp
                  do 7 k=0,t
                    do 8 kp=0,t-k
                      result=(-1)**(L+v+qp)
                      result=result/4.0**(L+sl-p-pp)
                      tmp=L-M
                      result=result*alpha(L,tmp,p)
                      tmp=L+M-2*p
                      result=result*beta(tmp,q)
                      tmp=L+M-2*p-q
                      result=result*beta(tmp,v)
                      tmp=sl+M
                      result=result*alpha(sl,tmp,pp)
                      tmp=sl-M-2*pp
                      result=result*beta(tmp,qp)
                      tmp=sl-M-2*pp-qp
                      result=result*beta(tmp,vp)
                      tmp=N-L+2*p+2*q+2*qp
                      result=result*efunc(tmp,k)
                      tmp=N-L+2*p+2*q+2*qp-k
                      result=result*beta(tmp,kp)
                 
                      result2=N-L+2*p+2*q+2*qp-k-kp
                      result2=MOD(result2,2)

                      result3=N+L+2*sl-2*pp-2*vp-2*v-k-kp

                      result4=2*pp+2*v+2*vp+kp
                     
                      print 99, result,result2,result3,result4

 8                  continue
 7                continue
 6              continue
 5            continue
 4          continue
 3        continue
 2      continue
 1    continue

 98   format(5x,"tttt ",4(I4,2x))
 99   format(5x,f20.8,3(I4,2x))

      end
       

      double precision function alpha(sl,n,p)
      implicit none
      integer sl,n,p,tmpa,tmpb,tmpc
      double precision fact(0:14)
      common /factorial/ fact

      tmpa=2*sl-2*p
      tmpb=2*sl-2*p-n
      tmpc=sl-p

      alpha=(-1)**p*fact(tmpa)/(fact(tmpb)*fact(tmpc)*fact(p))

      return
      end

      double precision function beta(n,k)
      implicit none
      integer n,k,tmpa
      double precision fact(0:14)
      common /factorial/ fact

      tmpa=n-k

      beta=fact(n)/(fact(tmpa)*fact(k))

      return
      end

      double precision function efunc(t,k)
      implicit none
      integer t,k,tmpa
      double precision fact(0:14)
      common /factorial/ fact

      tmpa=t-k

      efunc=fact(t)/fact(tmpa)

      return
      end

      BLOCK DATA
      double precision fact(0:14)
      common /factorial/ fact
      data fact/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
     &          362880.0,3628800.0,39916800.0,479001600.0,622020800.0,
     &          87178291200.0/
      end

C Suppose the executable file of this program is a.out and you are in a
C UNIX system, the following commands let you sort and collect the sums of 
C all coefficients of all like powers C of r and a.
C     a.out | sort +2 -4 | awk -f awk.prog 
C awk.prog is a file to do some special job for unix command awk.
C File awk.prog begins.
C /ttt/ { sl=$2;n=$3;l=$4;m=$5}
C {pop[$3,$4] += $1}
C END { printf("\nC Those are data for Lprime,N,L,M %d%d%d%d.\n\n",sl,n,l,m);
C   for (i=0;i<=(sl+n+l);i++) 
C     for (j=0;j<=(sl+n);j++)  
C       printf("      DATA cmatrx(%2d,%2d,%2d,%2d,%2d,%2d)/%20.8f/\n",sl,n,l,
C         m,i,j,pop[i,j]);
C }
C File awk.prog ends.
C This is how the all above data are produced.
