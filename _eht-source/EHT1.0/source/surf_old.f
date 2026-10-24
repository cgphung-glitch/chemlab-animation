C 888888888888888888888888888888888888888888888888888888888888888888
c The subroutine receives two arguments (coordinate and the minumum 
c number of atoms the program can work with) and returns two values
c (the total energy and the first derivative for each atom computed 
c numerically.
c 888888888888888888888888888888888888888888888888888888888888888888
c
       SUBROUTINE SURF(V,X,DX,N3TM)
C      program tb
C
C This subroutine can be a subroutine of POLYRATE Program.
c
      IMPLICIT NONE
C
c V energy of the current configuration, in atomic unit
c X an array of Cartesian coordinates of the nuclei
c DX first derivatives of energy with respect to the nuclear position
c N3tm=3*NATOMS --- maximum number of atoms in POLYRATE code
C
      integer n3tm
c
      integer maxatm,max3,max9,nchoice
c
c maxatm --- maximum number of atoms (Could be set as N3TM/3) 
c max3   --- for the dimension of Cartesian coordinates or derivatives
c max9   --- for the dimension of Hamiltonian matrix
c
      parameter (maxatm=60,max3=3*maxatm,max9=9*maxatm)
c
      double precision V,X(max3),DX(max3)
c
c symb --- elememt symbol
c stemp --- tempary variable to hold a symbol
c
      character*2 symb(maxatm),stemp
c
c The following variables are loop counters.
c
      integer i,j,k,ii,jj,ip,jp
c
c natom --- The number of atoms in the current system.
c Natmic --- The temporary atomic number
c idtype --- An array to hold the atomic number for each atom.
c
      integer natom,Natmic,idtype(maxatm)
c
c vsip --- valence state ionilization potentials
c
      double precision vsip(103,3)
c
c ndim --- The dimension of the Hamiltonian matrix.
c nrdim, ncdim --- temporary dimension number for row and column.
c     
      integer ndim,nrdim,ncdim
c
c nElect --- The number of electrons in the system.
c nValen --- The number of valence electrons for an atom.
c
      integer nElect,nValen(103)
c
c nOrbit --- an array holding the number of orbitals for an atom.
c
      integer nOrbit(103)
c
c nsd --- The principal quantum number of s orbital
c npd --- The principal quantum number of p orbital
c ndd --- The principal quantum number of d orbital
c
      integer nsd(103),npd(103),ndd(103)
c
c hamilt --- an array holding all the elements of Hamiltonian matrix.
c hamil  --- an array holding a submatrix reslting from two elements.
c dhamlt --- an array holding all the elements of first derivatives of hamilt.
c dhaml  --- an array holding its first derivative.
c     St --- an array holding the overlap integrals
c     S  --- an array holding submatrix resulting from two atoms
c
       double precision hamilt(max9,max9)  , hamil(9,9)
       double precision dhamlt(max9,max9,3),dhamil(9,9,3)
       double precision     St(max9,max9)  ,     S(9,9)
c
c Some parameters for diagonalization of Hamiltonian matrix.
c
      integer info
      double precision ap(max9*(max9+1)/2)
      double precision bp(max9*(max9+1)/2)
      double precision work(3*max9)
c
c  eva --- an array to hold eigenvalues resulting from HC=ESC.
c evec --- an array to hold eigenvectors resulting from HC=ESC.
c
      double precision eva(max9),evec(max9,max9)
c
c  occ --- an array holding occupation number for each molecular orbital.
c
      double precision occ(maxatm)
c
c
c q --- an array holding net AO population due to one electron.
c p --- an array holding Mulliken overlap population.
c qcharg --- a array holding gross AO population on one atom.
c
      double precision q(max9),p(max9,max9),qcharg(maxatm)
c
c forcex,forcey,forcez --- temporary variables computing first derivatives.
c tmp --- a temporary variable to calculate first derivatives.
c 
      double precision forcex,forcey,forcez,tmp
c 
c first --- logical variable representing whether or not the symbols are read.
c
      logical first
c
c Common block "firsRD" will be saved after finishing reading.
c Common block "coulom" will share with subroutine "ehuckl,data block"
c Common block "nspdd" will share with subroutine "ehuckl,data block"
c Common block "elctrn" will share with subroutine "data block"
c
      common /firsRD/ natom,ndim,nElect,idtype,first
      common /coulom/ vsip
      common /nspdd/ nsd,npd,ndd
      common /elctrn/ nValen 
      common /orbit/ nOrbit 
      common /choice/ nchoice
c
      SAVE /firsRD/
c
      if (first) then
        ndim=0
        nElect=0
c
c You should provide a file "fort.99" to list all the atoms in the system.
c
        do 20 i=1, maxatm
          read(99,*,end=21) symb(i)
          stemp=symb(i)
c
c Subroutine symnum converts an element symbol to an atomic number.
c
          call symnum(stemp,Natmic)
          idtype(i)=Natmic
          nElect=nElect+nValen(Natmic)
          ndim=ndim+nOrbit(Natmic)
 20     continue
        stop "The system has two many atoms"
 21     natom=i-1
        first = .false.
      endif
c
c Convert the unit of the coordinate from bohr to angstrom.
c
      do 11 i = 1, 3*natom
         X(i)=X(i)*0.52917706d0
 11   continue
c
c Fill the diagonal elements of the Hamiltonian matrix.
c
      nrdim=0
c
      do 30 i=1,natom
        ip=idtype(i)
        ncdim=0
c
        do 31 j=1,natom
          jp=idtype(j)
c
          if (j.eq.i) then
            do 32 ii=1,nOrbit(ip)
              do 33 jj=1,nOrbit(jp)
                dhamlt(nrdim+ii,ncdim+jj,1)=0.0
                dhamlt(nrdim+ii,ncdim+jj,2)=0.0
                dhamlt(nrdim+ii,ncdim+jj,3)=0.0
                if (ii.eq.jj) then
                        St(nrdim+ii,ncdim+jj)=1.0
                  if (ii .eq. 1) then
                    hamilt(nrdim+ii,ncdim+jj)=vsip(ip,1)
                  else if (ii .ge. 2 .and. ii .le. 4) then
                    hamilt(nrdim+ii,ncdim+jj)=vsip(ip,2)
                  else
                    hamilt(nrdim+ii,ncdim+jj)=vsip(ip,3)
                  endif
                else 
                      St(nrdim+ii,ncdim+jj)=0.0
                  hamilt(nrdim+ii,ncdim+jj)=0.0
                endif
 33           continue
 32         continue
            goto 39
          else
c 
c Subroutine lmndlm calculate directional cosines and its first derivatives.
c
            call lmndlm(i,j,X)
c
c Subroutine ehuckl calculate the submatrice.
c
            call exhukl(ip,jp,hamil,dhamil,S)
c
            do 34 ii=1,nOrbit(ip)
              do 35 jj=1,nOrbit(jp)
                    St(nrdim+ii,ncdim+jj)=    S(ii,jj)
                hamilt(nrdim+ii,ncdim+jj)=hamil(ii,jj)
                dhamlt(nrdim+ii,ncdim+jj,1)=dhamil(ii,jj,1)
                dhamlt(nrdim+ii,ncdim+jj,2)=dhamil(ii,jj,2)
                dhamlt(nrdim+ii,ncdim+jj,3)=dhamil(ii,jj,3)
 35           continue
 34         continue
          endif
 39       ncdim=ncdim+nOrbit(jp)
 31     continue
        nrdim=nrdim+nOrbit(ip)
 30   continue
c
c Set up the uppermatrices to do eigenvalue calculations.
c
      do 40 i =1,ndim 
        do 41 j=i,ndim
          ap(i+(j-1)*j/2)=hamilt(i,j)
          bp(i+(j-1)*j/2)=    St(i,j)
 41     continue 
 40   continue
c
c Call library routine: lapack.  In t1 the path is:
c /usr/local/lib -llapack -lblas
c
      call dspgv(1,'v','u',ndim,ap,bp,eva,evec,max9,work,info)
c
c Call a subroutine to arrange all the valence electrons to the necessary
c molecular orbitals so that the potential energy of the system is minimized.
c
      call filloc(eva,V,occ,nElect,ndim)
c
c Calculate HOMO orbitals
c
      if (nchoice.eq.0) then
      do i = 1, ndim
        if (abs(occ(i)-2.0).lt.0.00001) then
          if (abs(occ(i+1)-1.0).lt.0.000001) then
            write(*,*) i,eva(i)," This is highest doubly occupied."
            write(*,*) i+1,eva(i+1)," This is lowest singly occupied."
          endif 
        else if (abs(occ(i)-1.0).lt.0.00001) then
          if (abs(occ(i+1)).lt.0.000001) then
            write(*,*) i,eva(i)," This is HOMO."
            write(*,*) i+1,eva(i+1)," This is LUMO. ",ndim
          endif
        endif
      enddo
      endif
c
c Calculate Mulliken net AO Population for each atomic orbital
c
      do j = 1, ndim
        q(j)=0.0
        do i = 1, ndim
          q(j) = q(j) + occ(i)*evec(j,i)*evec(j,i)
        enddo
      enddo
c
c Calculate a Mulliken overlap population.
c
      do j = 1, ndim
        do k = 1, ndim
          if (j .eq. k) then
            p(j, k)=q(j)
          else
            do i = 1, ndim
              p(j,k)=p(j,k)+2.d0*occ(i)*evec(j,i)*evec(k,i)*St(j,k)
            enddo
          endif
        enddo
      enddo
c
c Calculate the net charge.
c
      nrdim=0.0
      do i = 1, natom
        qcharg(i)=0.0
        ip=idtype(i)
        ncdim=0.0
        do j = 1, natom
          jp=idtype(j)
          do ii = 1, nOrbit(ip)
            if (i .eq. j) then
              qcharg(i)=qcharg(i)+q(nrdim+ii)
            else
              do jj = 1, nOrbit(jp)
                qcharg(i)=qcharg(i)+p(nrdim+ii,ncdim+jj)/2.0
              enddo
            endif
          enddo      
          ncdim=ncdim+nOrbit(jp)
        enddo
        nrdim=nrdim+nOrbit(ip)
      enddo
c
c Print Out eva, occ, q
c
      if (nchoice.eq.0) then
      print *, "--------------------------------------------------"
      do i=0,ndim/8
        if (i.lt.ndim/8) then
          print *, "Row ",i+1," for eigenvalues, occupation numbers, q"
          print 99, (eva(i*8+j),j=1,8)
          print 99, (occ(i*8+j),j=1,8)
          print 99, (q(i*8+j),j=1,8)
          print *, "--------------------------------------------------"
          do k=1,ndim
            print 99, (    St(i*8+j,k),j=1,8)
          enddo
          print *, "--------------------------------------------------"
          print *, "--------------------------------------------------"
          do k=1,ndim
            print 99, (hamilt(i*8+j,k),j=1,8)
          enddo
          print *, "--------------------------------------------------"
        else if (mod(ndim,8).ne.0) then
          print *, "Row ",i+1," for eigenvalues, occupation numbers, q"
          print 99, (eva(j),j=8*i+1,ndim)
          print 99, (occ(j),j=8*i+1,ndim)
          print 99, (q(j),j=8*i+1,ndim)
          print *, "--------------------------------------------------"
          do k=1,ndim
            print 99, (    St(j,k),j=8*i+1,ndim)
          enddo
          print *, "--------------------------------------------------"
          print *, "--------------------------------------------------"
          do k=1,ndim
            print 99, (hamilt(j,k),j=8*i+1,ndim)
          enddo
          print *, "--------------------------------------------------"
        endif
      enddo   
      print *, "--------------------------------------------------"
      endif
 99   format(8(f9.4,1x))
c
c Calculate the first derivative analytically.
c
      nrdim=0
c
      do 51 i=1,natom
        ip=idtype(i)
        ncdim=0
c
        dx(3*(i-1)+1)=0.0
        dx(3*(i-1)+2)=0.0
        dx(3*(i-1)+3)=0.0

        do 52 j=1,natom
          jp=idtype(j)
c
          if (j.eq.i) goto 59
              
          forcex=0.d0
          forcey=0.d0
          forcez=0.d0

          do 53 k=1,ndim
            do 54 ii=1,nOrbit(ip)
              do 55 jj=1,nOrbit(jp)
                  tmp=evec(nrdim+ii,k)*evec(ncdim+jj,k)*occ(k)
                  forcex=forcex+tmp*dhamlt(nrdim+ii,ncdim+jj,1)
                  forcey=forcey+tmp*dhamlt(nrdim+ii,ncdim+jj,2)
                  forcez=forcez+tmp*dhamlt(nrdim+ii,ncdim+jj,3)
 55             continue
 54           continue
 53         continue
            dx(3*(i-1)+1)=dx(3*(i-1)+1)+2.d0*forcex
            dx(3*(i-1)+2)=dx(3*(i-1)+2)+2.d0*forcey
            dx(3*(i-1)+3)=dx(3*(i-1)+3)+2.d0*forcez
 59         ncdim=ncdim+nOrbit(jp)
 52       continue
          nrdim=nrdim+nOrbit(ip)
c
 51     continue
c
c The binding energy is given in unit of eV, which is converted to atomic unit.
c Convert Cartesian coordinates from Angstrom to bohr, 
c first derivative from eV/(Angstrom) to Hartree/bohr.
c
      V=V/27.21161
c
      do 60 i = 1, 3*natom
        X(i)=X(i)/0.52917706d0
        DX(i)=DX(i)*0.019446738d0
 60   continue

c
c print out some information about the configuration.
c
      if (nchoice.eq.0) then
      print *, "--------------------------------------------------"
      print *," Element,  Electrons, Charges,     dx          dy
     *  dz(a.u./bohr)"
      do i=1,natom
        ip=idtype(i)
        print 101,symb(i),qcharg(i),nValen(ip)-qcharg(i),
     *         (dx(3*(i-1)+j),j=1,3)
      enddo
      print *, "--------------------------------------------------"
 101  format(7x,A2,5(f10.5,2x))
c
      print *,"INFO = ",info
      endif
      return
      end 
