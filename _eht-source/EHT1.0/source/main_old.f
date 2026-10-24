       program main

      implicit none

      integer i,j,nchoice,n3tm
      character*2 symb(60)
      double precision X(180),V,DX(180)

      common /choice/ nchoice

      n3tm = 180

      do i = 1, 60
        read (98,*,end=21) symb(i), (X(3*(i-1)+j),j=1,3)
      enddo

 21   continue

      print *, "---nothing special 0---"
      print *, "---nic 1---"
      print *, "---nih 2---"
      print *, "---nich3 3---"
      print *, "---ni4ch3 4---"
      print *, "---ch 5---"
      print *, "---cc 6---"
      print *, "---hh 7---"
      print *, "---nini 8---"

      read *, nchoice
      if (nchoice.eq.0) then
        call surf(V,X,DX,N3TM)
        print *," V= ",V
      endif
c --------------------------------------------------------
      do i = 1, 100
        if (nchoice.eq.1) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+6.19294484964319)*27.21161
        else if (nchoice.eq.2) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+4.564228283442219)*27.21161
        else if (nchoice.eq.3) then
          X(4)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+8.105693363214016)*27.21161
        else if (nchoice.eq.4) then
          X(15)=i*0.05/0.52917706d0
          X(18)=X(15)+0.68912
          X(21)=X(15)+0.68912
          X(24)=X(15)+0.68912
          call surf(V,X,DX,N3TM) 
          print *, i*0.05, (v+24.44264273974898d0)*27.21161
        else if (nchoice.eq.41) then
          X(15)=i*0.05/0.52917706d0
          X(18)=X(15)+0.68912
          X(21)=X(15)+0.68912
          X(24)=X(15)+0.68912
          call surf(V,X,DX,N3TM) 
          print *, i*0.05, (v+24.35521057809603d0)*27.21161
        else if (nchoice.eq.5) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+2.991370227634434)*27.21161
        else if (nchoice.eq.6) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+4.821471423410773)*27.21161
        else if (nchoice.eq.7) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+0.926075303886833)*27.21161
        else if (nchoice.eq.8) then
          X(6)=i*0.05/0.52917706
          call surf(V,X,DX,N3TM)
          print *,i*0.05,(v+9.914885594788359)*27.21161
        endif
      enddo 

      end
