c This subroutine symnum converts an element symbol to an atomic number.
c
      subroutine symnum(stemp,Natmic)
c
      implicit none

      integer Natmic,i        ! i --- temporary loop number
c
c sym_st --- an array holding 103 elements
c
      character*2 stemp,sym_st(103)
c
      DATA (sym_st(i),i=1,103)/'H ','HE',
     2 'LI','BE','B ','C ','N ','O ','F ','NE',
     3 'NA','MG','AL','SI','P ','S ','CL','AR',
     4 'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG',
     5 'CD','IN','SN','SB','TE','I ','XE',
     6 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     6 'HO','ER','TM','YB','LU','HF','TA','W ','RE','OS','IR','PT',
     6 'AU','HG','TL','PB','BI','PO','AT','RN',
     7 'FR','RA','AC','TH','PA','U ','NP','PU','AM','CM','BK','CF','ES',
     8 'FM','MD','NO','LR'/
c
      do 1 i = 1, 103
        if (stemp .eq. sym_st(i)) then
          Natmic=i
          return
        endif
 1    continue
c
      stop "Unknown element symbol, make sure you put symbol correct!"
c
      return
      end
