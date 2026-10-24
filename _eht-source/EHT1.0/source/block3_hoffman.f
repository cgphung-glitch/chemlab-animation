      BLOCK DATA set3

      implicit none
c
c Temporary variable.
c
      integer i
c
c nOrbit --- The number of valence orbitals for an atom.
c nValen --- The number of valence electrons for an atom.
c nsd --- The principal quantum number of s orbital
c npd --- The principal quantum number of p orbital
c ndd --- The principal quantum number of d orbital
c vsip --- valence state ionilization potentials
c
      integer nOrbit(103),nValen(103),nsd(103),npd(103),ndd(103)
      double precision vsip(103,3),cdd1(103),cdd2(103)
      double precision expss(103),expps(103),expdd1(103),expdd2(103)
c
c Common block "orbit" will share with subroutine "BLOCK DATA"
c Common block "elctrn" will share with subroutine "BLOCK DATA"
c Common block "nspdd" will share with subroutine "ehuckl,BLOCK DATA"
c Common block "coulom" will share with subroutine "ehuckl,BLOCK DATA"
c Common block "dcoeff" will share with subroutine "ehuckl,BLOCK DATA"
c Common block "dtzeta" will share with subroutine "ehuckl,BLOCK DATA"
c
      common /orbit/ nOrbit
      common /elctrn/ nValen
      common /nspdd/ nsd,npd,ndd
      common /coulom/ vsip
      common /dcoeff/ cdd1,cdd2
      common /dtzeta/ expss,expps,expdd1,expdd2

c
c The data here are taken from the code of QCPE 571.
c
c ********************* valence electrons ************************
c
c                                 H HE LI BE  B  C  N  O  F NE
      data (nValen(i),i=1,103)  / 1, 2, 1, 2, 3, 4, 5, 6, 7, 8,
c     NA MG AL SI  P  S CL AR  K CA SC TI  V CR MN FE CO NI CU ZN
     1 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 2,
c     GA GE AS SE BR KR RB SR  Y ZR NB MO TC RU RH PD AG CD IN SN
     2 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 3, 4,
c     SB TE  I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB
     3 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
c     LU HF TA  W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH
     417, 4, 5, 6, 7, 8, 9,10,11,12, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4,
c     PA  U NP PU AM CM BK CF ES FM MD NO LR
     5 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17/
c
c *******************  valence orbitals  *****************************
c
c         s:
c                            H HE LI BE  B  C  N  O  F NE
      data (nsd(i),i=1,103)/ 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
c     NA MG AL SI  P  S CL AR  K CA SC TI  V CR MN FE CO NI CU ZN
     1 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
c     GA GE AS SE BR KR RB SR  Y ZR NB MO TC RU RH PD AG CD IN SN
     2 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
c     SB TE  I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB
     3 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
c     LU HF TA  W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH
     4 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7,
c     PA  U NP PU AM CM BK CF ES FM MD NO LR
     5 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7/
c    
c         P:
c                            H HE LI BE  B  C  N  O  F NE
      data (npd(i),i=1,103)/ 0, 0, 2, 2, 2, 2, 2, 2, 2, 2,
c     NA MG AL SI  P  S CL AR  K CA SC TI  V CR MN FE CO NI CU ZN
     1 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
c     GA GE AS SE BR KR RB SR  Y ZR NB MO TC RU RH PD AG CD IN SN
     2 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
c     SB TE  I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB
     3 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
c     LU HF TA  W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH
     4 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7,
c     PA  U NP PU AM CM BK CF ES FM MD NO LR
     5 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7/
c    
c          d:
c                            H HE LI BE  B  C  N  O  F NE
      data (ndd(i),i=1,103)/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c     NA MG AL SI  P  S CL AR  K CA SC TI  V CR MN FE CO NI CU ZN
     1 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0,
c     GA GE AS SE BR KR RB SR  Y ZR NB MO TC RU RH PD AG CD IN SN
     2 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0,
c     SB TE  I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB
     3 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
c     LU HF TA  W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH
     4 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6,
c     PA  U NP PU AM CM BK CF ES FM MD NO LR
     5 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/
c    
c Number of orbitals  *********************  Valence orbitals
c
c                                 H HE LI BE  B  C  N  O  F NE
      data (nOrbit(i),i=1,103)  / 1, 1, 4, 4, 4, 4, 4, 4, 4, 4,
c     NA MG AL SI  P  S CL AR  K CA SC TI  V CR MN FE CO NI CU ZN
     1 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
c     GA GE AS SE BR KR RB SR  Y ZR NB MO TC RU RH PD AG CD IN SN
     2 4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4,
c     SB TE  I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB
     3 4, 4, 4, 4, 4, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
c     LU HF TA  W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH
     4 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 9, 9,
c     PA  U NP PU AM CM BK CF ES FM MD NO LR
     5 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9/
c
c

c IP's of s,p,d and f *********************
c
c    s:
c
c                               H       HE      LI      BE      B
      data (vsip(i,1),i=1,69)/-13.600,-23.400, -5.400,-10.000,-15.200,
c      C       N       O       F       NE      NA      MG      AL
     1-21.400,-26.000,-32.300,-40.000,-43.200, -5.100, -9.000,-12.300,
c      SI      P       S       CL      AR      K       CA      SC
     2-17.300,-18.600,-20.000,-26.300, -0.000, -4.340, -7.000, -8.870,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 -8.970, -8.810, -8.660, -9.750, -9.100, -9.210, -7.800,-11.400,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4-12.410,-14.580,-16.000,-16.220,-20.500,-22.070, -0.000, -4.180,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 -6.620, -0.000, -8.000,-10.100, -8.340,-10.070,-10.400, -8.090,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 -7.320, -0.000, -0.000,-12.600,-16.160,-18.800,-20.800,-18.000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 -0.000, -3.880, -0.000, -7.670, -0.000, -0.000, -0.000, -0.000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 -4.860, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c                                YB      LU      HF      TA
      data (vsip(i,1),i=70,103)/ -5.350, -6.050, -0.000,-10.100,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 -8.260, -9.360, -8.170,-11.360, -9.077,-10.920,-13.680,-11.600,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2-15.700,-15.190, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 -5.390, -0.000, -5.500, -0.000, -0.000, -0.000, -0.000, -0.000,
c      CF      ES      FM      MD      NO      LR
     4 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c
c
c    p:
c
c                               H       HE      LI      BE      B
      data (vsip(i,2),i=1,69)/ -0.000, -0.000, -3.500, -6.000, -8.500,
c      C       N       O       F       NE      NA      MG      AL
     1-11.400,-13.400,-14.800,-18.100,-20.000, -3.000, -4.500, -6.500,
c      SI      P       S       CL      AR      K       CA      SC
     2 -9.200,-14.000,-13.300,-14.200, -0.000, -2.730, -4.000, -2.750,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 -5.440, -5.520, -5.240, -5.890, -5.320, -5.290, -3.700, -6.060,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 -6.530, -6.750, -9.000,-12.160,-14.400,-13.100, -0.000, -2.600,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 -3.920, -0.000, -5.400, -6.860, -5.240, -5.400, -6.870, -4.570,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 -3.750, -0.000, -0.000, -6.190, -8.320,-11.700,-13.200,-12.700,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 -0.000, -2.490, -0.000, -5.010, -0.000, -0.000, -0.000, -0.000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 -4.860, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c                                YB      LU      HF      TA
      data (vsip(i,2),i=70,103)/ -5.350, -6.050, -0.000, -6.860,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 -5.170, -5.960, -4.810, -4.500, -5.475, -5.550, -8.470, -5.800,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 -8.000, -7.790, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 -5.390, -0.000, -5.500, -0.000, -0.000, -0.000, -0.000, -0.000,
c      CF      ES      FM      MD      NO      LR
     4 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c
c
c    d:
c
c                               H       HE      LI      BE      B
      data (vsip(i,3),i=1,69)/ -0.000, -0.000, -0.000, -0.000, -0.000,
c      C       N       O       F       NE      NA      MG      AL
     1 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      SI      P       S       CL      AR      K       CA      SC
     2 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -8.510,
c      TI      V       CR      MN      FE      CO      NI      CU
     3-10.810,-11.000,-11.220,-11.670,-12.600,-13.180,-9.9000,-14.000,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 -0.000, -0.000,-10.200,-12.100,-10.500,-12.820,-14.900,-12.500,
c      PD      AG      CD      IN      SN      SB      TE      I
     6-12.020, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 -0.000, -0.000, -0.000, -8.210, -0.000, -0.000, -0.000, -0.000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 -6.060, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c                                YB      LU      HF      TA
      data (vsip(i,3),i=70,103)/ -5.210, -5.120, -0.000,-12.100,
c      W       RE      OS      IR      PT      AU      HG      TL
     1-10.370,-12.660,-11.840,-12.170,-12.590,-15.076,-17.500, -0.000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000, -0.000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3-10.110, -0.000, -9.190, -0.000, -0.000, -0.000, -0.000, -0.000,
c      CF      ES      FM      MD      NO      LR
     4 -0.000, -0.000, -0.000, -0.000, -0.000, -0.000/
c
c      d double zeta coefficient #1
c
c                            H       HE      LI      BE      B
      data (cdd1(i),i=1,69)/0.00000,0.00000,0.00000,0.00000,0.00000,
c      C       N       O       F       NE      NA      MG      AL
     10.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      SI      P       S       CL      AR      K       CA      SC
     20.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.42280,
c      TI      V       CR      MN      FE      CO      NI      CU
     30.42060,0.47550,0.50600,0.51390,0.55050,0.56800,0.56830,0.59330,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     40.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     50.00000,0.00000,0.62100,0.64010,0.60970,0.57150,0.53400,0.58070,
c      PD      AG      CD      IN      SN      SB      TE      I
     60.55350,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     70.00000,0.00000,0.00000,0.77650,0.00000,0.00000,0.00000,0.00000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     80.71840,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000/
c                              YB      LU      HF      TA
      data (cdd1(i),i=70,103)/0.70630,0.70440,0.00000,0.68150,
c      W       RE      OS      IR      PT      AU      HG      TL
     10.69400,0.66620,0.63720,0.66980,0.63340,0.68510,0.64380,0.00000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     20.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      TH      PA      U       NP      PU      AM      CM      BK
     30.76120,0.00000,0.76080,0.00000,0.00000,0.00000,0.00000,0.00000,
c      CF      ES      FM      MD      NO      LR
     40.00000,0.00000,0.00000,0.00000,0.00000,0.00000/
c
c
c
c      d double zeta coefficient #2
c
c                            H       HE      LI      BE      B
      data (cdd2(i),i=1,69)/0.00000,0.00000,0.00000,0.00000,0.00000,
c      C       N       O       F       NE      NA      MG      AL
     10.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      SI      P       S       CL      AR      K       CA      SC
     20.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.72760,
c      TI      V       CR      MN      FE      CO      NI      CU
     30.78390,0.70520,0.67500,0.69290,0.62600,0.60600,0.62920,0.57440,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     40.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     50.00000,0.00000,0.57690,0.55160,0.60970,0.60120,0.63650,0.56850,
c      PD      AG      CD      IN      SN      SB      TE      I
     60.67010,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     70.00000,0.00000,0.00000,0.45860,0.00000,0.00000,0.00000,0.00000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     80.44470,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000/
c                              YB      LU      HF      TA
      data (cdd2(i),i=70,103)/0.48340,0.48800,0.00000,0.68150,
c      W       RE      OS      IR      PT      AU      HG      TL
     10.56310,0.59100,0.55980,0.58600,0.55130,0.56960,0.52150,0.00000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     20.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
c      TH      PA      U       NP      PU      AM      CM      BK
     30.40710,0.00000,0.41260,0.00000,0.00000,0.00000,0.00000,0.00000,
c      CF      ES      FM      MD      NO      LR
     40.00000,0.00000,0.00000,0.00000,0.00000,0.00000/
c
c
c      s single zeta exponent
c
c                              H       HE      LI      BE      B
      data (expss(i),i=1,69)/ 1.3000, 1.6880, 0.6500, 0.9750, 1.3000,
c      C       N       O       F       NE      NA      MG      AL
     1 1.6250, 1.9500, 2.2750, 2.4250, 2.8790, 0.7330, 1.1000, 1.1670,
c      SI      P       S       CL      AR      K       CA      SC
     2 1.3830, 1.7500, 2.1220, 2.1830, 0.0000, 0.8740, 1.2000, 1.3000,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 1.0750, 1.3000, 1.7000, 0.9700, 1.9000, 2.0000, 2.1000, 2.2000,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 2.0100, 1.7700, 2.1600, 2.2300, 2.4400, 2.5880, 0.0000, 0.9970,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 1.2140, 0.0000, 1.8170, 1.8900, 1.9600, 2.0180, 2.0800, 2.1350,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 2.1900, 0.0000, 0.0000, 1.9030, 2.1200, 2.3230, 2.5100, 2.6790,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 0.0000, 1.0600, 0.0000, 2.1400, 0.0000, 0.0000, 0.0000, 0.0000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 1.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c                                YB      LU      HF      TA
      data (expss(i),i=70,103)/ 1.5400, 1.6660, 0.0000, 2.2800,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 2.3410, 2.3980, 2.4520, 2.5000, 2.5540, 2.6020, 2.6490, 2.3000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 2.3500, 2.5600, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 1.8340, 0.0000, 1.9140, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      CF      ES      FM      MD      NO      LR
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c
c                              H       HE      LI      BE      B
      data (expps(i),i=1,69)/ 0.0000, 0.0000, 0.6500, 0.9750, 1.3000,
c      C       N       O       F       NE      NA      MG      AL
     1 1.6250, 1.9500, 2.2750, 2.4250, 2.8790, 0.7330, 1.1000, 1.1670,
c      SI      P       S       CL      AR      K       CA      SC
     2 1.3830, 1.3000, 1.8270, 1.7330, 0.0000, 0.8740, 1.2000, 1.3000,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 0.6750, 1.3000, 1.7000, 0.9700, 1.9000, 2.0000, 2.1000, 2.2000,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 1.7000, 1.5500, 1.8500, 1.8900, 2.0700, 2.1310, 0.0000, 0.9970,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 1.2140, 0.0000, 1.7760, 1.8500, 1.9000, 1.9840, 2.0400, 2.1000,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 2.1520, 0.0000, 0.0000, 1.6770, 1.8200, 1.9990, 2.1600, 2.3220,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 0.0000, 1.0600, 0.0000, 2.0800, 0.0000, 0.0000, 0.0000, 0.0000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 1.4000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c                                YB      LU      HF      TA
c                               H       HE      LI      BE      B
      data (expdd1(i),i=1,69)/ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      C       N       O       F       NE      NA      MG      AL
     1 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      SI      P       S       CL      AR      K       CA      SC
     2 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 4.3500,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 4.5500, 4.7500, 4.9500, 5.1500, 5.3500, 5.5500, 5.7500, 5.9500,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 0.0000, 0.0000, 3.8350, 4.0800, 4.5400, 4.9000, 5.3800, 4.2900,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 5.9830, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 0.0000, 0.0000, 0.0000, 3.7800, 0.0000, 0.0000, 0.0000, 0.0000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 2.7470, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c                                 YB      LU      HF      TA
      data (expps(i),i=70,103)/ 1.5400, 1.6660, 0.0000, 2.2410,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 2.3090, 2.3720, 2.4290, 2.2000, 2.5540, 2.5840, 2.6310, 1.6000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 2.0600, 2.0720, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 1.8340, 0.0000, 1.9140, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      CF      ES      FM      MD      NO      LR
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c
      data (expdd1(i),i=70,103)/ 2.8100, 2.8130, 0.0000, 4.7620,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 4.9820, 5.3430, 5.5710, 5.7960, 6.0130, 6.1630, 6.4360, 0.0000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 2.4610, 0.0000, 2.5810, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      CF      ES      FM      MD      NO      LR
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c
c
c      d double zeta exponent #2
c
c                               H       HE      LI      BE      B
      data (expdd2(i),i=1,69)/ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      C       N       O       F       NE      NA      MG      AL
     1 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      SI      P       S       CL      AR      K       CA      SC
     2 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.7000,
c      TI      V       CR      MN      FE      CO      NI      CU
     3 1.4000, 1.7000, 1.8000, 1.7000, 2.0000, 2.1000, 2.0000, 2.3000,
c      ZN      GA      GE      AS      SE      BR      KR      RB
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      SR      Y       ZR      NB      MO      TC      RU      RH
     5 0.0000, 0.0000, 1.5050, 1.6400, 1.9000, 2.0940, 2.3000, 1.9700,
c      PD      AG      CD      IN      SN      SB      TE      I
     6 2.6130, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      XE      CS      BA      LA      CE      PR      ND      PM
     7 0.0000, 0.0000, 0.0000, 1.3810, 0.0000, 0.0000, 0.0000, 0.0000,
c      SM      EU      GD      TB      DY      HO      ER      TM
     8 1.2670, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c                                 YB      LU      HF      TA
      data (expdd2(i),i=70,103)/ 1.2160, 1.2100, 0.0000, 1.9380,
c      W       RE      OS      IR      PT      AU      HG      TL
     1 2.0680, 2.2770, 2.4160, 2.5570, 2.6960, 2.7940, 3.0320, 0.0000,
c      PB      BI      PO      AT      RN      FR      RA      AC
     2 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      TH      PA      U       NP      PU      AM      CM      BK
     3 1.1650, 0.0000, 1.2070, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
c      CF      ES      FM      MD      NO      LR
     4 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
c

      end
