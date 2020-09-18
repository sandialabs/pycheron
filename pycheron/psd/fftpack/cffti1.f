!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
! certain rights in this software.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTICE:
! For five (5) years from 10/21/2019 the United States Government is granted for
! itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
! license in this data to reproduce, prepare derivative works, and perform publicly and
! display publicly, by or on behalf of the Government. There is provision for the
! possible extension of the term of this license. Subsequent to that period or any
! extension granted, the United States Government is granted for itself and others
! acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
! data to reproduce, prepare derivative works, distribute copies to the public,
! perform publicly and display publicly, and to permit others to do so. The specific
! term of the license can be identified by inquiry made to National Technology and
! Engineering Solutions of Sandia, LLC or DOE. NEITHER THE UNITED STATES GOVERNMENT,
! NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING
! SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
! IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
! USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
! THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. Any licensee of this software
! has the obligation and responsibility to abide by the applicable export control laws,
! regulations, and general prohibitions relating to the export of technical data.
! Failure to obtain an export control license or other authority from the Government
! may result in criminal liability under U.S. laws.
! (End of Notice)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE CFFTI1 (N,WA,IFAC)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J.le.4) GO TO 102
      GO TO 103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRYH(4)+2*(J-4)
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.eq.0) GO TO 105
      GO TO 101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD+L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
