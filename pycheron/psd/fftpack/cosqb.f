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
      
      SUBROUTINE COSQB (N,X,WSAVE)
      DIMENSION       X(*)       ,WSAVE(*)
      DATA TSQRT2 /2.82842712474619/
      IF (N.lt.2) GO TO 101
      IF (N.eq.2) GO TO 102
      GO TO 103
  101 X(1) = 4.*X(1)
      RETURN
  102 X1 = 4.*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE COSQB1 (N,X,W,XH)
      DIMENSION       X(1)       ,W(1)       ,XH(1)
      NS2 = (N+1)/2
      NP2 = N+2
      DO 101 I=3,N,2
         XIM1 = X(I-1)+X(I)
         X(I) = X(I)-X(I-1)
         X(I-1) = XIM1
  101 CONTINUE
      X(1) = X(1)+X(1)
      MODN = MOD(N,2)
      IF (MODN .EQ. 0) X(N) = X(N)+X(N)
      CALL RFFTB (N,X,XH)
      DO 102 K=2,NS2
         KC = NP2-K
         XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)
  102 CONTINUE
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO 103 K=2,NS2
         KC = NP2-K
         X(K) = XH(K)+XH(KC)
         X(KC) = XH(K)-XH(KC)
  103 CONTINUE
      X(1) = X(1)+X(1)
      RETURN
      END
