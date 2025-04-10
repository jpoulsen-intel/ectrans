! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_MOD
CONTAINS
SUBROUTINE FSC(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PSCALAR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PEWDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PUVDERS(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = REAL(F%RACTHE(IGLG),JPRB)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2 = REAL(F%RACTHE(IGLG),JPRB)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = REAL(F%RACTHE2(IGLG),JPRB)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(PSCALAR(JF,IR)**2 + PSCALAR(JF,II)**2)
        ZPHASE = ATAN2(PSCALAR(JF,II),PSCALAR(JF,IR)) + REAL(JM,JPRB)*ZSHIFT
        
        PSCALAR(JF,IR) =  ZAMP*COS(ZPHASE)
        PSCALAR(JF,II) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS(JF,IR)**2 + PNSDERS(JF,II)**2)
          ZPHASE = ATAN2(PNSDERS(JF,II),PNSDERS(JF,IR)) + REAL(JM,JPRB)*ZSHIFT
          PNSDERS(JF,IR) =  ZAMP*COS(ZPHASE)
          PNSDERS(JF,II) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(PUV(JF,IR)**2 + PUV(JF,II)**2)
        ZPHASE = ATAN2(PUV(JF,II),PUV(JF,IR)) + REAL(JM,JPRB)*ZSHIFT
        PUV(JF,IR) =  ZAMP*COS(ZPHASE)
        PUV(JF,II) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,2*KF_UV
      PUV(JF,JLON) = PUV(JF,JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,KF_SCALARS
      PNSDERS(JF,JLON) = PNSDERS(JF,JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,2*KF_UV
      PUVDERS(JF,IR) = -PUV(JF,II)*ZMUL
      PUVDERS(JF,II) =  PUV(JF,IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,KF_SCALARS
      PEWDERS(JF,IR) = -PSCALAR(JF,II)*ZMUL
      PEWDERS(JF,II) =  PSCALAR(JF,IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC

#ifdef CPU_FFTW_INPLACE

SUBROUTINE FSC_TR(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV_TR,PSCALAR_TR,PNSDERS_TR,PEWDERS_TR,PUVDERS_TR)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV_TR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PSCALAR_TR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS_TR(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PEWDERS_TR(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PUVDERS_TR(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(PSCALAR_TR(IR,JF)**2 + PSCALAR_TR(II,JF)**2)
        ZPHASE = ATAN2(PSCALAR_TR(II,JF),PSCALAR_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        
        PSCALAR_TR(IR,JF) =  ZAMP*COS(ZPHASE)
        PSCALAR_TR(II,JF) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS_TR(IR,JF)**2 + PNSDERS_TR(II,JF)**2)
          ZPHASE = ATAN2(PNSDERS_TR(II,JF),PNSDERS_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
          PNSDERS_TR(IR,JF) =  ZAMP*COS(ZPHASE)
          PNSDERS_TR(II,JF) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(PUV_TR(IR,JF)**2 + PUV_TR(II,JF)**2)
        ZPHASE = ATAN2(PUV_TR(II,JF),PUV_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        PUV_TR(IR,JF) =  ZAMP*COS(ZPHASE)
        PUV_TR(II,JF) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
   DO JF=1,2*KF_UV
      DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
         PUV_TR(JLON,JF) = PUV_TR(JLON,JF)*ZACHTE2
      END DO
   END DO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
   DO JF=1,KF_SCALARS
      DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
         PNSDERS_TR(JLON,JF) = PNSDERS_TR(JLON,JF)*ZACHTE2
      END DO
   END DO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
   DO JF=1,2*KF_UV
      DO JM=0,IMEN
         IR = ISTAGTF+2*JM+1
         II = IR+1
         ZMUL = ZACHTE*JM
         PUVDERS_TR(IR,JF) = -PUV_TR(II,JF)*ZMUL
         PUVDERS_TR(II,JF) =  PUV_TR(IR,JF)*ZMUL
      END DO
   END DO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
   DO JF=1,KF_SCALARS
      DO JM=0,IMEN
         IR = ISTAGTF+2*JM+1
         II = IR+1
         ZMUL = ZACHTE*JM
         PEWDERS_TR(IR,JF) = -PSCALAR_TR(II,JF)*ZMUL
         PEWDERS_TR(II,JF) =  PSCALAR_TR(IR,JF)*ZMUL
      END DO
   END DO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_TR

SUBROUTINE FSC_TR_SCALAR(KGL,KF_SCALARS,KF_SCDERS,PSCALAR_TR,PNSDERS_TR,PEWDERS_TR,JF)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_SCALARS,KF_SCDERS,JF
! REAL(KIND=JPRB) , INTENT(INOUT) :: PUV_TR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PSCALAR_TR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PNSDERS_TR(:,:)
REAL(KIND=JPRB) , INTENT(  OUT), CONTIGUOUS :: PEWDERS_TR(:,:)
! REAL(KIND=JPRB) , INTENT(  OUT) :: PUVDERS_TR(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    ! DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(PSCALAR_TR(IR,JF)**2 + PSCALAR_TR(II,JF)**2)
        ZPHASE = ATAN2(PSCALAR_TR(II,JF),PSCALAR_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        
        PSCALAR_TR(IR,JF) =  ZAMP*COS(ZPHASE)
        PSCALAR_TR(II,JF) = ZAMP*SIN(ZPHASE)
      ENDDO
    ! ENDDO
    IF(KF_SCDERS > 0)THEN
      ! DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS_TR(IR,JF)**2 + PNSDERS_TR(II,JF)**2)
          ZPHASE = ATAN2(PNSDERS_TR(II,JF),PNSDERS_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
          PNSDERS_TR(IR,JF) =  ZAMP*COS(ZPHASE)
          PNSDERS_TR(II,JF) = ZAMP*SIN(ZPHASE)
        ENDDO
      ! ENDDO
    ENDIF
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
   ! DO JF=1,KF_SCALARS
      DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
         PNSDERS_TR(JLON,JF) = PNSDERS_TR(JLON,JF)*ZACHTE2
      END DO
   ! END DO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
   ! DO JF=1,KF_SCALARS
      DO JM=0,IMEN
         IR = ISTAGTF+2*JM+1
         II = IR+1
         ZMUL = ZACHTE*JM
         PEWDERS_TR(IR,JF) = -PSCALAR_TR(II,JF)*ZMUL
         PEWDERS_TR(II,JF) =  PSCALAR_TR(IR,JF)*ZMUL
      END DO
   ! END DO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_TR_SCALAR

SUBROUTINE FSC_SCALAR_1D(KGL,KF_SCALARS,KF_SCDERS,PSCALAR,PNSDERS,PEWDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PSCALAR(:)
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PNSDERS(:)
REAL(KIND=JPRB) , INTENT(  OUT), CONTIGUOUS :: PEWDERS(:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF

INTEGER(KIND=JPIM) :: JLON,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF ( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JM=0,IMEN
       IR = ISTAGTF+2*JM+1
       II = IR+1
       
       ! calculate amplitude and add phase shift then reconstruct A,B
       ZAMP = SQRT(PSCALAR(IR)**2 + PSCALAR(II)**2)
       ZPHASE = ATAN2(PSCALAR(II),PSCALAR(IR)) + REAL(JM,JPRB)*ZSHIFT
       
       PSCALAR(IR) =  ZAMP*COS(ZPHASE)
       PSCALAR(II) = ZAMP*SIN(ZPHASE)
    END DO
    IF (KF_SCDERS > 0)THEN
       DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS(IR)**2 + PNSDERS(II)**2)
          ZPHASE = ATAN2(PNSDERS(II),PNSDERS(IR)) + REAL(JM,JPRB)*ZSHIFT
          PNSDERS(IR) =  ZAMP*COS(ZPHASE)
          PNSDERS(II) = ZAMP*SIN(ZPHASE)
        END DO
    END IF
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

!*      1.2      N-S DERIVATIVES

IF (KF_SCDERS > 0)THEN
   DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
      PNSDERS(JLON) = PNSDERS(JLON)*ZACHTE2
   END DO
END IF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
   DO JM=0,IMEN
      IR = ISTAGTF+2*JM+1
      II = IR+1
      ZMUL = ZACHTE*JM
      PEWDERS(IR) = -PSCALAR(II)*ZMUL
      PEWDERS(II) =  PSCALAR(IR)*ZMUL
   END DO
ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE FSC_SCALAR_1D

SUBROUTINE FSC_TR_UV(KGL,KF_UV,PUV_TR,PUVDERS_TR, JF)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,JF
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PUV_TR(:,:)
REAL(KIND=JPRB) , INTENT(  OUT), CONTIGUOUS :: PUVDERS_TR(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    ! DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(PUV_TR(IR,JF)**2 + PUV_TR(II,JF)**2)
        ZPHASE = ATAN2(PUV_TR(II,JF),PUV_TR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        PUV_TR(IR,JF) =  ZAMP*COS(ZPHASE)
        PUV_TR(II,JF) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ! ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
   ! DO JF=1,2*KF_UV
      DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
         PUV_TR(JLON,JF) = PUV_TR(JLON,JF)*ZACHTE2
      END DO
   ! END DO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
   ! DO JF=1,2*KF_UV
      DO JM=0,IMEN
         IR = ISTAGTF+2*JM+1
         II = IR+1
         ZMUL = ZACHTE*JM
         PUVDERS_TR(IR,JF) = -PUV_TR(II,JF)*ZMUL
         PUVDERS_TR(II,JF) =  PUV_TR(IR,JF)*ZMUL
      END DO
   ! END DO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_TR_UV


SUBROUTINE FSC_UV_1D(KGL,KF_UV,PUV,PUVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV
REAL(KIND=JPRB) , INTENT(INOUT), CONTIGUOUS :: PUV(:)
REAL(KIND=JPRB) , INTENT(  OUT), CONTIGUOUS :: PUVDERS(:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JM=0,IMEN
       IR = ISTAGTF+2*JM+1
       II = IR+1
       ! calculate amplitude and phase shift and reconstruct A,B
       ZAMP = SQRT(PUV(IR)**2 + PUV(II)**2)
       ZPHASE = ATAN2(PUV(II),PUV(IR)) + REAL(JM,JPRB)*ZSHIFT
       PUV(IR) =  ZAMP*COS(ZPHASE)
       PUV(II) =  ZAMP*SIN(ZPHASE)
    END DO
  END IF
END IF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
   DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
      PUV(JLON) = PUV(JLON)*ZACHTE2
   END DO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (LUVDER)THEN
   DO JM=0,IMEN
      IR = ISTAGTF+2*JM+1
      II = IR+1
      ZMUL = ZACHTE*JM
      PUVDERS(IR) = -PUV(II)*ZMUL
      PUVDERS(II) =  PUV(IR)*ZMUL
   END DO
END IF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_UV_1D
#endif
END MODULE FSC_MOD
