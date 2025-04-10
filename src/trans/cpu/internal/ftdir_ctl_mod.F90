! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_CTL_MOD
  USE, INTRINSIC :: IEEE_EXCEPTIONS
  LOGICAL :: IEEE_HALT_FLAGS(SIZE(IEEE_USUAL))
#ifdef IFS_RAPS
  LOGICAL, PARAMETER :: IFS_CONTEXT = .TRUE.
#else
  LOGICAL, PARAMETER :: IFS_CONTEXT = .FALSE.
#endif
CONTAINS
SUBROUTINE FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS, &
 & KVSETUV,KVSETSC,KPTRGP,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!      R. El Khatib 01-Jun-2022 contiguous pointer
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : FOUBUF_IN, REALLOCATE_ARRAY
USE TPM_DISTR       ,ONLY : D

USE TRGTOL_MOD      ,ONLY : TRGTOL
USE FOURIER_OUT_MOD ,ONLY : FOURIER_OUT,FOURIER_OUT_SETUP
USE FTDIR_MOD       ,ONLY : FTDIR
!

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KPTRGP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), CONTIGUOUS :: KVSETSC2(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN), CONTIGUOUS :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

! Local variables
REAL(KIND=JPRB) :: PIN(D%NLENGTF,KF_FS)

INTEGER(KIND=JPIM) :: IST,JGL,IBLEN,JK
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3

!     ------------------------------------------------------------------
! Field distribution in Spectral/Fourier space

IF (IFS_CONTEXT) THEN
  ! Workaround MKL bug: FFTW triggers spurious FP exceptions, especially when AVX512 is in use
  CALL IEEE_GET_HALTING_MODE(IEEE_USUAL,IEEE_HALT_FLAGS)
  IF (ANY(IEEE_HALT_FLAGS == .TRUE.)) CALL IEEE_SET_HALTING_MODE(IEEE_USUAL,.FALSE.)
ENDIF

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:) = -1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
ENDIF

! Transposition

CALL GSTATS(158,0)
CALL TRGTOL(PIN,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
CALL GSTATS(158,1)
CALL GSTATS(106,0)

! Fourier transform

IBLEN=D%NLENGT0B*2*KF_FS
CALL REALLOCATE_ARRAY(FOUBUF_IN, MAX(1,IBLEN))

CALL GSTATS(1640, 0)

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) PRIVATE(JGL,JK)
DO JGL = 1, D%NDGL_FS
  DO JK = 1, KF_FS
    ! Fourier transform
    CALL FTDIR(PIN, JGL, JK)

    ! Save Fourier data in FOUBUF_IN
    CALL FOURIER_OUT(PIN, KF_FS, JGL, JK)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL GSTATS(1640, 1)

CALL GSTATS(106,1)

IF (IFS_CONTEXT) THEN
  ! Workaround MKL bug: FFTW triggers spurious FP exceptions, especially when AVX512 is in use
  IF (ANY(IEEE_HALT_FLAGS == .TRUE.)) CALL IEEE_SET_HALTING_MODE(IEEE_USUAL,IEEE_HALT_FLAGS)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTL
END MODULE FTDIR_CTL_MOD



