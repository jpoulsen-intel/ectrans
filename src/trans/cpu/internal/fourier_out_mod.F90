! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!                           KGL - local index of latitude we are currently on
!                           JF - index of field we are currently on
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------


MODULE FOURIER_OUT_MOD
USE PARKIND1   ,ONLY : JPIM
PRIVATE
PUBLIC FOURIER_OUT, FOURIER_OUT_SETUP
INTERFACE FOURIER_OUT
  MODULE PROCEDURE FOURIER_OUT_2D_NEW
END INTERFACE

INTEGER(KIND=JPIM), ALLOCATABLE :: ISTAIDX(:,:)

CONTAINS

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
SUBROUTINE FOURIER_OUT_2D_NEW(PREEL, KFIELDS, KGL, JF)

USE PARKIND1,     ONLY : JPIM, JPRB
USE TPM_DISTR,    ONLY : D, MYSETW
USE TPM_TRANS,    ONLY : FOUBUF_IN
USE TPM_GEOMETRY, ONLY : G

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(IN), CONTIGUOUS :: PREEL(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN) :: KGL
INTEGER(KIND=JPIM), INTENT(IN) :: JF

INTEGER(KIND=JPIM) :: JM, IR, II, IGLG
INTEGER(KIND=JPIM) :: OFFSET1, OFFSET2

IGLG   = D%NPTRLS(MYSETW) + KGL - 1
OFFSET1 = D%NSTAGTF(KGL)+1
OFFSET2 = 2*JF-1

DO JM = 0, G%NMEN(IGLG)
  IR = 2 * JM + OFFSET1
  II = ISTAIDX(JM,KGL)*KFIELDS + OFFSET2
  FOUBUF_IN(II)   = PREEL(IR  ,JF)
  FOUBUF_IN(II+1) = PREEL(IR+1,JF)
ENDDO

END SUBROUTINE FOURIER_OUT_2D_NEW

SUBROUTINE FOURIER_OUT_SETUP()

USE PARKIND1,     ONLY : JPIM
USE TPM_DISTR,    ONLY : D, MYSETW
USE TPM_GEOMETRY, ONLY : G

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JM, IGLG, IPROC, ISTA, NMENMAX, KGL

NMENMAX=MAXVAL(G%NMEN(:))

IF (.not. ALLOCATED (ISTAIDX) ) THEN
  ALLOCATE(ISTAIDX(0:NMENMAX,D%NDGL_FS))
ENDIF

DO KGL=1,D%NDGL_FS
  ! Determine global latitude index corresponding to local latitude index KGL
  IGLG = D%NPTRLS(MYSETW) + KGL - 1

  ! Loop over all zonal wavenumbers relevant for this latitude
  DO JM = 0, G%NMEN(IGLG)
    ! Get the member of the W-set responsible for this zonal wavenumber in the "m" representation
    IPROC = D%NPROCM(JM)

    ! Compute offset for insertion of the fields in the l-to-m transposition buffer, FOUBUF_IN
    ISTA = (D%NSTAGT1B(D%MSTABF(IPROC)) + D%NPNTGTB0(JM,KGL)) * 2 

    ISTAIDX(JM,KGL)=ISTA

  ENDDO
ENDDO

END SUBROUTINE FOURIER_OUT_SETUP

END MODULE FOURIER_OUT_MOD
