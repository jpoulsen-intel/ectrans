! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_MOD
INTERFACE FTINV
   MODULE PROCEDURE FTINV_1D ! CALL FTINV_1D(ZGTF_TR(:,JF),KF_FS,IGL) ! inplace
   MODULE PROCEDURE FTINV_2D ! CALL FTINV(ZGTF,KF_FS,IGL) ! no inplace
END INTERFACE
CONTAINS
SUBROUTINE FTINV_2D(PREEL,KFIELDS,KGL)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFTW - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFTW        ,ONLY : EXEC_FFTW
USE TPM_DIM         ,ONLY : R

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT), CONTIGUOUS  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE

!     ------------------------------------------------------------------

ITYPE = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST
IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0

DO JJ=IST1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

IF (G%NLOEN(IGLG)>1) THEN
  IOFF=D%NSTAGTF(KGL)+1
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2

  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,PREEL)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_2D

SUBROUTINE FTINV_1D(PREEL,KFIELDS,KGL)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFTW - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFTW        ,ONLY : EXEC_FFTW
USE TPM_DIM         ,ONLY : R

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT), CONTIGUOUS  :: PREEL(:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,JJ,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE

!     ------------------------------------------------------------------

ITYPE = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST
IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0

DO JJ=IST1,ILEN
   PREEL(IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
END DO

IF (G%NLOEN(IGLG)>1) THEN
  IOFF=D%NSTAGTF(KGL)+1
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2

  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,PREEL) 
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_1D

!     ------------------------------------------------------------------
END MODULE FTINV_MOD
