! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_TRANS

! Module to contain variables "local" to a specific call to a transform

!
USE PARKIND1  ,ONLY : JPIM,   JPRB

IMPLICIT NONE

SAVE

!INTEGER_M :: NF_UV      ! Number of u-v fields (spectral/fourier space)
!INTEGER_M :: NF_SCALARS ! Number of scalar fields (spectral/fourier space)
!INTEGER_M :: NF_SCDERS  ! Number of fields for derivatives of scalars
                        ! (inverse transform, spectral/fourier space)
!INTEGER_M :: NF_OUT_LT  ! Number of fields that comes out of Inverse
                        ! Legendre transform
INTEGER(KIND=JPIM) :: NF_SC2  ! Number of fields in "SPSC2" arrays.
INTEGER(KIND=JPIM) :: NF_SC3A ! Number of fields in "SPSC3A" arrays.
INTEGER(KIND=JPIM) :: NF_SC3B ! Number of fields in "SPSC3B" arrays.

!LOGICAL   :: LUV        ! uv fields requested
!LOGICAL   :: LSCALAR    ! scalar fields requested
LOGICAL   :: LVORGP     ! vorticity requested
LOGICAL   :: LDIVGP     ! divergence requested
LOGICAL   :: LUVDER     ! E-W derivatives of U and V requested
LOGICAL   :: LSCDERS    ! derivatives of scalar variables are req.
LOGICAL   :: LATLON     ! lat-lon output requested

!INTEGER_M :: NLEI2 ! 8*NF_UV + 2*NF_SCALARS + 2*NF_SCDERS (dimension in
                   ! inverse  Legendre transform)
!INTEGER_M :: NLED2 ! 2*NF_FS (dimension in direct Legendre transform)

!INTEGER_M :: NF_FS    ! Total number of fields in Fourier space

!INTEGER_M :: NF_GP        ! Total number of field in grid-point space
!INTEGER_M :: NF_UV_G      ! Global version of NF_UV (grid-point space)
!INTEGER_M :: NF_SCALARS_G ! Global version of NF_SCALARS (grid-point space)

REAL(KIND=JPRB), ALLOCATABLE :: FOUBUF_IN(:)  ! Fourier buffer
REAL(KIND=JPRB), ALLOCATABLE :: FOUBUF(:)     ! Fourier buffer

INTEGER(KIND=JPIM) :: NPROMA  ! Blocking factor for gridpoint input/output
INTEGER(KIND=JPIM) :: NGPBLKS ! Number of NPROMA blocks

LOGICAL :: LGPNORM = .FALSE.  ! indicates whether transform is being done for gpnorm

INTERFACE REALLOCATE_ARRAY
   MODULE PROCEDURE REALLOCATE_ARRAY_JPRB_R1, REALLOCATE_ARRAY_JPRB_R2
END INTERFACE REALLOCATE_ARRAY

INTERFACE REALLOCATE_ARRAY_PAD64BYTE
   MODULE PROCEDURE REALLOCATE_ARRAY_JPRB_PAD64BYTE_R1, REALLOCATE_ARRAY_JPRB_PAD64BYTE_R2
END INTERFACE REALLOCATE_ARRAY_PAD64BYTE

CONTAINS

  SUBROUTINE REALLOCATE_ARRAY_JPRB_R1(ARRAY, DIM1, INIT_VALUE)
    IMPLICIT NONE
    REAL(KIND=JPRB), ALLOCATABLE, INTENT(INOUT) :: ARRAY(:)
    INTEGER, INTENT(IN) :: DIM1
    REAL(KIND=JPRB), OPTIONAL :: INIT_VALUE
    INTEGER :: JI
    LOGICAL :: INIT_NEEDED

    INIT_NEEDED = .FALSE.
    
    IF (.NOT. ALLOCATED(ARRAY)) THEN
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1))
       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    ELSE IF (SIZE(ARRAY,1) < DIM1) THEN
       DEALLOCATE(ARRAY)
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1))
       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    END IF

   IF (INIT_NEEDED) THEN
      !$omp parallel do
      DO JI=1,DIM1
         ARRAY(JI)=INIT_VALUE
      END DO
      !$omp end parallel do 
   END IF
    
  END SUBROUTINE REALLOCATE_ARRAY_JPRB_R1
  
  SUBROUTINE REALLOCATE_ARRAY_JPRB_R2(ARRAY, DIM1, DIM2, INIT_VALUE)
    IMPLICIT NONE
    REAL(KIND=JPRB), ALLOCATABLE, INTENT(INOUT) :: ARRAY(:,:)
    INTEGER, INTENT(IN) :: DIM1, DIM2
    REAL(KIND=JPRB), OPTIONAL :: INIT_VALUE
    INTEGER :: JI, JJ
    LOGICAL :: INIT_NEEDED

    INIT_NEEDED = .FALSE.
    
    IF (.NOT. ALLOCATED(ARRAY)) THEN
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1, DIM2))
       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    ELSE IF (SIZE(ARRAY,1) < DIM1 .OR. SIZE(ARRAY,2) < DIM2) THEN
       DEALLOCATE(ARRAY)
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1, DIM2))
       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    END IF
    IF (INIT_NEEDED) THEN
       !$omp parallel do private(JI)
       DO JJ=1,DIM2
          DO JI=1,DIM1
             ARRAY(JI,JJ)=INIT_VALUE
          END DO
       END DO
       !$omp end parallel do 
    END IF
    
  END SUBROUTINE REALLOCATE_ARRAY_JPRB_R2

  FUNCTION PAD64BYTE_JPRB(NDIM) RESULT(NDIM_PAD)
    USE ISO_FORTRAN_ENV, ONLY : INT64
    INTEGER, INTENT(IN) ::NDIM
    INTEGER :: NDIM_PAD
    
    INTEGER, PARAMETER :: ELSIZE = STORAGE_SIZE(REAL(1,JPRB))/8
    INTEGER(KIND=INT64) :: NBYTES
    
    NBYTES = NDIM * ELSIZE ! Compute number of bytes
    NDIM_PAD = (( (NBYTES + 64 - 1) / 64) * 64 ) / ELSIZE ! Pad to next 64 bytes in elements
  END FUNCTION PAD64BYTE_JPRB
  
  SUBROUTINE REALLOCATE_ARRAY_JPRB_PAD64BYTE_R1(ARRAY, DIM1, INIT_VALUE)
    IMPLICIT NONE
    REAL(KIND=JPRB), ALLOCATABLE, INTENT(INOUT) :: ARRAY(:)
    INTEGER, INTENT(IN) :: DIM1
    REAL(KIND=JPRB), OPTIONAL :: INIT_VALUE
    
    INTEGER :: DIM1_PAD
    INTEGER :: JI
    LOGICAL :: INIT_NEEDED

    INIT_NEEDED = .FALSE.
    DIM1_PAD = PAD64BYTE_JPRB(DIM1)
    
    IF (.NOT. ALLOCATED(ARRAY)) THEN
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1_PAD))

       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    ELSE IF (SIZE(ARRAY,1) < DIM1_PAD) THEN
       DEALLOCATE(ARRAY)
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1_PAD))

       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    END IF

    IF (INIT_NEEDED) THEN
       !$omp parallel do private(JI)
       DO JI=1,DIM1_PAD
          ARRAY(JI)=INIT_VALUE
       END DO
       !$omp end parallel do 
    END IF    
  END SUBROUTINE REALLOCATE_ARRAY_JPRB_PAD64BYTE_R1
  
  SUBROUTINE REALLOCATE_ARRAY_JPRB_PAD64BYTE_R2(ARRAY, DIM1, DIM2, INIT_VALUE)
    USE ISO_FORTRAN_ENV, ONLY : INT64
    IMPLICIT NONE
    REAL(KIND=JPRB), ALLOCATABLE, INTENT(INOUT) :: ARRAY(:,:)
    INTEGER, INTENT(IN) :: DIM1, DIM2
    REAL(KIND=JPRB), OPTIONAL :: INIT_VALUE
    
    INTEGER :: DIM1_PAD
    INTEGER :: JI, JJ
    LOGICAL :: INIT_NEEDED

    INIT_NEEDED = .FALSE.
    DIM1_PAD = PAD64BYTE_JPRB(DIM1)
    
    IF (.NOT. ALLOCATED(ARRAY)) THEN
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1_PAD, DIM2))
       
       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    ELSE IF (SIZE(ARRAY,1) < DIM1_PAD .OR. SIZE(ARRAY,2) < DIM2) THEN
       DEALLOCATE(ARRAY)
#ifdef NO_OMP_ALLOCATOR_SUPPORT
#else
       !$omp allocators allocate(align(64): array)
#endif
       ALLOCATE(ARRAY(DIM1_PAD, DIM2))

       IF (PRESENT(INIT_VALUE)) INIT_NEEDED = .TRUE.
    END IF

    IF (INIT_NEEDED) THEN
       !$omp parallel do private(JI)
       DO JJ=1,DIM2
          DO JI=1,DIM1
             ARRAY(JI,JJ)=INIT_VALUE
          END DO
       END DO
       !$omp end parallel do 
    END IF
  END SUBROUTINE REALLOCATE_ARRAY_JPRB_PAD64BYTE_R2
    

END MODULE TPM_TRANS
