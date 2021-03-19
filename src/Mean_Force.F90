!--------------------------------------------------------------------------------
! THIS SUBROUTINE COMPUTE'S MEAN FORCE AT EACH UMBRELLA WINDOW AND THEN 
! INTEGRATE THAT MEAN FORCE TO CALCULATE FREE ENERGY (1D) ALONG UMBRELLA CV
! u --> UMBRELLA CV INDEX ; nr --> NUMBER OF UMBRELLA WINDOWS
! pcons --> POSITION OF UMBRELL ; kcons --> Kappa VALUE IN EACH UMBRELLA WINDOW
! dfds = delF/delS ; fes --> FREE ENERGY ALONG UMBRELLA CV
! WRITTEN BY : Rahul Verma
!--------------------------------------------------------------------------------
MODULE MeanForce
USE GetSteps
CONTAINS
SUBROUTINE mean_force(u,ncv,nr,kt,gridmin,gridmax,griddif,nbin,t_min,t_max,pcons,fes)
IMPLICIT NONE
INTEGER                 :: i,j,i_md,dummy1,n,t_min,t_max,i_s1,i_s2,ir,nr,ios,narg
INTEGER                 :: ncv,w_cv,w_hill,md_steps,mtd_steps,indx1
INTEGER                 :: indx,nbin(*),u
REAL*8                  :: diff_s,den,num,kt,dum
REAL*8                  :: gridmin(*),gridmax(*),griddif(*)
REAL*8,ALLOCATABLE      :: cv(:,:,:),dummy(:,:,:),ct(:,:),prob(:)
REAL*8,ALLOCATABLE      :: dfds(:,:),av_dfds(:),pcons(:),kcons(:),norm(:),fes(:)
REAL*8, PARAMETER       :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER       :: au_to_kcal = 627.51
REAL*8, PARAMETER       :: kj_to_kcal = 0.239006

LOGICAL                 :: pmf,inpgrid,read_ct, read_vbias
CHARACTER(LEN=50)       :: filename_loc
CHARACTER(LEN=50),ALLOCATABLE :: filename(:),filename_mtd(:)

kt = kb*kt
IF(nbin(1) .eq. 0 ) STOP "ERROR : NUMBER OF BINS CAN NOT BE ZERO"

ALLOCATE(pcons(nr))
ALLOCATE(kcons(nr))
ALLOCATE(filename(nr))
ALLOCATE(filename_mtd(nr))

OPEN(10,FILE='replica.inp',STATUS='old',IOSTAT=ios)

IF (ios .lt. 0) STOP "ERROR : replica.inp doesn't exist..!" 

DO ir = 1,nr
   READ(10,*)pcons(ir),kcons(ir)
   kcons(ir)=kcons(ir)*au_to_kcal
   READ(10,'(A)')filename(ir)
ENDDO

101 FORMAT (A8,1X,A14,4X,A20,1X,A20,20X,A10,1X,A8,1X,A8)
102 FORMAT (I4,5X,F8.2,11X,F10.2,10X,A30,10X,I8,2X,I8,2X,I8)
write(*,101)'#replica','umbrella_mean','umbrella_k(kcal/mol)','CV_VAL_file',"MD Steps","t_min","t_max"

ALLOCATE(dfds(nr,9999999))
ALLOCATE(cv(nr,ncv,9999999))
ALLOCATE(av_dfds(nr))

DO i = 1,nr
OPEN(11,FILE=filename(i),status='old')
IF (t_max .eq. 0) t_max=t_min
CALL get_steps(11,md_steps)
t_max=md_steps
WRITE(*,102)i, pcons(i), kcons(i), filename(i),md_steps,t_min,t_max
ALLOCATE(dummy(nr,ncv,md_steps))

DO i_md=1,md_steps
     READ(11,*)dummy1,dummy1,(dummy(i,j,i_md),j=1,ncv),(cv(i,j,i_md),j=1,ncv)
!     WRITE(*,*)dummy1,dummy1,(dummy(i,j,i_md),j=1,ncv),(cv(i,j,i_md),j=1,ncv)
     diff_s = cv(i,u,i_md) - pcons(i)
     dfds(i,i_md) = -diff_s*kcons(i)
ENDDO
den = 0.d0 ; num = 0.d0 ; dum = 0.d0
   DO i_md = t_min,t_max
     num = num + dfds(i,i_md)
     den = den + DEXP(dum/kt)
   ENDDO
   av_dfds(i) = num/den
DEALLOCATE(dummy)
ENDDO
PRINT*,"av_dfds is computed"

OPEN(12,FILE='av_dfds.dat')
DO i = 1,nr
WRITE(12,*)pcons(i),av_dfds(i)
ENDDO

ALLOCATE(prob(nbin(1)))
ALLOCATE(norm(nr))

!!prob = 0.d0
!!DO i = 1,nr
!!   norm(i) = 0.d0 ; den = 0.d0
!!   DO i_md = 1,md_steps
!!     IF((i_md .gt. t_min) .and. (i_md .lt. t_max)) THEN
!!       indx = NINT((cv(i,1,i_md) - gridmin(1))/griddif(1)) + 1
!!         IF((indx .gt. 0) .and. (indx .le. nbin(1))) THEN
!!         prob(indx) = prob(indx) + 1.0
!!         ENDIF
!!     ENDIF
!!   ENDDO
!!      DO indx1 = 1,nbin(1)
!!        den = den + prob(indx1)
!!      ENDDO
!!norm(i) = den*griddif(1)
!!DO indx1 = 1,nbin(1)
!!WRITE(22,*)pcons(i),norm(i),prob(indx1)/norm(i)
!!ENDDO
!!ENDDO
!!PRINT*,"Probability Computed"

OPEN(13,FILE="free_energy.dat")
ALLOCATE(fes(nr))
fes(1) = 0.d0 ; num = 0.d0

DO i = 1,nr-1
   dum = pcons(i+1) - pcons(i)
   num = num + dum*(av_dfds(i+1) + av_dfds(i))
   fes(i+1) = num*0.5d0
   WRITE(13,*)pcons(i+1),fes(i+1)
ENDDO
PRINT*,"Free Energy is Computed"
!DEALLOCATE(fes)
DEALLOCATE(av_dfds)
DEALLOCATE(dfds)
CLOSE(12)

DEALLOCATE(cv)
!DEALLOCATE(pcons,kcons)
DEALLOCATE(filename,filename_mtd)
END SUBROUTINE 
END MODULE MeanForce
