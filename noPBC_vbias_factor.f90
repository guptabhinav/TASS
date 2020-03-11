PROGRAM vbias_factor
IMPLICIT NONE
REAL*8 :: gridmin, gridmax, griddif
REAL*8 :: kt0, kt, ktb, bias_fact, alpha_cv, deltaT, den
REAL*8 :: diff_s2, ds2, ss, hh, dum, num, gamma_
REAL*8, ALLOCATABLE :: cv(:), ht(:), vbias(:), hill(:)
REAL*8, ALLOCATABLE ::  width(:)
INTEGER :: mtd_steps, md_steps, i, t_min, t_max 
INTEGER :: mtd_max, w_cv, w_hill, i_mtd, i_md
INTEGER :: i_s1, i_s2, i_s3, i_s4, nbin

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

OPEN(1,FILE='input',STATUS='old')
OPEN(11,FILE='COLVAR',STATUS='old')
OPEN(12,FILE='HILLS',STATUS='old')
OPEN(21,FILE='data_vbias.dat',STATUS='replace')


CALL get_steps(11,md_steps)
CALL get_steps(12,mtd_steps)

PRINT *, 'md_steps=', md_steps, 'mtd_steps=', mtd_steps

READ(1,*) kt0, kt, bias_fact
READ(1,*) t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'

READ(1,*) gridmin, gridmax, nbin 
READ(1,*) w_cv, w_hill

! deltaT=(bias_fact - 1.d0)*kt
! deltaT = (bias_fact - 1.d0)*kt0
! alpha_cv = (kt + deltaT)/deltaT


kt = kb*kt
! gamma_ = 1/alpha
gamma_ = (bias_fact - 1.0)/bias_fact
write(*,*) 'gamma_=', gamma_


WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max


ALLOCATE(cv(md_steps))
ALLOCATE(vbias(md_steps))
ALLOCATE(ht(mtd_steps))
ALLOCATE(hill(mtd_steps),width(mtd_steps))


DO i_md=1,md_steps
   READ(11,*) dum,dum,dum,dum,dum,dum,cv(i_md),dum
END DO

griddif = (gridmax-gridmin)/dfloat(nbin -1)
!nbin = NINT((gridmax-gridmin)/griddif)+1 

write(*,*) 'griddif =', griddif 

DO i_mtd=1,mtd_steps
   READ(12,*) dum,hill(i_mtd),width(i_mtd),ht(i_mtd)
   ht(i_mtd)=ht(i_mtd)*kj_to_kcal
END DO

WRITE(*,*) 'calculating vbias'

DO i_md=1,md_steps
   mtd_max=((i_md-1)*w_cv/w_hill)      ! bias is added at 11th,21st,... step
   ss=cv(i_md)
   dum=0.d0
   DO i_mtd=1,mtd_max
      ds2=width(i_mtd)*width(i_mtd)
      hh=ht(i_mtd)*gamma_
      diff_s2=ss-hill(i_mtd)
      diff_s2=diff_s2*diff_s2*0.5D0
      dum=dum+hh*DEXP(-diff_s2/ds2)
   END DO
   vbias(i_md)=dum
   WRITE(21,*) i_md, vbias(i_md)
END DO

write(*,*) " vbias written in data_vbias.dat file. "

close(1)
close(11)
close(21)
close(12)

DEALLOCATE(cv, ht, vbias, hill)
DEALLOCATE(width)

END PROGRAM

SUBROUTINE get_steps(iunit,nsteps)
 IMPLICIT NONE
 INTEGER :: iunit, nsteps
 INTEGER :: ios
  nsteps=0
   REWIND(iunit)
     Read_Loop: DO
       READ(iunit,*,IOSTAT=ios)
       IF(ios.ne.0)EXIT Read_Loop
       nsteps=nsteps+1
       END DO Read_Loop
   REWIND(iunit)
END SUBROUTINE





















