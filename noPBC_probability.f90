PROGRAM probability
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddiff1, gridmin2, gridmax2, griddiff2
REAL*8 :: gridmin3, gridmax3, griddiff3, gridmin4, gridmax4, griddiff4
REAL*8 :: gridmin5, gridmax5, griddiff5, gridmin6, gridmax6, griddiff6
REAL*8 :: gridmin7, gridmax7, griddiff7
REAL*8 :: kt0, kt, ktb, bias_fact, alpha_cv, deltaT, den
REAL*8 :: diff_s2, ds2, ss, hh, dum, num, gamma_
REAL*8 :: s1, s2, s3, s4, s5, s6, s7

REAL*8, ALLOCATABLE :: prob(:,:), ct(:)
REAL*8, ALLOCATABLE :: cv1(:), cv2(:), cv5(:), vbias(:)
REAL*8, ALLOCATABLE :: cv3(:), cv4(:), cv6(:), cv7(:)

INTEGER :: mtd_steps, md_steps, i, t_min, t_max 
INTEGER :: mtd_max, w_cv, w_hill, i_mtd, i_md
INTEGER :: i_s1, i_s2, i_s3, i_s4, i_s5, i_s6, i_s7
INTEGER :: nbin1, nbin2, nbin3, nbin4, nbin5, nbin6, nbin7
INTEGER :: index1, index2, index3, index4, index5, index6, index7

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

OPEN(1,FILE='input',STATUS='old')
OPEN(2,FILE='1_PROB_2D.dat',STATUS='replace')
OPEN(11,FILE='COLVAR',STATUS='old')
OPEN(12,FILE='HILLS',STATUS='old')
OPEN(21,file='data_ct.dat',status='old')
OPEN(22,file='data_vbias.dat',status='old')


CALL get_steps(11,md_steps)
CALL get_steps(12,mtd_steps)

PRINT *, 'md_steps=', md_steps, 'mtd_steps=', mtd_steps

READ(1,*) kt0, kt, bias_fact
READ(1,*) t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
WRITE(*,*)' working for ', t_max, 'number of steps'

READ(1,*) gridmin1, gridmax1, nbin1
READ(1,*) gridmin2, gridmax2, nbin2
READ(1,*) gridmin3, gridmax3, nbin3
READ(1,*) gridmin4, gridmax4, nbin4
READ(1,*) gridmin5, gridmax5, nbin5
READ(1,*) gridmin6, gridmax6, nbin6
READ(1,*) gridmin7, gridmax7, nbin7
READ(1,*) w_cv, w_hill

! deltaT=(bias_fact - 1.d0)*kt
! deltaT = (bias_fact - 1.d0)*kt0
! alpha_cv = (kt + deltaT)/deltaT

kt = kb*kt
! gamma_ = 1/alpha
gamma_ = (bias_fact - 1.0)/bias_fact
WRITE(*,*) 'gamma_=', gamma_


WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max


ALLOCATE(cv1(md_steps),cv2(md_steps))
ALLOCATE(cv3(md_steps),cv4(md_steps))
ALLOCATE(cv5(md_steps),cv6(md_steps))
ALLOCATE(cv7(md_steps))
ALLOCATE(vbias(md_steps),ct(mtd_steps))


DO i_md=1,md_steps
   READ(11,*) dum,cv1(i_md),cv2(i_md),cv3(i_md),cv4(i_md),cv5(i_md), cv6(i_md), cv7(i_md)
                ! time, Ah, Bs, Cg, Dc, Hb, Rg, Rmsd
END DO

griddiff1 = (gridmax1-gridmin1)/dfloat(nbin1 -1)
griddiff2 = (gridmax2-gridmin2)/dfloat(nbin2 -1)
griddiff3 = (gridmax3-gridmin3)/dfloat(nbin3 -1)
griddiff4 = (gridmax4-gridmin4)/dfloat(nbin4 -1)
griddiff5 = (gridmax5-gridmin5)/dfloat(nbin5 -1)
griddiff6 = (gridmax6-gridmin6)/dfloat(nbin6 -1)
griddiff7 = (gridmax7-gridmin7)/dfloat(nbin7 -1)

!nbin = NINT((gridmax-gridmin)/griddif)+1 

WRITE(*,*) 'griddiff =', griddiff1, griddiff2, griddiff3, griddiff4, griddiff5, griddiff6, griddiff7

WRITE(*,*) "reading vbias.dat file"
DO i_md=1,md_steps
   READ(22,*) dum, vbias(i_md)
END DO

WRITE(*,*) "reading ct.dat file"
DO i_mtd=1,mtd_steps
   READ(21,*) dum, ct(i_mtd)
END DO

ALLOCATE(prob(nbin7,nbin1))

WRITE(*,*) 'calculating  probability'

den=0.d0
prob=0.d0

DO i_md=1,md_steps
   IF((i_md.GE.t_min).AND.(i_md.LE.t_max))THEN
      index7 = nint((cv7(i_md)-gridmin7)/griddiff7) +1
      index1 = nint((cv1(i_md)-gridmin1)/griddiff1) +1
      IF(index1.gt.0.and.index7.gt.0.and.index1.le.nbin1.and.index7.le.nbin7) then
      i_mtd=((i_md-1)*w_cv/w_hill)
      dum=vbias(i_md) - ct(i_mtd)
      prob(index7,index1) = prob(index7,index1) + DEXP(dum/kt)
      den=den+DEXP(dum/kt)
      END IF
   END IF
END DO

  den=den*griddiff1*griddiff7


DO i_s7=1,nbin7
s7=DFLOAT(i_s7-1)*griddiff7 + gridmin7
   DO i_s1=1,nbin1
   s1=DFLOAT(i_s1-1)*griddiff1 + gridmin1
       prob(i_s7,i_s1)=prob(i_s7,i_s1)*(1.d0/den)
         
       WRITE(2,*) s7, s1, prob(i_s7,i_s1)
    END DO
END DO
      WRITE(*,'(A)')'Unbiased distribution written in PROB_2D.dat'

close(1)
close(2)
close(11)
close(12)
close(21)
close(22)

DEALLOCATE(cv1, cv2, vbias, ct)
DEALLOCATE(cv3, cv4, cv5, cv6, cv7)
DEALLOCATE(prob)



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


