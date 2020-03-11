PROGRAM WSMTD_rw_2D
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddif1, dummy2,v, &
             gridmin2, gridmax2, griddif2,num,num2, num3, &
             gridmin3, gridmax3, griddif3,&
             den,alpha,&
            kt,kt0,ktb,bias_fact, deltaT,&
            diff_s2,diff_s,ds2,ss,hh,dum,s1,s2,dummy11

REAL*8, ALLOCATABLE ::cv1(:,:),cv2(:,:),ht(:,:),vbias(:,:),ct(:,:),hill(:,:),width(:,:),&
           prob(:,:),fes1(:),fes2(:),dfds(:,:),av_dfds1(:),pcons(:),&
kcons(:), norm(:)




INTEGER :: dummy1,i,j,index1,index2,k,&
           t_min,t_max, &
           nbin1, nbin2,&
           i_mtd, i_md, i_s2, i_s1, ir, &
           mtd_max, nr, &
           narg,ios


INTEGER ::  w_cv,w_hill,mtd_steps, md_steps 
CHARACTER(LEN=50) :: filename_loc
CHARACTER(LEN=50), ALLOCATABLE :: filename(:), filename_mtd(:)

LOGICAL :: pmf,inpgrid,read_ct, read_vbias

CHARACTER*120 :: arg 
REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 
REAL*8, PARAMETER :: deg_to_rad = 4.d0*atan(1.d0)/180.d0
REAL*8, PARAMETER :: pi=4.d0*atan(1.d0)


 OPEN(11,FILE='COLVAR_1',STATUS='unknown')
 OPEN(12,FILE='HILLS_1',STATUS='unknown')
 CALL get_steps(11,md_steps)
 CALL get_steps(12,mtd_steps)
print *, "MD steps =", md_steps
print *, "MTD steps =", mtd_steps
close(11)
close(12)

kt0=300.d0
kt=2700.D0
deltaT=2700.D0
t_min=1
t_max=md_steps
pmf=.FALSE.
inpgrid=.false.
narg = IARGC()
nr=1
!md_steps=0
!mtd_steps=0
w_hill=0
w_cv=0
read_vbias=.false.
read_ct=.false.

DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-T0').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSE IF(INDEX(arg,'-dT').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)deltaT
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
  ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_max
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)gridmin1
      CALL GETARG(i+2,arg)
      READ(arg,*)gridmax1
      CALL GETARG(i+3,arg)
      READ(arg,*)griddif1
      CALL GETARG(i+4,arg)
      READ(arg,*)gridmin2
      CALL GETARG(i+5,arg)
      READ(arg,*)gridmax2
      CALL GETARG(i+6,arg)
      READ(arg,*)griddif2
      inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_hill
  ELSE IF(INDEX(arg,'-nr').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)nr
  ELSE IF(INDEX(arg,'-mtd_steps').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)mtd_steps
  ELSE IF(INDEX(arg,'-md_steps').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)md_steps
  ELSE IF(INDEX(arg,'-read_ct').NE.0)THEN
      read_ct=.true.
  ELSE IF(INDEX(arg,'-read_vbias').NE.0)THEN
      read_vbias=.true.
  END IF
END DO

bias_fact=(kt0+deltaT)/kt0  !gamma factor in PLUMED manual

if(nr.eq.1)print *, 'Warning! nr=1 is taken!'
IF(.NOT.inpgrid) print *,"error grids have to be mentioned in the input!"

if(md_steps.eq.0)stop 'md_steps cannot be zero'
if(mtd_steps.eq.0)stop 'mtd_steps cannot be zero'
if(w_hill.eq.0)stop 'w_hill cannot be zero'
if(w_cv.eq.0)stop 'w_cv cannot be zero'

WRITE(*,'(A,F9.2)')'System Temp (K)        =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'DeltaT Factor (K)      =',deltaT
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
WRITE(*,'(A,I10)')'No of replicas          =',nr

nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1

WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'BLM  COORD:', &
               gridmin1,gridmax1,griddif1,nbin1
WRITE(*,'(A10,3F8.4,I10)')'MTD COORD:', &
               gridmin2,gridmax2,griddif2,nbin2

if(nbin1.eq.0.or.nbin2.eq.0) stop 'number of bins cannot be zero'
if(nr.ne.nbin1)stop 'nr and nbin1 have to be the same'

ALLOCATE(pcons(nr))
ALLOCATE(kcons(nr))
ALLOCATE(filename(nr))
ALLOCATE(filename_mtd(nr))
open(44,file='replica.inp',status='old')
do ir=1,nr
  read(44,*)pcons(ir),kcons(ir)
  kcons(ir)=kcons(ir)*kj_to_kcal
  read(44,'(a)')filename(ir)
  read(44,'(a)')filename_mtd(ir)
end do
write(*,'(A10,2A16,2A20)')'#replica','rest_mean',&
'rest_k''CV_VAL_file','MTD_file'

do ir=1,nr
!  pcons(ir)= gridmin1+DFLOAT(ir-1)*griddif1
  write(*,'(i10,2f16.6,2A30)') ir,pcons(ir), kcons(ir), filename(ir),&
filename_mtd(ir)
end do


alpha=(kt0+deltaT)/deltaT
!alpha=(kt+deltaT)/deltaT   !USE THIS IF TEMP keyword is used with METAD keyword
!in plumed.dat
kt=kb*kt             ! kB T is calculed

!ktb=(alpha-1.d0)/kt  ! (alpha-1)/kT factor for c(t) computation

!alpha=bias_fact/(bias_fact-1.D0) !this is gamma in the paper   !using
!kt0 as PLUMED uses kt0 instead of kt
!bias_fact=(kt+bias_fact)/kt ! using correct bias_fact for c(t)
!calculation

print *, "NOTE: alpha parameter         = ", alpha
print *, "NOTE: kB T0 (kcal/mol)        = ", kt
print *, "NOTE: (alpha-1)/kT            = ", ktb

ALLOCATE(cv1(nr,md_steps),cv2(nr,md_steps),dfds(nr,md_steps))
DO ir=1,nr
  open(11,file=filename(ir),status='old')
  DO i_md=1,md_steps
! time, Ah, Bs, Cg, Dc, Hb, Rg, Rmsd
     READ(11,*)dummy11,dummy11,dummy11,dummy11,dummy11,dummy11,cv2(ir,i_md),cv1(ir,i_md) 
     diff_s=cv1(ir,i_md)-pcons(ir)  !TODO pcons should be in radians
     dfds(ir,i_md)=-diff_s*kcons(ir) !TODO kcons should be in kcal/mol
     !TODO in the above, the sign is not clear; need to change if
     !necessary
  END DO
  close(11)
END DO
      

ALLOCATE(ht(nr,mtd_steps),ct(nr,mtd_steps),hill(nr,mtd_steps),width(nr,mtd_steps))
do ir=1,nr
  OPEN(12,file=filename_mtd(ir),status='old')
  DO i_mtd=1,mtd_steps
     !READ(12,*) dummy11,ht(ir,i_mtd), width(ir,i_mtd),hill(ir,i_mtd)
     READ(12,*) dummy11,hill(ir,i_mtd), width(ir,i_mtd),ht(ir,i_mtd)
     ht(ir,i_mtd)=ht(ir,i_mtd)*kj_to_kcal
  END DO
  close(12)
END DO

if(read_ct)then
open(77,file='ct.rst',form='unformatted')
print *, '...reading ct.rst'
do i_mtd=1,mtd_steps
 read(77,iostat=ios)ct(1:nr,i_mtd) 
 if(ios.ne.0) then
    print *, 'error reading ct.rst ; i_mtd=', i_mtd, '  nr=',nr
     stop
 end if
end do
close(77)
else
!calculate c(t)
ALLOCATE(fes1(nbin2))
DO ir=1,nr
  call get_filename('ct.dat_',filename_loc,ir)
  OPEN(21,FILE=filename_loc,STATUS='unknown')
  fes1=0.d0
  DO i_mtd=1,mtd_steps
    ds2=width(ir,i_mtd)*width(ir,i_mtd)
    ss=hill(ir,i_mtd)
    hh=ht(ir,i_mtd)/alpha
    num=0.D0
    den=0.D0
    DO i_s2=1,nbin2
       diff_s=gridmin2+DFLOAT(i_s2-1)*griddif2-ss
       diff_s2=0.5d0*diff_s**2
       fes1(i_s2)=fes1(i_s2)+hh*DEXP(-diff_s2/ds2)
       num=num+DEXP(alpha*fes1(i_s2)/kt)
       den=den+DEXP((alpha-1.d0)*fes1(i_s2)/kt)
    END DO
    ct(ir,i_mtd)=kt*DLOG(num/den)
    WRITE(21,'(I10,F16.8)')i_mtd,ct(ir,i_mtd)
  END DO
  CLOSE(21)
  WRITE(*,'(A)')'CV values written in ',filename_loc
END DO
DEALLOCATE(fes1)

open(77,file='ct.rst',form='unformatted')
do i_mtd=1,mtd_steps
 write(77)ct(1:nr,i_mtd) 
end do
close(77)
end if
print *, "allocating..vbias"

allocate(vbias(nr,md_steps))
if(read_vbias)then
open(77,file='vbias.rst',form='unformatted')
do i_md=1,md_steps
 read(77)vbias(1:nr,i_md) 
end do
close(77)
else
print *, "computing vbias"
!calculate v(s,t)
DO ir=1,nr
   DO i_md=1,md_steps
     mtd_max=(i_md*w_cv/w_hill)+1
     dum=0.d0
     DO i_mtd=1,mtd_max-1 !till the previous hill added till the current
!time
           ds2=width(ir,i_mtd)*width(ir,i_mtd)
           hh=ht(ir,i_mtd)/alpha
           diff_s=cv2(ir,i_md)-hill(ir,i_mtd) 
           diff_s2=diff_s*diff_s*0.5D0
           dum=dum+hh*DEXP(-diff_s2/ds2)
     END DO
     vbias(ir,i_md)=dum
   END DO
   print *, "done...vbias ir=", ir
END DO

open(77,file='vbias.rst',form='unformatted')
do i_md=1,md_steps
 write(77)vbias(1:nr,i_md) 
end do
close(77)
end if

DEALLOCATE(hill)
print *, 'hill deallocated'


!calculate reweighted mean force < dF/ds exp(+beta [V^b (s2,t)-c(t)] ) >
!bug <totaling should be done within each replica>
!den=0.d0
!num=0.d0
ALLOCATE(av_dfds1(nr))
DO ir=1,nr
!bug <totaling should be done within each replica>
  den=0.d0
  num=0.d0
!  DO i_md=1,md_steps
  DO i_md=t_min,t_max
     i_mtd=(i_md*w_cv/w_hill) + 1
!since the bias is only felt when md_step is greater than first mtd
!step, the following is added !bug
     if(i_mtd*w_cv.gt.w_hill)then  !bug
       dum=vbias(ir,i_md) - ct(ir,i_mtd)
       num=num+dfds(ir,i_md)*DEXP(dum/kt)
       den=den+DEXP(dum/kt)
!  write(44,*) ir, i_md, dfds(ir,i_md), dum
     end if !bug
  END DO
  av_dfds1(ir)=num/den
!  write(43,*) ir,av_dfds1(ir),num,den
END DO
print *, 'av_dfds1 is computed'
open(42,file="av_dfds.dat")
do ir=1,nr
  write(42,*)pcons(ir), av_dfds1(ir)
end do
close(42)

DEALLOCATE(dfds)


!calculate projected free energy along s1
open(42,file="free_energy_1.dat")
ALLOCATE(fes1(nr))
ALLOCATE(fes2(nr))
fes1(1)=0.d0
num=0.d0
DO ir=1,nr-1
  dum=pcons(ir+1)-pcons(ir)
!  if (dum .gt. dmax ) dum =dum - drange  
!  if (dum .lt. dmin ) dum =dum + drange   
    num=num+dum* &
          (av_dfds1(ir+1)+av_dfds1(ir))
  fes1(ir)=num*0.5d0
  WRITE(42,*) pcons(ir+1),fes1(ir)
END DO
print *, 'fes1 computed'

DEALLOCATE(pcons)
DEALLOCATE(av_dfds1)

print *, 'deallocated pcons avdfds1'

ALLOCATE(prob(nbin1,nbin2))
allocate(norm(nr))
!calculate prob (unbiased from MTD potential)       
prob=0.d0
DO ir=1,nr
   norm(ir)=0.d0 !normalizing each slices of probability
   den=0.d0
   DO i_md=1,md_steps
      IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
        index1 = nint((cv1(ir,i_md)-gridmin1)/griddif1) +1
        index2 = nint((cv2(ir,i_md)-gridmin2)/griddif2) +1
        if(index1.gt.0.and.index2.gt.0.and. &
                  index1.le.nbin1.and.index2.le.nbin2)then
           !i_mtd=(i_md*w_cv/w_hill) + 1  !bug
           i_mtd=(i_md*w_cv/w_hill)  
           if(i_md*w_cv.gt.w_hill)then !bug <only after the first mtd
!step, the bias is experienced>
             dum=vbias(ir,i_md) - ct(ir,i_mtd)
             prob(index1,index2)=prob(index1,index2)+DEXP(dum/kt)
             den=den+DEXP(dum/kt)
           end if !bug
        end if
      END IF
   END DO
   norm(ir)=1.d0/(den*griddif1*griddif2)
END DO
print *, 'prob computed'

DEALLOCATE(cv1)
DEALLOCATE(cv2)
DEALLOCATE(vbias)
DEALLOCATE(ht)
DEALLOCATE(ct)

OPEN(2,FILE='free_energy.dat',STATUS='unknown')
DO i_s1=1,nbin1
   s1=DFLOAT(i_s1-1)*griddif1+gridmin1
   DO i_s2=1,nbin2
      s2=DFLOAT(i_s2-1)*griddif2+gridmin2
      num=prob(i_s1,i_s2)*norm(i_s1) !TODO: here it is assumed that i_s1
!and i_r are correctly tallied
      WRITE(2,'(5E16.8)')s1,s2,-kt*dlog(max(num,1E-32))+fes1(i_s1-1)
                               
                                   
                        
   END DO
   WRITE(2,*)
END DO
DEALLOCATE(prob)
DEALLOCATE(fes1)
DEALLOCATE(fes2)

WRITE(*,'(A)')'Free-energy and Unbiased distribution are in free_energy.dat'

END PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
      SUBROUTINE get_steps(iunit,nsteps)
      IMPLICIT NONE
      INTEGER iunit, nsteps
      INTEGER ios1,ios2
      nsteps=0
      REWIND(iunit)
      Read_Loop: DO
         READ(iunit,*,IOSTAT=ios1)
         IF(ios1.ne.0)EXIT Read_Loop
!         READ(iunit,*,IOSTAT=ios2)
!         IF(ios2.ne.0)stop 'error! expecting even no. of steps in
!         trajectory of cv'
         nsteps=nsteps+1
      END DO Read_Loop 
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE check_files(iunit,dt)
      IMPLICIT NONE
      INTEGER iunit, dt
      INTEGER ii, jj,i,ios
      dt=0
      i=2
      REWIND(iunit)
      READ(iunit,*)ii
      READ(iunit,*)jj
      dt=jj-ii
      ii=jj
      RLoop: DO 
        i=i+1
        READ(iunit,*,IOSTAT=ios)jj
        IF(ios.ne.0)EXIT RLoop
        IF(jj.ne.ii+dt)THEN
           print *, '!!ERROR: Steps are not at constant stride!!'
           print *, '!!       Unit No:',iunit,'!!'
           print *, '!!       Line No:',i,'!!'
           print *, '!! Expected stride =', dt,'!!'
           print *, '!! Actual stride =', jj-ii,'!!'
           STOP
        END IF
        ii=jj
      END DO RLoop
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1,&
                                gridmin2,gridmax2,griddif2)
      IMPLICIT NONE 
      INTEGER :: iunit 
      REAL*8  :: gridmin1, gridmax1, griddif1, & 
                gridmin2, gridmax2, griddif2
      INTEGER :: ii, ios
      REAL*8  :: cv1, cv2
      INTEGER, PARAMETER :: Def_Grid_Size=101
      REWIND(iunit)
      READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
      if(ios.ne.0)stop 'ERROR reading CV.dat'
      gridmin1=cv1
      gridmax1=cv1
      gridmin2=cv2
      gridmax2=cv2
      RLoop: DO 
        READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
        if(ios.ne.0)EXIT RLoop
        gridmin1=MIN(gridmin1,cv1)
        gridmin2=MIN(gridmin2,cv2)
        gridmax1=MAX(gridmax1,cv1)
        gridmax2=MAX(gridmax2,cv2)
      END DO RLoop
      griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
      griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
      END
!---------------------------------------------------------------------!
subroutine get_filename(head,filename,ir)
implicit none
character (len=*) :: head
character (len=*) :: filename
integer :: ir

 if(ir.lt.10)then
   write(filename,'(A,i1)')trim(head),ir
 else if(ir.lt.100)then
   write(filename,'(A,i2)')trim(head),ir
 else
   print *,'ERROR! get_filaname error: ir=',ir
   stop 
 end if
end subroutine get_filename
