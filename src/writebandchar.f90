!Written by A. D. N. James 

! !INTERFACE:
subroutine writebandchar
use modmain
implicit none
! local variables
integer lm,lmmax,irp,task_
integer is,ia,ias,ist,ispn,ik
!file extension
character(256) fileext
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:,:)

task_=task
if(task.eq.802) then
  task=20
  fileext='_BAND.OUT'
else 
  task=10
  fileext='.OUT'
endif
call init0
call init1
task=task_
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr

lmmax=(lmaxdos+1)**2
allocate(bc(lmmax,nspinor,natmtot,nstsv,nkpt))
call getbandchar(.false.,.false.,lmmax,bc)
!output the irep used to generate the characters
irp=0
if(lmirep) irp=1
!set up output file
open(50,file='BC'//trim(fileext),form='FORMATTED')
write(50,*) lmmax,nspinor,natmtot,nstsv,nkpt,irp
!write output 
do ik=1,nkpt
  do ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    do ispn=1,nspinor
      !output key for bc entries
      write(50,*) 
      write(50,*) ispn,ias,is,ia,ik
      do ist=1,nstsv
        !output nstsv \times lmmax bc matrix
        write(50,*) (dble(bc(lm,ispn,ias,ist,ik)), lm=1,lmmax)
      enddo
    enddo
  enddo
enddo
close(50)
write(*,*)
write(*,'("Info(writebandchar):")')
write(*,'(" Band characters of states written to BC.OUT")')
deallocate(bc)
return
end subroutine
!EOC

