!Written by A. D. N. James

subroutine wanolp(ias,nst,l,lmmax,u,wfmt,wantemp)
!Calculates the overlap matrix elements between the wavfunctions and (apw)
!radial function
use modmain
implicit none
!arguments
integer, intent(in) :: ias
integer, intent(in) :: nst
integer, intent(in) :: l
integer, intent(in) :: lmmax
real(8), intent(in) :: u(nrmtmax)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(out) :: wantemp(lmmax,nst,nspinor)
! local variables
integer is,iro,ir,jr,nr,nri,nro,npi,ispn
integer i,idx,lm,m,ist,p
complex(8) z1
!correlated orb ias index
is=idxis(ias)
iro=nrmti(is)+lradstp
if (lradstp.eq.1) then
  nr=nrmt(is)
  nri=nrmti(is)
  npi=npmti(is)
else 
  nr=nrcmt(is)
  nri=nrcmti(is)
  npi=npcmti(is)
end if
nro=nr-nri
!Calculate the temporary Wannier matrix elements
do ispn=1,nspinor
  do ist=1,nst
    idx=1
    do lm=l**2+1,(l+1)**2
      z1=zzero
      i=npi+lm
! calculating the overlap between wfmt and the (apw) radial function
      if (l.le.lmaxi) then
        p=lm
        jr=1
        do ir=1,nri
          z1=z1+wrmt(jr,is)*u(jr)* &
               wfmt(p,ias,ispn,ist)
          p=p+lmmaxi
          jr=jr+lradstp
        end do
      end if
      jr=iro
      p=i
      do ir=nri+1,nr
        z1=z1+wrmt(jr,is)*u(jr)* &
             wfmt(p,ias,ispn,ist)
        p=p+lmmaxo
        jr=jr+lradstp
      end do
! put overlap into temporary projector array   
      wantemp(idx,ist,ispn)=z1
      idx=idx+1
    end do
  end do
end do
return
end subroutine
!EOC

