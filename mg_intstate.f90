!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_intstate
!***********************************************************************
!                                                                      !
! Contains declarations and allocations of internal state variables    !
! use for filtering                                                    !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds, only: r_kind,i_kind
use jp_pkind2, only: fpi
use mpimod, only: mype
!SA use mg_mppstuff, only: mype
use mg_parameter, only: im,jm,nh,hx,hy,pasp0,nm,mm,ib,jb,nb,mb
use mg_parameter, only: lm,hz,p,km,km2,km3,lm_all
use mg_parameter, only: lmf,lmh,lmf_all,lmh_all
!SA use mg_parameter, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
use berror, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
use mg_mppstuff, only: my_hgen,finishMPI,barrierMPI
use jp_pbfil,only: cholaspect
use jp_pbfil,only: getlinesum
use jp_pbfil3, only: inimomtab,t22_to_3,tritform,t33_to_6,hextform
!TEST
!use gridmod, only: lat1,lon1
!TEST
implicit none

real(r_kind), allocatable,dimension(:,:,:):: V
!
! Composite control variable on first generation o filter grid
!
real(r_kind), allocatable,dimension(:,:,:):: VALL
real(r_kind), allocatable,dimension(:,:,:):: HALL
!
! Composite control variable on high generations of filter grid
!
!
!FOR ADJOINT TEST
!
!real(r_kind), allocatable,dimension(:,:):: A
!real(r_kind), allocatable,dimension(:,:):: B
!real(r_kind), allocatable,dimension(:,:):: A0
!real(r_kind), allocatable,dimension(:,:):: B0
!
real(r_kind), allocatable,dimension(:,:,:,:):: a2_diff
real(r_kind), allocatable,dimension(:,:,:,:):: b2_diff
real(r_kind), allocatable,dimension(:,:,:,:,:):: a3_diff
real(r_kind), allocatable,dimension(:,:,:,:,:):: b3_diff
real(r_kind), allocatable,dimension(:,:,:):: a_all_diff_f
real(r_kind), allocatable,dimension(:,:,:):: a_all_diff_h
real(r_kind), allocatable,dimension(:,:,:):: b_all_diff_f
real(r_kind), allocatable,dimension(:,:,:):: b_all_diff_h

real(r_kind), allocatable,dimension(:,:):: p_eps
real(r_kind), allocatable,dimension(:,:):: p_del
real(r_kind), allocatable,dimension(:,:):: p_sig
real(r_kind), allocatable,dimension(:,:):: p_rho

real(r_kind), allocatable,dimension(:,:,:):: paspx
real(r_kind), allocatable,dimension(:,:,:):: paspy
real(r_kind), allocatable,dimension(:,:,:):: pasp1
real(r_kind), allocatable,dimension(:,:,:,:):: pasp2
real(r_kind), allocatable,dimension(:,:,:,:,:):: pasp3

real(r_kind), allocatable,dimension(:,:,:):: vpasp2
real(r_kind), allocatable,dimension(:,:,:):: hss2
real(r_kind), allocatable,dimension(:,:,:,:):: vpasp3
real(r_kind), allocatable,dimension(:,:,:,:):: hss3

real(r_kind), allocatable,dimension(:):: ssx
real(r_kind), allocatable,dimension(:):: ssy
real(r_kind), allocatable,dimension(:):: ss1
real(r_kind), allocatable,dimension(:,:):: ss2
real(r_kind), allocatable,dimension(:,:,:):: ss3

integer(fpi), allocatable,dimension(:,:,:):: dixs
integer(fpi), allocatable,dimension(:,:,:):: diys
integer(fpi), allocatable,dimension(:,:,:):: dizs

integer(fpi), allocatable,dimension(:,:,:,:):: dixs3
integer(fpi), allocatable,dimension(:,:,:,:):: diys3
integer(fpi), allocatable,dimension(:,:,:,:):: dizs3

integer(fpi), allocatable,dimension(:,:,:,:):: qcols

!real(r_kind), allocatable,dimension(:,:,:,:):: r_vol
!
!
! Composite stacked variable
!

real(r_kind), allocatable,dimension(:,:,:):: WORKA


integer(i_kind),allocatable,dimension(:):: iref,jref
integer(i_kind),allocatable,dimension(:):: Lref,Lref_h
real(r_kind),allocatable,dimension(:):: cvf1,cvf2,cvf3,cvf4
real(r_kind),allocatable,dimension(:):: cvh1,cvh2,cvh3,cvh4

real(r_kind),allocatable,dimension(:):: cx0,cx1,cx2,cx3
real(r_kind),allocatable,dimension(:):: cy0,cy1,cy2,cy3

real(r_kind),allocatable,dimension(:,:):: cf00,cf01,cf02,cf03           &
                                         ,cf10,cf11,cf12,cf13           &
                                         ,cf20,cf21,cf22,cf23           &
                                         ,cf30,cf31,cf32,cf33

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine allocate_mg_intstate
!***********************************************************************
!                                                                      !
! Allocate internal state variables                                    !
!                                                                      !
!***********************************************************************

allocate(V(-hx:im+hx,-hy:jm+hy,lm))                          ; V=0.
allocate(VALL(lmf_all,-hx:im+hx,-hy:jm+hy))                  ; VALL=0.
allocate(HALL(lmh_all,-hx:im+hx,-hy:jm+hy))                  ; HALL=0.

!
!FOR ADJOINT TEST
!
!allocate(A(-hx:im+hx,-hy:jm+hy))                              ; A=0.
!allocate(B(-hx:im+hx,-hy:jm+hy))                              ; B=0.
!allocate(A0(-hx:im+hx,-hy:jm+hy))                             ; A0=0.
!allocate(B0(-hx:im+hx,-hy:jm+hy))                             ; B0=0.



allocate(a2_diff(-hx:im+hx,-hy:jm+hy,1:km2,1:2))                ; a2_diff=0. 
allocate(b2_diff(-hx:im+hx,-hy:jm+hy,1:km2,1:2))                ; b2_diff=0. 

allocate(a3_diff(-hx:im+hx,-hy:jm+hy,1:lm,1:km3,1:2))           ; a3_diff=0. 
allocate(b3_diff(-hx:im+hx,-hy:jm+hy,1:lm,1:km3,1:2))           ; b3_diff=0. 

allocate(a_all_diff_f(lmf_all,-hx:im+hx,-hy:jm+hy))             ; a_all_diff_f=0. 
allocate(a_all_diff_h(lmh_all,-hx:im+hx,-hy:jm+hy))             ; a_all_diff_h=0. 
allocate(b_all_diff_f(lmf_all,-hx:im+hx,-hy:jm+hy))             ; b_all_diff_f=0. 
allocate(b_all_diff_h(lmf_all,-hx:im+hx,-hy:jm+hy))             ; b_all_diff_h=0. 

allocate(p_eps(-hx:im+hx,-hy:jm+hy))                            ; p_eps=0.
allocate(p_del(-hx:im+hx,-hy:jm+hy))                            ; p_del=0.
allocate(p_sig(-hx:im+hx,-hy:jm+hy))                            ; p_sig=0.
allocate(p_rho(-hx:im+hx,-hy:jm+hy))                            ; p_rho=0.
allocate(paspx(1,1,0:im))                                       ; paspx=0.
allocate(paspy(1,1,0:jm))                                       ; paspy=0.
allocate(pasp1(1,1,1:lm))                                       ; pasp1=0.
allocate(pasp2(2,2,0:im,0:jm))                                  ; pasp2=0.
allocate(pasp3(3,3,0:im,0:jm,1:lm))                             ; pasp3=0.

allocate(vpasp2(0:2,0:im,0:jm))                                 ; pasp2=0.
allocate(hss2(0:im,0:jm,1:3))                                   ; hss2= 0.

allocate(vpasp3(1:6,0:im,0:jm,1:lm))                            ; vpasp3= 0.
allocate(hss3(0:im,0:jm,1:lm,1:6))                              ; hss3= 0.

allocate(ssx(0:im))                                             ; ssx=0.
allocate(ssy(0:jm))                                             ; ssy=0.
allocate(ss1(1:lm))                                             ; ss1=0.
allocate(ss2(0:im,0:jm))                                        ; ss2=0.
allocate(ss3(0:im,0:jm,1:lm))                                   ; ss3=0.

allocate(dixs(0:im,0:jm,3))                                     ; dixs=0
allocate(diys(0:im,0:jm,3))                                     ; diys=0

allocate(dixs3(0:im,0:jm,1:lm,6))                               ; dixs3=0
allocate(diys3(0:im,0:jm,1:lm,6))                               ; diys3=0
allocate(dizs3(0:im,0:jm,1:lm,6))                               ; dizs3=0

allocate(qcols(0:7,0:im,0:jm,1:lm))                             ; qcols=0

!
! In stnadalone version
!
!allocate(r_vol(lm_all,0:nm,0:mm,2))                             ; r_vol=0.
!
! ... but in GSI there will be 
!     r_vol2 and r_vol3 for 2d and 3d variables
! and r_vol3 will need to be given vertical structure
!

!
allocate(WORKA(lm_all,0:nm,0:mm))                               ; WORKA=0.

!
! for re-decomposition
!

allocate(iref(0:nm))                                     ; iref=0
allocate(jref(0:mm))                                     ; jref=0

allocate(cx0(0:nm))                                      ; cx0=0.
allocate(cx1(0:nm))                                      ; cx1=0.
allocate(cx2(0:nm))                                      ; cx2=0.
allocate(cx3(0:nm))                                      ; cx3=0.

allocate(cy0(0:mm))                                      ; cy0=0.
allocate(cy1(0:mm))                                      ; cy1=0.
allocate(cy2(0:mm))                                      ; cy2=0.
allocate(cy3(0:mm))                                      ; cy3=0.

allocate(cf00(0:nm,0:mm))                            ; cf00=0.
allocate(cf01(0:nm,0:mm))                            ; cf01=0.
allocate(cf02(0:nm,0:mm))                            ; cf02=0.
allocate(cf03(0:nm,0:mm))                            ; cf03=0.
allocate(cf10(0:nm,0:mm))                            ; cf10=0.
allocate(cf11(0:nm,0:mm))                            ; cf11=0.
allocate(cf12(0:nm,0:mm))                            ; cf12=0.
allocate(cf13(0:nm,0:mm))                            ; cf13=0.
allocate(cf20(0:nm,0:mm))                            ; cf20=0.
allocate(cf21(0:nm,0:mm))                            ; cf21=0.
allocate(cf22(0:nm,0:mm))                            ; cf22=0.
allocate(cf23(0:nm,0:mm))                            ; cf23=0.
allocate(cf30(0:nm,0:mm))                            ; cf30=0.
allocate(cf31(0:nm,0:mm))                            ; cf31=0.
allocate(cf32(0:nm,0:mm))                            ; cf32=0.
allocate(cf33(0:nm,0:mm))                            ; cf33=0.

allocate(Lref(1:lm))                                 ; Lref=0
allocate(Lref_h(1:lmf))                              ; Lref_h=0

allocate(cvf1(1:lm))                                 ; cvf1=0.
allocate(cvf2(1:lm))                                 ; cvf2=0.
allocate(cvf3(1:lm))                                 ; cvf3=0.
allocate(cvf4(1:lm))                                 ; cvf4=0.

allocate(cvh1(1:lmf))                                ; cvh1=0.
allocate(cvh2(1:lmf))                                ; cvh2=0.
allocate(cvh3(1:lmf))                                ; cvh3=0.
allocate(cvh4(1:lmf))                                ; cvh4=0.


!-----------------------------------------------------------------------
                        endsubroutine allocate_mg_intstate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_mg_weights
!***********************************************************************
!                                                                      !
! Define weights and scales                                            !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L
real(r_kind):: gen_fac
!-----------------------------------------------------------------------

      p_eps(:,:)=0.0
      p_del(:,:)=0.0
      p_sig(:,:)=0.0
      p_rho(:,:)=0.0

!--------------------------------------------------------
      gen_fac=1.
      a2_diff(:,:,:,1)=mg_weig1 
      a3_diff(:,:,:,:,1)=mg_weig1 
      a_all_diff_f(:,:,:)=mg_weig1 
      a_all_diff_h(:,:,:)=mg_weig1 

      b2_diff(:,:,:,1)=0.
      b3_diff(:,:,:,:,1)=0.
      b_all_diff_f(:,:,:)=0.
      b_all_diff_h(:,:,:)=0.

!      r_vol(:,:,:,1)=1.


      select case(my_hgen)
        case(2) 
!          r_vol(:,:,:,2)=0.25             ! In standalone case
!          gen_fac=0.25
          a2_diff(:,:,:,2)=mg_weig2 
          a3_diff(:,:,:,:,2)=mg_weig2
          b2_diff(:,:,:,2)=0.
          b3_diff(:,:,:,:,2)=0.
          a_all_diff_h(:,:,:)=mg_weig2
          b_all_diff_h(:,:,:)=0.
        case(3) 
!          r_vol(:,:,:,2)=0.0625           ! In standalone case
!          gen_fac=0.0625
          a2_diff(:,:,:,2)=mg_weig3 
          a3_diff(:,:,:,:,2)=mg_weig3 
          b2_diff(:,:,:,2)=0
          b3_diff(:,:,:,:,2)=0
          a_all_diff_h(:,:,:)=mg_weig3 
          b_all_diff_h(:,:,:)=0.
        case default 
!          r_vol(:,:,:,2)=0.015625         ! In standalone case
!          gen_fac=0.015625
          a2_diff(:,:,:,2)=mg_weig4 
          a3_diff(:,:,:,:,2)=mg_weig4
          b2_diff(:,:,:,2)=0.
          b3_diff(:,:,:,:,2)=0.
          a_all_diff_h(:,:,:)=mg_weig4
          b_all_diff_h(:,:,:)=0.
      end select


          do L=1,lm
           pasp1(1,1,L)=pasp0
          enddo

          do i=0,im
            paspx(1,1,i)=pasp0
          enddo  
          do j=0,jm
            paspy(1,1,j)=pasp0
          enddo  

          do j=0,jm
          do i=0,im
            pasp2(1,1,i,j)=pasp0*(1.+p_del(i,j))
            pasp2(2,2,i,j)=pasp0*(1.-p_del(i,j))
            pasp2(1,2,i,j)=pasp0*p_eps(i,j)     
            pasp2(2,1,i,j)=pasp0*p_eps(i,j)     
          end do
          end do

        do L=1,lm
          do j=0,jm
          do i=0,im
            pasp3(1,1,i,j,l)=pasp0*(1+p_del(i,j))
            pasp3(2,2,i,j,l)=pasp0
            pasp3(3,3,i,j,l)=pasp0*(1-p_del(i,j))
            pasp3(1,2,i,j,l)=pasp0*p_eps(i,j)
            pasp3(2,1,i,j,l)=pasp0*p_eps(i,j)
            pasp3(2,3,i,j,l)=pasp0*p_sig(i,j)
            pasp3(3,2,i,j,l)=pasp0*p_sig(i,j)
            pasp3(1,3,i,j,l)=pasp0*p_rho(i,j)
            pasp3(3,1,i,j,l)=pasp0*p_rho(i,j)
          end do
          end do
        end do


        call cholaspect(1,lm,pasp1)
        call cholaspect(0,im,0,jm,pasp2)
        call cholaspect(0,im,0,jm,1,lm,pasp3)


        call getlinesum(hx,0,im,paspx,ssx)
        call getlinesum(hy,0,jm,paspy,ssy)
        call getlinesum(hz,1,lm,pasp1,ss1)
        call getlinesum(hx,0,im,hy,0,jm,pasp2,ss2)
        call getlinesum(hx,0,im,hy,0,jm,hz,1,lm,pasp3,ss3)
!-----------------------------------------------------------------------
                        endsubroutine def_mg_weights

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_line
!***********************************************************************
!                                                                      !
! Inititate line filters                                               !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L,icol
logical:: ff
!-----------------------------------------------------------------------

  do j=0,jm
  do i=0,im
    call t22_to_3(pasp2(:,:,i,j),vpasp2(:,i,j))
  enddo
  enddo

  do l=1,lm
  do j=0,jm
  do i=0,im
    call t33_to_6(pasp3(:,:,i,j,l),vpasp3(:,i,j,l))
  enddo
  enddo
  enddo



  call inimomtab(p,nh,ff)

  call tritform(0,im,0,jm,vpasp2, dixs,diys, ff)

  do icol=1,3
    hss2(:,:,icol)=vpasp2(icol-1,:,:)
  enddo  


  call hextform(0,im,0,jm,1,lm,vpasp3,qcols,dixs3,diys3,dizs3, ff)


  do icol=1,6
    hss3(:,:,:,icol)=vpasp3(icol,:,:,:)
  enddo
 

!-----------------------------------------------------------------------
                        endsubroutine init_mg_line

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine deallocate_mg_intstate
!***********************************************************************
!                                                                      !
! Deallocate internal state variables                                  !
!                                                                      !
!***********************************************************************

deallocate(V)

deallocate(HALL,VALL)

!FOR ADJOINT TEST
!
!deallocate(A,B,A0,B0)
!

deallocate(a2_diff,b2_diff)
deallocate(a3_diff,b3_diff)
deallocate(a_all_diff_f,b_all_diff_f)
deallocate(a_all_diff_h,b_all_diff_h)
deallocate(p_eps,p_del,p_sig,p_rho,pasp1,pasp2,pasp3,ss1,ss2,ss3)
deallocate(dixs,diys)
deallocate(dixs3,diys3,dizs3)
deallocate(qcols)
!
! for testing
!
deallocate(WORKA)

!
! for re-decomposition
!
deallocate(iref,jref)

deallocate(cf00,cf01,cf02,cf03,cf10,cf11,cf12,cf13)
deallocate(cf20,cf21,cf22,cf23,cf30,cf31,cf32,cf33)

deallocate(Lref,Lref_h)

deallocate(cvf1,cvf2,cvf3,cvf4)

deallocate(cvh1,cvh2,cvh3,cvh4)




!-----------------------------------------------------------------------
                        endsubroutine deallocate_mg_intstate


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_intstate
