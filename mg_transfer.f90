!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_transfer 
!***********************************************************************
!                                                                      !
!  Transfer data between analysis and filter grid                      !
!                                                                      !
! Modules: kinds, mg_parameter, mg_intstate, mg_bocos, mg_interpolate, !
!          mg_timers, mg_mppstuff                                      !
!                                                     M. Rancic (2021) !
!***********************************************************************
use kinds, only: r_kind,i_kind
use mg_parameter
use mg_intstate, only: VALL,WORKA
use mg_intstate   , only: cvf1,cvf2,cvf3,cvf4,lref
use mg_interpolate, only: lsqr_adjoint_xyk,lsqr_forward_xyk
use mg_timers
use mg_mppstuff, only:  mype,ierror,mpi_comm_world
use mg_mppstuff, only: nx,my,ns,ms,ninc,minc,ninc2,minc2,mpi_comm_comp
use mg_mppstuff, only: finishMPI

implicit none
integer(i_kind):: n,m,l,k,i,j

public anal_to_filt_all
public filt_to_anal_all

public stack_to_composite
public composite_to_stack

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt_all
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_adjoint_all,lwq_vertical_adjoint
use mg_bocos, only:  bocoT_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: WORKB
real(r_kind),allocatable,dimension(:,:,:):: VFILT

!----------------------------------------------------------------------
!
! Apply vertical adjoint interpolation if needed
!
    allocate(WORKB(lmf_all,0:nm,0:mm))                      


     if(lmf_all==lm_all) then
        WORKB = WORKA
     else
        do k=1,km3
            call lwq_vertical_adjoint (lm,lmf,0,nm,0,mm,cvf1,cvf2,cvf3,cvf4,lref &
                                      ,WORKA((k-1)*lm+1:k*lm,:,:)                &
                                      ,WORKB((k-1)*lmf+1:k*lmf,:,:) )
        enddo
        do k=1,km2
          WORKB(km3*lmf+k,:,:) = WORKA(km3*lm+k,:,:)     
        enddo
     endif

!TEST
!     call finishMPI
!TEST

    allocate(VFILT(lmf_all,-ib:im+ib,-jb:jm+jb))                      


!T                                                 call btim(  aintp_tim)

      VFILT=0.
         call lsqr_adjoint_all(WORKB,VFILT,lmf_all)


    deallocate(WORKB)
!T                                                 call etim(  aintp_tim)

!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    

         call bocoT_2d(VFILT,lmf_all,im,jm,ib,jb,Fimax,Fjmax)
 
       VALL=0.
       VALL(1:lmf_all,0:im,0:jm)=VFILT(1:lmf_all,0:im,0:jm)
      

    deallocate(VFILT)

!                                            call etim(   btrns1_tim)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal_all
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: 
use mg_interpolate, only: lsqr_forward_all,lwq_vertical_direct
use mg_bocos, only:  boco_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: WORKB
real(r_kind),allocatable,dimension(:,:,:):: VFILT


!----------------------------------------------------------------------

!T                                            call btim(   btrns2_tim)

!***
!***  Define VFILT
!***

    allocate(VFILT(1:lmf_all,-ib:im+ib,-jb:jm+jb))                     

      VFILT=0.
      VFILT(1:lmf_all,0:im,0:jm)=VALL(1:lmf_all,0:im,0:jm)
        

!***
!***  Supply boundary conditions for VFILT
!***
         call boco_2d(VFILT,lmf_all,im,jm,ib,jb,Fimax,Fjmax)


!***
!*** Interpolate to analysis grid composite variables
!***

    allocate(WORKB(lmf_all,0:nm,0:mm))                      

!T                                                 call btim(   intp_tim)

         call lsqr_forward_all(VFILT,WORKB,lmf_all)

!                                                 call etim(   intp_tim)
    deallocate(VFILT)

     if(lmf_all==lm_all) then
        WORKA = WORKB
     else
        do k=1,km3
            call lwq_vertical_direct(lmf,lm,0,nm,0,mm,cvf1,cvf2,cvf3,cvf4,lref  &
                                    ,WORKB((k-1)*lmf+1:k*lmf,:,:)               &
                                    ,WORKA((k-1)*lm+1:k*lm,:,:) )
        enddo
        do k=1,km2
          WORKA(km3*lm+k,:,:)   =  WORKB(km3*lmf+k,:,:)  
        enddo
     endif


     deallocate(WORKB)



!                                                 call etim(   btrns2_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite                   &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D,lmax,lmax_all)
!----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: lmax,lmax_all
real(r_kind),dimension(lmax_all,-hx:im+hx,-hy:jm+hy),intent(in):: ARR_ALL
real(r_kind),dimension(km3,-hx:im+hx,-hy:jm+hy,lmax),intent(out):: A3D
real(r_kind),dimension(km2,-hx:im+hx,-hy:jm+hy)     ,intent(out):: A2D
!----------------------------------------------------------------------
    do L=1,lmax
      do j=-hy,jm+hy
      do i=-hx,im+hx
        do k=1,km3
          A3D(k,i,j,L)=ARR_ALL((k-1)*lmax+L,i,j)
        enddo
      enddo
      enddo
    enddo

   do k=1,km2
    A2D(k,:,:)=ARR_ALL(km3*lmax+k,:,:)
   enddo

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL,lmax,lmax_all)
!----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: lmax,lmax_all
real(r_kind),dimension(km2     ,-hx:im+hx,-hy:jm+hy     ),intent(in):: A2D
real(r_kind),dimension(km3     ,-hx:im+hx,-hy:jm+hy,lmax),intent(in):: A3D
real(r_kind),dimension(lmax_all,-hx:im+hx,-hy:jm+hy),intent(out):: ARR_ALL
integer(i_kind):: i,j,L
!----------------------------------------------------------------------
    do L=1,lmax
      do j=-hy,jm+hy
      do i=-hx,im+hx
        do k=1,km3
          ARR_ALL((k-1)*lmax+L,i,j)= A3D(k,i,j,L)
        enddo
      enddo
      enddo
    enddo

      do k=1,km2
        ARR_ALL(km3*lmax+k,:,:)= A2D(k,:,:)
      enddo

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_transfer
