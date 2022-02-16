!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_fast
!***********************************************************************
!                                                                      !
! Fast multigrid filtering procedure:                                  !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 1d+1d horizontal filter + 1d vertical filter                   !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,paspx,paspy,ss1,ssx,ssy
use mg_intstate, only: VALL,HALL
implicit none

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

integer(i_kind) L,i,j
!-----------------------------------------------------------------------

allocate(VM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; VM3D=0.
allocate(VM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; VM2D=0.
allocate(HM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; HM3D=0.
allocate(HM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; HM2D=0.



!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontally
!

    do j=0,jm
      call rbetaT(lm_all,hx,0,im,paspx,ssx,VALL(:,:,j))
    enddo
       call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
!T       call bocoTx_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)


    do i=0,im
      call rbetaT(lm_all,hy,0,jm,paspy,ssy,VALL(:,i,:))
    enddo
      call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
!T      call bocoTy_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)

      call stack_to_composite(VALL,VM2D,VM3D)

  if(l_hgen)  then
    do j=0,jm
      call rbetaT(lm_all,hx,0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
      call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
!T      call bocoTx_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

  if(l_hgen)  then
    do i=0,im
      call rbetaT(lm_all,hy,0,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
      call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
!T      call bocoTy_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
!
! Vertically 
!
      call stack_to_composite(HALL,HM2D,HM3D)
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL)
   endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!        call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
!        call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call differencing_all(VALL,HALL)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( bfilt_tim)

!      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
!      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
! Horizonatally

!T      call bocox_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
    do j=0,jm
      call rbeta(lm_all,hx,0,im,paspx,ssx,VALL(:,:,j))
    enddo

!T      call bocoy_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
    do i=0,im
      call rbeta(lm_all,hy,0,im,paspx,ssy,VALL(:,i,:))
    enddo
      call stack_to_composite(VALL,VM2D,VM3D)

!T      call bocox_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do j=0,jm
      call rbeta(lm_all,hx,0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
!T      call bocoy_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do i=0,im
      call rbeta(lm_all,hx,0,im,paspx,ssx,HALL(:,i,:))
    enddo
  endif
  if(l_hgen)  then
    call stack_to_composite(HALL,HM2D,HM3D)
  endif

!
! Vertically
!

      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL)
   endif
       call barrierMPI


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)

deallocate(VM3D) 
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_fast
