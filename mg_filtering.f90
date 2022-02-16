!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_filtering
!***********************************************************************
!                                                                      !
! Contains all multigrid filtering prodecures                          ! 
!                                                                      ! 
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds, only: r_kind,i_kind
use mg_parameter, only: im,jm,hx,hy,hz,km2,km3,gm,Fimax,Fjmax
use mg_parameter, only: lm,lm_all,lmf,lmf_all,lmh,lmh_all
use mg_parameter, only: mgbf_line
!use mpimod, only: mype,ierror
use mg_mppstuff, only: mype,ierror
use mg_mppstuff, only: l_hgen,my_hgen,finishMPI,barrierMPI
use mg_generations, only: upsending_all,downsending_all,differencing_all
use mg_transfer, only: stack_to_composite,composite_to_stack
use mg_bocos, only: boco_2d
use mg_bocos, only: bocoT_2d
use mg_bocos, only: boco_3d
use mg_bocos, only: bocoT_3d
use mg_bocos, only: bocox,bocoy
use mg_bocos, only: bocoTx,bocoTy
use jp_pbfil, only: rbeta,rbetaT
use jp_pbfil3, only: dibetat,dibeta
use mg_output

use mg_timers

public mg_filtering_procedure 

private mg_filtering_rad1
private mg_filtering_rad2
private mg_filtering_rad3
private mg_filtering_lin1       
private mg_filtering_lin2
private mg_filtering_lin3
private mg_filtering_fast

private sup_vrbeta1
private sup_vrbeta1T
private sup_vrbeta3
private sup_vrbeta3T

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure(mg_filt) 
!***********************************************************************
!                                                                      !
! Driver for Multigrid filtering procedures with Helmholtz operator    !
!                                                                      !
!   1, 2, 3: Radial filter                                             !
!         1: 2d radial filter for all variables                        !
! ->      2: 2d radial filter with 1d in vertical for 3d variables     !
!         3: 3d radial filter for 3d variables                         !
!                                                                      !
!   4, 5, 6: Line filter                                               !
!         4: 2d line filter for all variables                          !
!         5: 2d line filter with 1d in vertical for 3d variables       !
!         6: 3d line filter for 3d variables                           !
!                                                                      !
!                                                                      !
!***********************************************************************
implicit none 

integer(i_kind),intent(in):: mg_filt
!-----------------------------------------------------------------------
  if(mgbf_line) then
    if(mg_filt<4) then
       print*,'("Line filters have options 4-6")'
       stop
    endif
  else
    if(mg_filt>3) then
       print*,'("Radial filters have options 1-3")'
       stop
    endif
  endif 
      select case(mg_filt)
        case(1)
          call mg_filtering_rad1
        case(2)
          call mg_filtering_rad2
        case(3)
          call mg_filtering_rad3
        case(4)
          call mg_filtering_lin1
        case(5)
          call mg_filtering_lin2
        case(6)
          call mg_filtering_lin3
        case default
          call mg_filtering_fast          
       end select

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_rad1
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 1:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter only for all variables                        !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp2,ss2
use mg_intstate, only: VALL,HALL
implicit none

integer(i_kind) L,i,j,g
!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
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

      call rbetaT(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbetaT(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


        call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
        call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call differencing_all(VALL,HALL)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***

                                                 call btim( bfilt_tim)

      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call rbeta(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbeta(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_rad1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_rad2
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter + 1d vertical filter                          !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
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

      call rbetaT(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,VALL(:,:,:))
      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
  if(l_hgen)  then
      call rbetaT(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,HALL(:,:,:))
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
  endif

      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)
  if(l_hgen)  then
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
   endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


        call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
        call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


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

      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call rbeta(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,VALL(:,:,:))
      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
  if(l_hgen)  then
      call rbeta(lm_all,hx,0,im,hy,0,jm,pasp2,ss2,HALL(:,:,:))
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
  endif

      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)
  if(l_hgen)  then
      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
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
                        endsubroutine mg_filtering_rad2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_rad3
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d radial filter 
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
use mg_intstate, only: pasp2,pasp3,ss2,ss3
use mg_intstate, only: VALL,HALL
implicit none


real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D


integer(i_kind) L,i,j

!----------------------------------------------------------------------
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


!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!
      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
        call rbetaT(km2,hx,0,im,hy,0,jm,pasp2,ss2,VM2D)
        call sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)

    if(l_hgen) then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
        call rbetaT(km2,hx,0,im,hy,0,jm,pasp2,ss2,HM2D)
        call sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
    endif 

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

        call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
        call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                 call btim( weight_tim)

      call differencing_all(VALL,HALL)

                                                 call etim( weight_tim)


!***
!*** Apply Beta filter at all generations 
!***

                                                 call btim( bfilt_tim)

      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
        call rbeta(km2,hx,0,im,hy,0,jm,pasp2,ss2,VM2D(:,:,:))
        call sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)
  if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
        call rbeta(km2,hx,0,im,hy,0,jm,pasp2,ss2,HM2D(:,:,:))
        call sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add 
!***  Then zero high generations 
!***

                                                 call btim(   dnsend_tim)

       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)
deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_rad3   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_lin1
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 4:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d line filter only for all variables                          !
!                                                                      !
!***********************************************************************
use mg_parameter, only: nfil
use mg_intstate, only: dixs,diys,hss2
use mg_intstate, only: VALL,HALL
implicit none

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff
!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
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

       do icol=3,1,-1
         call dibetat(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)

!         if(ff) then
!           print*,'("Failure in dietat at stage",i3,)',icol
!           print*,'(" failure occures at grid points i,j=",3i5)',iout,jout
!         endif

         call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
       enddo

     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif


         call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call differencing_all(VALL,HALL)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***


                                                 call btim( bfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
       do icol=1,3
         call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
         call dibeta(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)

!         if(ff) then
!           print*,'("Failure in dietat at stage",i3,)',icol
!           print*,'(" failure occures at grid points i,j=",3i5)',iout,jout
!         endif

       enddo

     do icol=1,3
         call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)

!         if(ff) then
!           print*,'("Failure in dietat at stage",i3,)',icol
!           print*,'(" failure occures at grid points i,j=",3i5)',iout,jout
!         endif

       endif
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_lin2
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 5:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter + 1d vertical filter
!                                                                      !
!***********************************************************************
use mg_parameter, only: nfil
use mg_intstate, only: dixs,diys,hss2
use mg_intstate, only: VALL,HALL
use mg_intstate, only: pasp1,ss1
implicit none

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

!----------------------------------------------------------------------

allocate(VM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; VM3D=0.
allocate(VM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; VM2D=0.
allocate(HM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; HM3D=0.
allocate(HM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; HM2D=0.


!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
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
! Horizontal
!

       do icol=3,1,-1
         call dibetat(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
         call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
       enddo

     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
         call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo
!
! Vertical
!

      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
        call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)

    if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
        call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
    endif

        call bocoT_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
        call bocoT_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call differencing_all(VALL,HALL)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***


                                                 call btim( bfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontal
!
       do icol=1,3
         call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
         call dibeta(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
       enddo

     do icol=1,3
         call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(lm_all,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
     enddo
!
! Vertical
!

      call boco_2d(VALL,lm_all,im,jm,hx,hy,Fimax,Fjmax)
      call boco_2d(HALL,lm_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)


      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
        call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)

    if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
        call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
    endif


       call barrierMPI
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)


deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_lin3
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 6:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d line filter                      
!                                                                      !
!***********************************************************************
!TEST
use, intrinsic :: ieee_arithmetic
!TEST
use mg_parameter, only: nfil
use mg_intstate, only: dixs,diys,dizs,hss2,vpasp3
use mg_intstate, only: qcols,dixs3,diys3,dizs3
use mg_intstate, only: VALL,HALL
use jp_pkind2, only: fpi
implicit none

integer(i_kind) k,i,j,L
integer(i_kind) icol,iout,jout,lout
logical:: ff

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

real(r_kind), allocatable, dimension(:,:,:,:):: W
real(r_kind), allocatable, dimension(:,:,:,:):: H

integer(fpi), allocatable, dimension(:,:,:):: JCOL


allocate(VM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; VM3D=0.
allocate(VM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; VM2D=0.
allocate(HM3D(km3,-hx:im+hx,-hy:jm+hy,lm))                      ; HM3D=0.
allocate(HM2D(km2,-hx:im+hx,-hy:jm+hy   ))                      ; HM2D=0.

allocate(W(km3,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz))                 ; W=0.
allocate(H(km3,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz))                 ; H=0.

allocate(JCOL(0:im,0:jm,1:Lm))                                  ; JCOL=0

!-----------------------------------------------------------------------


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

!
! From single stack to composite variables
!

      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
    if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
    endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply adjoint filter to 2D variables first
!

       do icol=3,1,-1
         call dibetat(km2,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
         call bocoT_2d(VM2D,km2,im,jm,hx,hy,Fimax,Fjmax)
       enddo
    
     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(km2,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
         call bocoT_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo

!
! Create and apply adjoint filter to extended 3D variables
!

         W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)

       do icol=7,1,-1
         do L=1,hz
           W(:,:,:,1-L )=0.
           W(:,:,:,LM+L)=0.
         end do
         call dibetat(km3,-hx,0,im,im+hx, -hy,0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
         call bocoT_3d(W,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax)
       enddo

     if(l_hgen)  then
          H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
     endif

     do icol=7,1,-1
       if(l_hgen)  then
         do L=1,hz
           H(:,:,:,1-L )=0.
           H(:,:,:,LM+L)=0.
         end do

         call dibetat(km3,-hx,0,im,im+hx, -hy,0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
       endif
         call bocoT_3d(H,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax,2,gm)
     enddo

!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

      VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)

    if(l_hgen)  then
      HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
    endif
!TEST
!      do L = 1,lm_all
!      do j=-hy,jm+hy
!      do i=-hx,im+hx
!        if(.not.ieee_is_finite(HALL(L,i,j))) then
!           print *,'mype,L,i,j is infinite',mype,L,i,j
!        endif
!      enddo
!      enddo
!      enddo
!TEST


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call differencing_all(VALL,HALL)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***
!
! From single stacked to composite variables
!
                                                 call btim( bfilt_tim)

      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)
    if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
    endif



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply filter to 2D variables first
!
       do icol=1,3
         call boco_2d(VM2D,km2,im,jm,hx,hy,Fimax,Fjmax)
         call dibeta(km2,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
       enddo

     do icol=1,3
         call boco_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(km2,-hx,0,im,im+hx, -hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
     enddo

!
! Create and apply filter to extended 3D variables
!

           W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)
        do L=1,hz
          do j=-hy,jm+hy
          do i=-hx,im+hx
              W(:,i,j,1-L )=VM3D(:,i,j, 1+L)
              W(:,i,j,LM+L)=VM3D(:,i,j,LM-L)
          end do
          end do
        end do

       do icol=1,7
         call boco_3d(W,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax)
         call dibeta(km3,-hx,0,im,im+hx, -hy,0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
        enddo 

     if(l_hgen)  then
           H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
        do L=1,hz
          do j=-hy,jm+hy
          do i=-hx,im+hx
              H(:,i,j,1-L )=HM3D(:,i,j, 1+L)
              H(:,i,j,LM+L)=HM3D(:,i,j,LM-L)
          end do
          end do
        end do
     endif
       do icol=1,7
         call boco_3d(H,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax,2,gm)
         if(l_hgen)  then
         call dibeta(km3,-hx,0,im,im+hx, -hy,0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
         endif
       enddo

!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

      VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)

    if(l_hgen)  then
      HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
    endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------

deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

deallocate(W)
deallocate(H)

deallocate(JCOL)

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin3

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

allocate(VM3D(km3,-hx:im+hx,-hy:jm+hy,lmf))                      ; VM3D=0.
allocate(VM2D(km2,-hx:im+hx,-hy:jm+hy    ))                      ; VM2D=0.
allocate(HM3D(km3,-hx:im+hx,-hy:jm+hy,lmh))                      ; HM3D=0.
allocate(HM2D(km2,-hx:im+hx,-hy:jm+hy    ))                      ; HM2D=0.



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
      call rbetaT(lmf_all,hx,0,im,paspx,ssx,VALL(:,:,j))
    enddo
       call bocoTx(VALL,lmf_all,im,jm,hx,hy)

    do i=0,im
      call rbetaT(lmf_all,hy,0,jm,paspy,ssy,VALL(:,i,:))
    enddo
      call bocoTy(VALL,lmf_all,im,jm,hx,hy)

      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)

  if(l_hgen)  then
    do j=0,jm
      call rbetaT(lmh_all,hx,0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
      call bocoTx(HALL,lmh_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

  if(l_hgen)  then
    do i=0,im
      call rbetaT(lmh_all,hy,0,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
      call bocoTy(HALL,lmh_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

      call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
!
! Vertically 
!
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lmf,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)
  if(l_hgen)  then
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lmh,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
   endif


       call barrierMPI
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


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


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
! Horizonatally

     call bocox(VALL,lmf_all,im,jm,hx,hy)
    do j=0,jm
      call rbeta(lmf_all,hx,0,im,paspx,ssx,VALL(:,:,j))
    enddo

      call bocoy(VALL,lmf_all,im,jm,hx,hy)
    do i=0,im
      call rbeta(lmf_all,hy,0,jm,paspy,ssy,VALL(:,i,:))
    enddo

      call stack_to_composite(VALL,VM2D,VM3D,lmf,lmf_all)

      call bocox(HALL,lmh_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do j=0,jm
      call rbeta(lmh_all,hx,0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
      call bocoy(HALL,lmh_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do i=0,im
      call rbeta(lmh_all,hy,0,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
  if(l_hgen)  then
    call stack_to_composite(HALL,HM2D,HM3D,lmh,lmh_all)
  endif

!
! Vertically
!

      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lmf,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL,lmf,lmf_all)
  if(l_hgen)  then
      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lmh,pasp1,ss1,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL,lmh,lmh_all)
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                          *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=0,jm
        do i=0,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L)=W(:,1+L)
            W(:,LM+L)=W(:,LM-L)
          end do
             call rbeta(kmax,hz,1,lm,  pasp,ss,W)
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

  
!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1T                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1T                                         *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm,  pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=0,jm
        do i=0,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L )=W(:,1+L )
            W(:,LM+L)=W(:,LM-L)
          end do
             call rbetaT(kmax,hz,1,lm, pasp,ss,W)
!
! Apply adjoint at the edges of domain
!
          do L=1,hz
            W(:,1+L)=W(:,1+L)+W(:,1-L)
            W(:,LM-L)=W(:,LM-L)+W(:,LM+L)
          enddo
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,0:im,0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(0:im,0:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

          do L=1,Lm
          do j=-hy,jm+hy
          do i=-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=-hy,jm+hy
          do i=-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j,1+L )
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call rbeta(kmax,hx,0,im, hy,0,jm, hz,1,lm, pasp,ss,W)

  
          do l=1,Lm
          do j=0,jm
          do i=0,im
            V(:,i,j,L)=W(:,i,j,L)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3T                       &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss,V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,0:im,0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(0:im,0:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,l

!----------------------------------------------------------------------

          do L=1,Lm
          do j=-hy,jm+hy
          do i=-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=-hy,jm+hy
          do i=-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j, 1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call rbetaT(kmax,hx,0,im, hy,0,jm, hz,1,lm, pasp,ss,W)

!
! Apply adjoint at the edges of domain
!
        do L=1,hz
          do j=-hy,jm+hy
          do i=-hx,im+hx
              W(:,i,j,1+L )=W(:,i,j, 1+L)+W(:,i,j, 1-L)
              W(:,i,j,LM-L)=W(:,i,j,LM-L)+W(:,i,j,LM+L)
          end do
          end do
         end do
  
          do l=1,lm
          do j=0,jm
          do i=0,im
            V(:,i,j,l)=W(:,i,j,l)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3T

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_filtering
