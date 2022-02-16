!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_generations
!***********************************************************************
!                                                                      !
!  Contains procedures that include differrent generations             !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds, only: r_kind,i_kind
use mg_parameter, only: im,jm,imL,jmL,imH,jmH,hx,hy,gm
use mg_parameter, only: km,km2,km3,Fimax,Fjmax,FimaxL,FjmaxL
use mg_parameter, only: lm,lm_all,lmf,lmf_all,lmh,lmh_all
!use mpimod, only: mype
use mg_mppstuff, only: mype
use mg_mppstuff, only: my_hgen,l_hgen,barrierMPI,finishMPI,Fimax,Fjmax
use mg_bocos, only: upsend,downsend
use mg_bocos, only: boco_2d
use mg_bocos, only: bocoT_2d
use mg_bocos, only: upsend_all,downsend_all
use mg_intstate, only: a2_diff,b2_diff,a3_diff,b3_diff
use mg_intstate, only: a_all_diff_f,b_all_diff_f
use mg_intstate, only: a_all_diff_h,b_all_diff_h
use mg_intstate, only: cvh1,cvh2,cvh3,cvh4,lref_h
use mg_interpolate, only: lwq_vertical_adjoint,lwq_vertical_direct
use mg_timers
!TEST
use, intrinsic :: ieee_arithmetic
!TEST

interface forward_k
module procedure forward_k2d,forward_k3d
endinterface forward_k

interface adjoint_k
module procedure adjoint_k2d,adjoint_k3d
endinterface adjoint_k


interface differencing
module procedure differencing_2d,differencing_3d
endinterface differencing

public upsending_all
public downsending_all
public differencing_all

private forward_all
private adjoint_all

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine forward_k2d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations   (2d case)                  !
!                                                                      !
!***********************************************************************
(L02D,H2D,g)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-1:imL+1,-1:jmL+1,km2), intent(in):: L02D
real(r_kind), dimension(0:im,0:jm,km2), intent(out):: H2D
integer(i_kind), intent(in):: g
real(r_kind), dimension(0:im,-1:jmL+1,km2):: HXLY
integer(i_kind):: i,iL,j,jL
!-----------------------------------------------------------------------
 
 if(g==my_hgen.or.g==1) then

    do jL=-1,jmL+1
    do i=0,im,2
      iL=i/2
      HXLY(i,jL,:)=L02D(iL,jL,:)
    enddo
    enddo

    do jL=-1,jmL+1
    do i=1,im-1,2
      iL=(i+1)/2
      HXLY(i,jL,:)=-0.0625*L02D(iL-2,jL,:) &
                   +0.5625*L02D(iL-1,jL,:) &
                   +0.5625*L02D(iL  ,jL,:) &
                   -0.0625*L02D(iL+1,jL,:)
     enddo
     enddo

     do j=0,jm,2
       jL = j/2
     do i=0,im
       H2D(i,j,:)=HXLY(i,jL,:)
     enddo
     enddo 

     do j=1,jm-1,2
       jL=(j+1)/2
     do i=0,im
      H2D(i,j,:)=-0.0625*HXLY(i,jL-2,:) &
                 +0.5625*HXLY(i,jL-1,:) &
                 +0.5625*HXLY(i,jL  ,:) &
                 -0.0625*HXLY(i,jL+1,:)
     enddo
     enddo

 endif
!-----------------------------------------------------------------------
                        endsubroutine forward_k2d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine forward_k3d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(L02D,H2D,g)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-1:imL+1,-1:jmL+1,lm,km3), intent(in):: L02D
real(r_kind), dimension(0:im,0:jm,lm,km3), intent(out):: H2D
integer(i_kind),intent(in):: g
real(r_kind), dimension(0:im,-1:jmL+1,km3):: HXLY
integer(i_kind):: i,iL,j,jL,L
!-----------------------------------------------------------------------

 if(g==my_hgen.or.g==1) then

  do L=1,lm

    do jL=-1,jmL+1
    do i=0,im,2
      iL=i/2
      HXLY(i,jL,:)=L02D(iL,jL,L,:)
    enddo
    enddo

    do jL=-1,jmL+1
    do i=1,im-1,2
      iL=(i+1)/2
      HXLY(i,jL,:)=-0.0625*L02D(iL-2,jL,L,:) &
                   +0.5625*L02D(iL-1,jL,L,:) &
                   +0.5625*L02D(iL  ,jL,L,:) &
                   -0.0625*L02D(iL+1,jL,L,:)
     enddo
     enddo

     do j=0,jm,2
       jL = j/2
     do i=0,im
       H2D(i,j,L,:)=HXLY(i,jL,:)
     enddo
     enddo 

     do j=1,jm-1,2
       jL=(j+1)/2
     do i=0,im
      H2D(i,j,L,:)=-0.0625*HXLY(i,jL-2,:) &
                   +0.5625*HXLY(i,jL-1,:) &
                   +0.5625*HXLY(i,jL  ,:) &
                   -0.0625*HXLY(i,jL+1,:)
     enddo
     enddo

  enddo

 endif
!-----------------------------------------------------------------------
                        endsubroutine forward_k3d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine forward_all                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(L02D,H2D,g,lmax_all)
!-----------------------------------------------------------------------
implicit none
integer(i_kind):: lmax_all
real(r_kind), dimension(lmax_all,-1:imL+1,-1:jmL+1), intent(in):: L02D
real(r_kind), dimension(lmax_all,0:im,0:jm), intent(out):: H2D
integer(i_kind), intent(in):: g
real(r_kind), dimension(lmax_all,0:im,-1:jmL+1):: HXLY
integer(i_kind):: i,iL,j,jL,L
!-----------------------------------------------------------------------

!T if(g==my_hgen.or.g==1) then
!TEST
!      do jL=-1,jmL+1
!      do iL=-1,imL+1
!      do L = 1,lm_all
!        if(.not.ieee_is_finite(L02D(L,iL,jL))) then
!           print *,'mype,L,iL,jL is infinite',mype,L,iL,jL
!        endif
!      enddo
!      enddo
!      enddo
!TEST

    do jL=-1,jmL+1
    do i=0,im,2
      iL=i/2
      HXLY(:,i,jL)=L02D(:,iL,jL)
    enddo
    enddo

    do jL=-1,jmL+1
    do i=1,im-1,2
      iL=(i+1)/2
      HXLY(:,i,jL)=-0.0625*L02D(:,iL-2,jL) &
                   +0.5625*L02D(:,iL-1,jL) &
                   +0.5625*L02D(:,iL  ,jL) &
                   -0.0625*L02D(:,iL+1,jL)
     enddo
     enddo

     do j=0,jm,2
       jL = j/2
     do i=0,im
       H2D(:,i,j)=HXLY(:,i,jL)
     enddo
     enddo 

     do j=1,jm-1,2
       jL=(j+1)/2
     do i=0,im
      H2D(:,i,j)=-0.0625*HXLY(:,i,jL-2) &
                 +0.5625*HXLY(:,i,jL-1) &
                 +0.5625*HXLY(:,i,jL  ) &
                 -0.0625*HXLY(:,i,jL+1)
     enddo
     enddo


!T endif

!T call barrierMPI
!-----------------------------------------------------------------------
                        endsubroutine forward_all


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_k2d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations   (2d case)                  !
!                                                                      !
!***********************************************************************
(H02D,L2D,g)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:im,0:jm,km2), intent(in):: H02D
real(r_kind), dimension(-1:imL+1,-1:jmL+1,km2), intent(out):: L2D
integer(i_kind), intent(in):: g
real(r_kind), dimension(0:im,-1:jmL+1,km2):: HXLY
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------

  if(g==1 .or. g==my_hgen) then

      L2D(:,:,:) = 0.
      HXLY(:,:,:)= 0.

    do j=jm-1,1,-2
      jL = (j+1)/2
    do i=im,0,-1
       iL=i/2
      HXLY(i,jL-2,:)=HXLY(i,jL-2,:)-0.0625*H02D(i,j,:)
      HXLY(i,jL-1,:)=HXLY(i,jL-1,:)+0.5625*H02D(i,j,:)
      HXLY(i,jL  ,:)=HXLY(i,jL  ,:)+0.5625*H02D(i,j,:)
      HXLY(i,jL+1,:)=HXLY(i,jL+1,:)-0.0625*H02D(i,j,:)
    enddo
    enddo
  
    do j=jm,0,-2
      jL = j/2
    do i=im,0,-1
      HXLY(i,jL,:)=HXLY(i,jL,:)+H02D(i,j,:)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i=im-1,1,-2
      iL = (i+1)/2
      L2D(iL-2,jL,:)=L2D(iL-2,jL,:)-0.0625*HXLY(i,jL,:)
      L2D(iL-1,jL,:)=L2D(iL-1,jL,:)+0.5625*HXLY(i,jL,:)
      L2D(iL  ,jL,:)=L2D(iL  ,jL,:)+0.5625*HXLY(i,jL,:)
      L2D(iL+1,jL,:)=L2D(iL+1,jL,:)-0.0625*HXLY(i,jL,:)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i =im,0,-2
      iL = i/2
      L2D(iL,jL,:)=L2D(iL,jL,:)+HXLY(i,jL,:)
    enddo  
    enddo  

  endif
!-----------------------------------------------------------------------
                        endsubroutine adjoint_k2d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_k3d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(H02D,L2D,g)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:im,0:jm,lm,km3), intent(in):: H02D
real(r_kind), dimension(-1:imL+1,-1:jmL+1,lm,km3), intent(out):: L2D
integer(i_kind), intent(in):: g
real(r_kind), dimension(0:im,-1:jmL+1,km3):: HXLY
integer(i_kind):: i,j,iL,jL,l
!-----------------------------------------------------------------------

  if(g==1 .or. g==my_hgen) then

      L2D(:,:,:,:) = 0.

    Do l=1,lm

      HXLY(:,:,:)= 0.

    do j=jm-1,1,-2
      jL = (j+1)/2
    do i=im,0,-1
       iL=i/2
      HXLY(i,jL-2,:)=HXLY(i,jL-2,:)-0.0625*H02D(i,j,L,:)
      HXLY(i,jL-1,:)=HXLY(i,jL-1,:)+0.5625*H02D(i,j,L,:)
      HXLY(i,jL  ,:)=HXLY(i,jL  ,:)+0.5625*H02D(i,j,L,:)
      HXLY(i,jL+1,:)=HXLY(i,jL+1,:)-0.0625*H02D(i,j,L,:)
    enddo
    enddo
  
    do j=jm,0,-2
      jL = j/2
    do i=im,0,-1
      HXLY(i,jL,:)=HXLY(i,jL,:)+H02D(i,j,L,:)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i=im-1,1,-2
      iL = (i+1)/2
      L2D(iL-2,jL,L,:)=L2D(iL-2,jL,L,:)-0.0625*HXLY(i,jL,:)
      L2D(iL-1,jL,L,:)=L2D(iL-1,jL,L,:)+0.5625*HXLY(i,jL,:)
      L2D(iL  ,jL,L,:)=L2D(iL  ,jL,L,:)+0.5625*HXLY(i,jL,:)
      L2D(iL+1,jL,L,:)=L2D(iL+1,jL,L,:)-0.0625*HXLY(i,jL,:)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i =im,0,-2
      iL = i/2
      L2D(iL,jL,L,:)=L2D(iL,jL,L,:)+HXLY(i,jL,:)
    enddo  
    enddo  

    enddo

  endif
!-----------------------------------------------------------------------
                        endsubroutine adjoint_k3d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_all                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(H02D,L2D,g,lmax_all)
!-----------------------------------------------------------------------
implicit none
integer(i_kind):: lmax_all
real(r_kind), dimension(lmax_all,0:im,0:jm), intent(in):: H02D
real(r_kind), dimension(lmax_all,-1:imL+1,-1:jmL+1), intent(out):: L2D
integer(i_kind), intent(in):: g
real(r_kind), dimension(lmax_all,0:im,-1:jmL+1):: HXLY
integer(i_kind):: i,j,iL,jL,l
!-----------------------------------------------------------------------

  if(g==1 .or. g==my_hgen) then

      L2D(:,:,:) = 0.

      HXLY(:,:,:)= 0.

    do j=jm-1,1,-2
      jL = (j+1)/2
    do i=im,0,-1
       iL=i/2
      HXLY(:,i,jL-2)=HXLY(:,i,jL-2)-0.0625*H02D(:,i,j)
      HXLY(:,i,jL-1)=HXLY(:,i,jL-1)+0.5625*H02D(:,i,j)
      HXLY(:,i,jL  )=HXLY(:,i,jL  )+0.5625*H02D(:,i,j)
      HXLY(:,i,jL+1)=HXLY(:,i,jL+1)-0.0625*H02D(:,i,j)
    enddo
    enddo
  
    do j=jm,0,-2
      jL = j/2
    do i=im,0,-1
      HXLY(:,i,jL)=HXLY(:,i,jL)+H02D(:,i,j)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i=im-1,1,-2
      iL = (i+1)/2
      L2D(:,iL-2,jL)=L2D(:,iL-2,jL)-0.0625*HXLY(:,i,jL)
      L2D(:,iL-1,jL)=L2D(:,iL-1,jL)+0.5625*HXLY(:,i,jL)
      L2D(:,iL  ,jL)=L2D(:,iL  ,jL)+0.5625*HXLY(:,i,jL)
      L2D(:,iL+1,jL)=L2D(:,iL+1,jL)-0.0625*HXLY(:,i,jL)
    enddo
    enddo

    do jL=jmL+1,-1,-1
    do i =im,0,-2
      iL = i/2
      L2D(:,iL,jL)=L2D(:,iL,jL)+HXLY(:,i,jL)
    enddo  
    enddo  


  endif
!-----------------------------------------------------------------------
                        endsubroutine adjoint_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsending_all                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(lmf_all,-hx:im+hx,-hy:jm+hy),intent(in):: V
real(r_kind),dimension(lmh_all,-hx:im+hx,-hy:jm+hy),intent(out):: H

real(r_kind),dimension(lmf_all,-1:imL+1,-1:jmL+1):: V_INT
real(r_kind),dimension(lmh_all,-1:imL+1,-1:jmL+1):: H_INT
integer(i_kind):: g,L,k
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

!                                                        call btim(adj1_tim)
        call adjoint_all(V(1:lmf_all,0:im,0:jm),V_INT,1,lmf_all) 
!                                                        call etim(adj1_tim)
!                                                        call btim(boco1_tim)
        call bocoT_2d(V_INT,lmf_all,imL,jmL,1,1,FimaxL,FjmaxL)
!                                                        call etim(boco1_tim)
!Place for transition to lmh
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if(lmh==lmf) then
        H_INT = V_INT
      else
        do k=1,km3
             call lwq_vertical_adjoint(lmf,lmh,0,imL,0,jmL                    &
                                      ,cvh1,cvh2,cvh3,cvh4,lref_h             &
                                      ,V_INT((k-1)*lmf+1:k*lmf,0:imL,0:jmL)   &
                                      ,H_INT((k-1)*lmh+1:k*lmh,0:imL,0:jmL))
        enddo
        do k=1,km2
          H_INT(km3*lmh+k,0:imL,0:jmL) = V_INT(km3*lmf+k,0:imL,0:jmL)     
        enddo
      endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Place for transition to lmh
!                                                        call btim(upsend1_tim)
        call upsend_all(H_INT(1:lmh_all,0:imL,0:jmL),H,lmh_all)
!                                                        call etim(upsend1_tim)

! Reset value for H_INT

       H_INT = 0.
!
! From generation 2 sequentially to higher generations
!
      do g=2,gm-1

!                                                        call btim(adj2_tim)
        call adjoint_all(H(1:lmh_all,0:im,0:jm),H_INT,g,lmh_all) 
!                                                        call etim(adj2_tim)
!                                                        call btim(boco2_tim)
        call bocoT_2d(H_INT,lmh_all,imL,jmL,1,1,FimaxL,FjmaxL,g,g)
!                                                        call etim(boco2_tim)
!                                                        call btim(upsend2_tim)
        call upsend_all(H_INT(1:lmh_all,0:imL,0:jmL),H,lmh_all,g,g+1)
!                                                        call etim(upsend2_tim)

      end do 


!-----------------------------------------------------------------------
                        endsubroutine upsending_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsending_all                      &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(H,V)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(lmh_all,-hx:im+hx,-hy:jm+hy),intent(inout):: H
real(r_kind),dimension(lmf_all,-hx:im+hx,-hy:jm+hy),intent(inout):: V
real(r_kind),dimension(lmh_all,-1:imL+1,-1:jmL+1):: H_INT
real(r_kind),dimension(lmf_all,-1:imL+1,-1:jmL+1):: V_INT
real(r_kind),dimension(lmh_all,0:im,0:jm):: H_PROX
real(r_kind),dimension(lmf_all,0:im,0:jm):: V_PROX
integer(i_kind):: g,L
integer(i_kind):: iL,jL,i,j,k
!-----------------------------------------------------------------------


      do g=gm,3,-1

      if(my_hgen==g) then
         H_PROX(1:lmh_all,0:im,0:jm)= H(1:lmh_all,0:im,0:jm)
      endif
!T                                                        call btim(dnsend1_tim)
        call downsend_all(H_PROX,H_INT,lmh_all,g,g-1)
!T                                                        call etim(dnsend1_tim)

!T                                                        call btim(boco1_tim)
        call boco_2d(H_INT,lmh_all,imL,jmL,1,1,FimaxL,FjmaxL,g-1,g-1)
!T                                                        call etim(boco1_tim)

!T                                                        call btim(forw1_tim)
      if(my_hgen==g-1) then
        call forward_all(H_INT,H_PROX,g-1,lmh_all)
        H(1:lmh_all,0:im,0:jm)=H     (1:lmh_all,0:im,0:jm)            &
                              +H_PROX(1:lmh_all,0:im,0:jm)
      endif
!T                                                       call etim(forw1_tim)

      enddo



      if(my_hgen==2) then
         H_PROX(1:lmh_all,0:im,0:jm)= H(1:lmh_all,0:im,0:jm)
      endif

!T                                                        call btim(dnsend2_tim)
        call downsend_all(H_PROX,H_INT,lmh_all)

          H(:,:,:)=0.

!Place for transition to lmf
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if(lmf==lmh) then
          V_INT = H_INT
        else
          do k=1,km3
              call lwq_vertical_direct (lmh,lmf,0,iml,0,jml                   &
                                       ,cvh1,cvh2,cvh3,cvh4,lref_h            &
                                       ,H_INT((k-1)*lmh+1:k*lmh,0:imL,0:jmL)  &
                                       ,V_INT((k-1)*lmf+1:k*lmf,0:imL,0:jmL) )
          enddo
          do k=1,km2
            V_INT(km3*lmf+k,0:imL,0:jmL)   =  H_INT(km3*lmh+k,0:imL,0:jmL)  
          enddo
        endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Place for transition to lmf

!T                                                        call etim(dnsend2_tim)
!T                                                        call btim(boco2_tim)
        call boco_2d(V_INT,lmf_all,imL,jmL,1,1,FimaxL,FjmaxL)
!T                                                        call etim(boco2_tim)

!T                                                        call btim(forw2_tim)
        call forward_all(V_INT,V_PROX,1,lmf_all)
!T                                                        call etim(forw2_tim)

!T                                                        call btim(add2_tim)
          V(1:lmf_all,0:im,0:jm)=V     (1:lmf_all,0:im,0:jm)         &
                                +V_PROX(1:lmf_all,0:im,0:jm)
!T                                                        call etim(add2_tim)

!-----------------------------------------------------------------------
                        endsubroutine downsending_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing_2d                      &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator                                      !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(-hx:im+hx,-hy:jm+hy,km2),intent(inout):: V
real(r_kind),dimension(-hx:im+hx,-hy:jm+hy,km2),intent(inout):: H
real(r_kind),dimension(-1:im,0:jm, km2):: DIFX
real(r_kind),dimension( 0:im,-1:jm,km2):: DIFY
integer(i_kind):: i,j,k
!-----------------------------------------------------------------------


     do j=0,jm
     do i=-1,im
       DIFX(i,j,:)=V(i+1,j,:)-V(i,j,:)
     enddo
     enddo
     do j=-1,jm
     do i=0,im
       DIFY(i,j,:)=V(i,j+1,:)-V(i,j,:)
     enddo
     enddo
     do j=0,jm
     do i=0,im
       V(i,j,:)=a2_diff(i,j,:,1)*V(i,j,:)                               &
               -b2_diff(i,j,:,1)*(DIFX(i,j,:)-DIFX(i-1,j,:)             &
                                 +DIFY(i,j,:)-DIFY(i,j-1,:))
     enddo
     enddo

 if(l_hgen) then

     do j=0,jm
     do i=-1,im
       DIFX(i,j,:)=H(i+1,j,:)-H(i,j,:)
     enddo
     enddo
     do j=-1,jm
     do i=0,im
       DIFY(i,j,:)=H(i,j+1,:)-H(i,j,:)
     enddo
     enddo
     do j=0,jm
     do i=0,im
       H(i,j,:)=a2_diff(i,j,:,2)*H(i,j,:)                               &
               -b2_diff(i,j,:,2)*(DIFX(i,j,:)-DIFX(i-1,j,:)             &
                                 +DIFY(i,j,:)-DIFY(i,j-1,:))
     enddo
     enddo

 endif


!-----------------------------------------------------------------------
                        endsubroutine differencing_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing_3d                      &
!***********************************************************************
!                                                                      !
!  Apply 3D differential operator                                      !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(-hx:im+hx,-hy:jm+hy,lm,km3),intent(inout):: V
real(r_kind),dimension(-hx:im+hx,-hy:jm+hy,lm,km3),intent(inout):: H
real(r_kind),dimension(-1:im,0:jm ,lm  ,km3):: DIFX
real(r_kind),dimension(0:im ,-1:jm,lm  ,km3):: DIFY
real(r_kind),dimension(0:im ,0:jm ,0:lm,km3):: DIFZ
integer(i_kind):: i,j,l,k
!-----------------------------------------------------------------------

   do l=1,lm
     do j=0,jm
     do i=-1,im
       DIFX(i,j,l,:)=V(i+1,j,l,:)-V(i,j,l,:)
     enddo
     enddo
     do j=-1,jm
     do i=0,im
       DIFY(i,j,l,:)=V(i,j+1,l,:)-V(i,j,l,:)
     enddo
     enddo
   enddo

   do l=1,lm-1
     do j=0,jm
     do i=0,im
       DIFZ(i,j,l,:)=V(i,j,l+1,:)-V(i,j,l,:)
     enddo
     enddo
   enddo
     do j=0,jm
     do i=0,im
       DIFZ(i,j,0 ,:)=-DIFZ(i,j,1   ,:)
       DIFZ(i,j,lm,:)=-DIFZ(i,j,lm-1,:)
     enddo
     enddo

   do l=1,lm
     do j=0,jm
     do i=0,im
        V(i,j,l,:)=a3_diff(i,j,l,:,1)*V(i,j,l,:)                        &
                  -b3_diff(i,j,l,:,1)*(DIFX(i,j,l,:)-DIFX(i-1,j,l,:)    &
                                      +DIFY(i,j,l,:)-DIFY(i,j-1,l,:)    &
                                      +DIFZ(i,j,l,:)-DIFZ(i,j,l-1,:))
     enddo
     enddo
   enddo

if(l_hgen) then

   do l=1,lm
     do j=0,jm
     do i=-1,im
       DIFX(i,j,l,:)=H(i+1,j,l,:)-H(i,j,l,:)
     enddo
     enddo
     do j=-1,jm
     do i=0,im
       DIFY(i,j,l,:)=H(i,j+1,l,:)-H(i,j,l,:)
     enddo
     enddo
   enddo

   do l=1,lm-1
     do j=0,jm
     do i=0,im
       DIFZ(i,j,l,:)=H(i,j,l+1,:)-H(i,j,l,:)
     enddo
     enddo
   enddo
     do j=0,jm
     do i=0,im
       DIFZ(i,j,0 ,:)=-DIFZ(i,j,1   ,:)
       DIFZ(i,j,lm,:)=-DIFZ(i,j,lm-1,:)
     enddo
     enddo

   do l=1,lm
     do j=0,jm
     do i=0,im
        H(i,j,l,:)=a3_diff(i,j,l,:,2)*H(i,j,l,:)                       &
                  -b3_diff(i,j,l,:,2)*(DIFX(i,j,l,:)-DIFX(i-1,j,l,:)   &
                                      +DIFY(i,j,l,:)-DIFY(i,j-1,l,:)   &
                                      +DIFZ(i,j,l,:)-DIFZ(i,j,l-1,:))
     enddo
     enddo
   enddo

endif

!-----------------------------------------------------------------------
                        endsubroutine differencing_3d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing_all                     &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(lmf_all,-hx:im+hx,-hy:jm+hy),intent(inout):: V
real(r_kind),dimension(lmh_all,-hx:im+hx,-hy:jm+hy),intent(inout):: H
real(r_kind),dimension(lmf_all,-1:im, 0:jm):: DIFX
real(r_kind),dimension(lmf_all,0:im ,-1:jm):: DIFY
real(r_kind),dimension(lmh_all,-1:im, 0:jm):: DIFXH
real(r_kind),dimension(lmh_all,0:im ,-1:jm):: DIFYH
integer(i_kind):: i,j,l,k,imx,jmx
!-----------------------------------------------------------------------

     do j=0,jm
     do i=-1,im
       DIFX(:,i,j)=V(:,i+1,j)-V(:,i,j)
     enddo
     enddo
     do j=-1,jm
     do i=0,im
       DIFY(:,i,j)=V(:,i,j+1)-V(:,i,j)
     enddo
     enddo


     do j=0,jm
     do i=0,im
       V(:,i,j)=a_all_diff_f(:,i,j)*V(:,i,j)                      &
               -b_all_diff_f(:,i,j)*(DIFX(:,i,j)-DIFX(:,i-1,j)    &
                                    +DIFY(:,i,j)-DIFY(:,i,j-1))   
     enddo
     enddo

if(l_hgen) then

!  imx = Fimax(my_hgen)
!  jmx = Fjmax(my_hgen)

   imx = im
   jmx = jm

     do j=0,jmx
     do i=-1,imx
       DIFXH(:,i,j)=H(:,i+1,j)-H(:,i,j)
     enddo
     enddo
     do j=-1,jmx
     do i=0,imx
       DIFYH(:,i,j)=H(:,i,j+1)-H(:,i,j)
     enddo
     enddo

     do j=0,jmx
     do i=0,imx
        H(:,i,j)=a_all_diff_h(:,i,j)*H(:,i,j)                          &
                -b_all_diff_h(:,i,j)*(DIFXH(:,i,j)-DIFXH(:,i-1,j)      &
                                     +DIFYH(:,i,j)-DIFYH(:,i,j-1))  
     enddo
     enddo

endif

!-----------------------------------------------------------------------
                        endsubroutine differencing_all

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_generations
