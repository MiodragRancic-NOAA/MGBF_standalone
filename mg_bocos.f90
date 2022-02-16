!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_bocos
!***********************************************************************
!                                                                      !
!  Provide communication between subdomains and supply halos on        !
!  filter grid                                                         !
!                                                                      !
! Libraries: mpi                                                       !
! Modules: kinds, mg_mppstuff, mg_parameter, mg_domain                 !
!                                                     M. Rancic (2021) !
!***********************************************************************
use kinds, only: r_kind,i_kind
!use mpimod, only: mype,mpi_comm_world
use mg_mppstuff, only: mype,mpi_comm_world,mpi_comm_work
use mg_mppstuff, only: itype,rtype,dtype,mpi_comm_comp,my_hgen     &
                      ,barrierMPI,finishMPI,l_hgen,mype_hgen
use mg_parameter, only: gm

implicit none

interface boco_2d
  module procedure boco_2d_g1 
  module procedure boco_2d_gh 
endinterface

interface bocoT_2d
  module procedure bocoT_2d_g1 
  module procedure bocoT_2d_gh 
endinterface

public:: upsend       
public:: downsend     

interface upsend_all
  module procedure upsend_all_g1
  module procedure upsend_all_gh
endinterface

interface downsend_all     
  module procedure downsend_all_gh
  module procedure downsend_all_g2
endinterface

interface boco_3d
 module procedure boco_3d_g1
 module procedure boco_3d_gh
endinterface

interface bocoT_3d
 module procedure bocoT_3d_g1
 module procedure bocoT_3d_gh
endinterface


public:: multiply_add

private:: bocosH1
private:: boco05      
private:: v02v   

interface bocox
  module procedure bocox_2d_g1
  module procedure bocox_2d_gh
endinterface

interface bocoy
  module procedure bocoy_2d_g1
  module procedure bocoy_2d_gh
endinterface

interface bocoTx
  module procedure bocoTx_2d_g1
  module procedure bocoTx_2d_gh
endinterface

interface bocoTy
  module procedure bocoTy_2d_g1
  module procedure bocoTy_2d_gh
endinterface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine multiply_add                         &
!***********************************************************************
!                                                                      *
!         Adjoint test                                                 *
!                                                                      *
!***********************************************************************
(a,b,im,jm,hx,hy,res_glob)
!-----------------------------------------------------------------------
use mpi
use mg_domain, only: Fleast,Flwest,Flsouth,Flnorth

implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: im,jm,hx,hy
real(r_kind), dimension(-hx:im+hx,-hy:jm+hy),intent(in):: a,b
real(r_kind), intent(out) :: res_glob

real(r_kind) :: res
integer(i_kind) i,j,imax,jmax,ierr
!-----------------------------------------------------------------------


      res=0.
      res_glob=0.


      if(Fleast(1)) then
        imax=im
      else
        imax=im-1
      endif     
      if(Flnorth(1)) then
        jmax=jm
      else
        jmax=jm-1
      endif     

      do j=0,jmax
      do i=0,imax
        res=res+a(i,j)*b(i,j)
      end do
      end do

!-----------------------------------------------------------------------

      call MPI_REDUCE(res,res_glob,1,dtype,MPI_SUM,0,mpi_comm_comp,ierr)

!-----------------------------------------------------------------------
                        endsubroutine multiply_add                         

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend                               &
!***********************************************************************
!                                                                      *
!         Upsending data from the one resolution pes                   *
!         to the next low-resolution pes                               *
!                                                                      *
!***********************************************************************
(Harray,Warray,km,Lm,mygen_dn,mygen_up)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km,Lm
real(r_kind), dimension(0:imL,0:jmL,Lm,km),intent(in):: Harray
real(r_kind), dimension(-hx:im+hx,-hy:jm+hy,Lm,km),intent(out):: Warray
integer(i_kind),intent(in):: mygen_dn,mygen_up

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:,:)::                          &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_SW
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_SE
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_NW
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
!
! Define generational flags
!
 
     if(mygen_dn==1) then
         g_ind=1
       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)
     else 
         g_ind=2 
       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)
     endif


       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:,:)=0.
   endif

     ndata =km*(imL+1)*(jmL+1)*Lm


      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_SW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                sBuf_SW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      end if

!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(i,j,L,:)=dBuf_SW(i,j,L,:)
             enddo
             enddo
             enddo
          
        else

        allocate( rBuf_SW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(i,j,L,:)=Rbuf_SW(i,j,L,:)
             enddo
             enddo
             enddo

        endif

      endif

      call barrierMPI

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_SE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               sBuf_SE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)

        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(imL+i,j,L,:)=dBuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_SE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(imL+i,j,L,:)=Rbuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

       

        endif

      endif

      call barrierMPI

!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               dBuf_NW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               sBuf_NW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
 
        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(i,jmL+j,L,:)=dBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(i,jmL+j,L,:)=rBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

        end if

      end if

      call barrierMPI
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               dBuf_NE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               sBuf_NE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(4), isend)

         call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(imL+i,jmL+j,L,:)=dBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Warray(imL+i,jmL+j,L,:)=rBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

          deallocate( rBuf_NE, stat = iderr)

        endif
      endif

      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend                             &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                       (MPI version)                                  *
!                                                                      *
!***********************************************************************
(Warray,Harray,km,Lm,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne     &
                    ,Fitarg_up                                           &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw       
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in)::km,Lm
real(r_kind), dimension(0:im,0:jm,Lm,1:km),intent(in):: Warray
real(r_kind), dimension(-1:imL+1,-1:jmL+1,Lm,1:km),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                          &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_SW
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_SE
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_NW
real(r_kind),dimension(0:imL,0:jmL,1:Lm,1:km):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
!-----------------------------------------------------------------------
!
! Define generational flags
!

     if(mygen_dn==1) then
         g_ind=1
       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)
     else 
         g_ind=2
       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)
     endif

       itarg_up=Fitarg_up(g_ind)
!TEST
!if(mygen_dn==1) then
!  write(100+mype,'(a,i3,l3,i3)')'mype,lsendup_sw,itarg_up=',mype,lsendup_sw,itarg_up
!  write(100+mype,'(a,i3,l3,i3)')'mype,lsendup_se,itarg_up=',mype,lsendup_se,itarg_up
!  write(100+mype,'(a,i3,l3,i3)')'mype,lsendup_nw,itarg_up=',mype,lsendup_nw,itarg_up
!  write(100+mype,'(a,i3,l3,i3)')'mype,lsendup_ne,itarg_up=',mype,lsendup_ne,itarg_up
!  write(100+mype,'(a      )')'  '
!  call finishMPI
!endif
!TEST

!
      ndata =km*(imL+1)*(jmL+1)*Lm

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation

if(l_hgen) THEN
 
  if(my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw

!TEST
!if(mygen_dn==1) then
!   write(100+mype,'(a,3i5)')'mype,mype_hgen,itargdn_sw=',mype,mype_hgen,itargdn_sw
!endif
!TEST

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_SW(i,j,L,:) = Warray(i,j,L,:)
             enddo
             enddo
             enddo
!TEST
!if(mygen_dn==1) then
!   write(200+mype,'(a,3i5)')'After give: mype,mype_hgen,itargdn_sw=',mype,mype_hgen,itargdn_sw
!endif
!TEST

        else

        allocate( sBuf_SW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                sBuf_SW(i,j,L,:) = Warray(i,j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )
!TEST
!if(mygen_dn==1) then
!   write(200+mype,'(a,3i5)')'After send: mype,mype_hgen,itargdn_sw=',mype,mype_hgen,itargdn_sw
!endif
!TEST

        endif

  endif

ENDIF
!
! --- Receive SW portion of data at lower generation


  if (lsendup_sw ) then
!!  if(my_hgen.ne.mygen_up .and. lsendup_sw  ) then

        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=dBuf_SW(i,j,L,:)
             enddo
             enddo
             enddo
!TEST
!if(mygen_dn==1) then
!   write(300+mype,'(a,3i5)')'After take: mype,mype_hgen,itargup=',mype,mype_hgen,itarg_up
!endif
!TEST

        else


        allocate( rBuf_SW(0:imL,0:jmL,Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=rBuf_SW(i,j,L,:)  
             enddo
             enddo
             enddo
!TEST
!if(mygen_dn==1) then
!   write(300+mype,'(a,3i5)')'After recv: mype,mype_hgen,itargup=',mype,mype_hgen,itarg_up
!endif
!TEST

        deallocate( rBuf_SW, stat = iderr)

        endif

      endif

!TEST
!if(mygen_dn==1) then
!   call finishMPI
!endif
!TEST
      call barrierMPI

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

if(l_hgen) THEN

  if(my_hgen==mygen_up .and.  itargdn_se >= 0 ) then
        nebpe = itargdn_se


        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_SE(i,j,L,:) = Warray(imL+i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               sBuf_SE(i,j,L,:) = Warray(imL+i,j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

        endif

  endif

ENDIF
!
! --- Receive SE portion of data at lower generation

     if( lsendup_se ) then
 
!!  if(my_hgen.ne.mygen_up .and. lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=dBuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_SE(0:imL,0:jmL,Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=Rbuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
       endif

     end if

     call barrierMPI

! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

if(l_hgen) THEN

  if(my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw


        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_NW(i,j,L,:) = Warray(i,jmL+j,L,:)
             enddo
             enddo
             enddo

        else


        allocate( sBuf_NW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                sBuf_NW(i,j,L,:) = Warray(i,jmL+j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )

        endif

  endif

ENDIF
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then
!!  if(my_hgen.ne.mygen_up .and. lsendup_nw ) then

        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=dBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NW(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=Rbuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

        endif

      end if

      call barrierMPI

! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

if(l_hgen) THEN

  if(my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne
        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                dBuf_NE(i,j,L,:) = Warray(imL+i,jmL+j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
                sBuf_NE(i,j,L,:) = Warray(imL+i,jmL+j,L,:)
             enddo
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

        endif

  endif

ENDIF
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
!!  if(my_hgen.ne.mygen_up .and. lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=dBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NE(0:imL,0:jmL,1:Lm,1:km), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do L=1,Lm
             do j=0,jmL
             do i=0,imL
               Harray(i,j,L,:)=rBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

        endif

      end if
   
      call barrierMPI
     

!-----------------------------------------------------------------------
                        endsubroutine downsend

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco05                               &
!***********************************************************************
!                                                                      *
! Half values at the edges to accomodate adjoint of beta filter        *
!                                                                      *
!***********************************************************************
(V,H,km,im,jm,Lm,nbx,nby)
!-----------------------------------------------------------------------
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby
real(r_kind), dimension(-nbx:im+nbx,-nby:jm+nby,Lm,1:km),intent(inout):: V,H
!-----------------------------------------------------------------------
logical lgen
integer(i_kind) L
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!
!
! Define new boundaries
!
     
           do L=1,LM
              V(0 ,0:jm,L,:)=V(0 ,0:jm,L,:)*0.5
              V(im,0:jm,L,:)=V(im,0:jm,L,:)*0.5
              V(0:im,0 ,L,:)=V(0:im,0 ,L,:)*0.5
              V(0:im,jm,L,:)=V(0:im,jm,L,:)*0.5
            end do

        if(l_hgen) then
               
           do L=1,LM
              H(0 ,0:jm,L,:)=H(0 ,0:jm,L,:)*0.5
              H(im,0:jm,L,:)=H(im,0:jm,L,:)*0.5
              H(0:im,0 ,L,:)=H(0:im,0 ,L,:)*0.5
              H(0:im,jm,L,:)=H(0:im,jm,L,:)*0.5
            end do

        endif
!-----------------------------------------------------------------------
                        endsubroutine boco05

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine v02v                                 &
!**********************************************************************!
!                                                                      !
!  Provide common edges for analysis subdomains                        !
!                                                                      !
!  Conversion from 'original' analysis array (which does not share     !
!  neighbors) to 'work' analysis array (which share edges with         !
!  neighbors                                                           !
!                                                                      !
!          o o o o o o          x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o x x x x o    ->    x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o o o o o o          x x x x x x                            !
!                                                                      !
!**********************************************************************!
(WA,nmax,mmax,lmax)
!-----------------------------------------------------------------------
use mg_domain, only: itarg_wA,itarg_eA,itarg_sA,itarg_nA                &
                    ,lwestA,leastA,lsouthA,lnorthA
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: nmax,mmax,lmax
real(r_kind), dimension(0:nmax,0:mmax,lmax),intent(inout):: WA
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:):: sBuf_N,rBuf_S
real(r_kind), allocatable, dimension(:,:):: sBuf_E,rBuf_W
integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,n,m
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical lgen
!-----------------------------------------------------------------------

  
!
! Define boundary conditions
!


!
!                           SEND boundaries toward North
!

      ndatay = nmax*lmax



      if( itarg_nA >= 0 ) then
        nebpe = itarg_nA

           allocate( sBuf_N(1:nmax,1:lmax), stat = iaerr )
           
              do n=1,nmax
                sBuf_N(n,:) = WA(n,mmax,:)
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(1), isend)
      end if

!
!                           RECEIVE boundaries from South
! 

      if( itarg_sA >= 0 ) then
        nebpe = itarg_sA

          allocate( rBuf_S(1:nmax,1:lmax), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

!
!                           ASSIGN boundaries received from South
! 

      if(  itarg_sA >= 0 ) then

          do n=1,nmax
            WA(n,0,:)=rBuf_S(n,:)
          enddo

      endif

                                  call barrierMPI

!----------------------------------------------------------------------

!
!                           SEND boundaries toward Eeast
!

      ndatax = (mmax+1)*lmax

      if( itarg_eA >= 0) then
        nebpe = itarg_eA

              allocate( sBuf_E(0:mmax,1:lmax), stat = iaerr )

              do m=0,mmax
                sBuf_E(m,1:lmax) = WA(nmax,m,1:lmax)
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                             mpi_comm_comp, sHandle(2), isend)

      end if

!
!                           RECEIVE boundaries from Weast
! 

      if( itarg_wA >= 0 ) then
        nebpe = itarg_wA

          allocate( rBuf_W(0:mmax,1:lmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if
!
!                           ASSIGN boundaries from east

      if( itarg_wA >= 0 ) then

          do m=0,mmax
            WA(0,m,1:lmax)=rBuf_W(m,1:lmax)
          enddo
      else 
          do m=0,mmax                            !  Do not need that
            WA(0,m,1:lmax)=0.                    !
          enddo                                  !

      end if


                                  call barrierMPI

!-----------------------------------------------------------------------
!
!                           DEALLOCATE Bufferes
!
      if(itarg_eA >= 0) deallocate( sBuf_E)
      if(itarg_wA >= 0) deallocate( rBuf_W)

      if(itarg_nA >= 0) deallocate( sBuf_N)
      if(itarg_sA >= 0) deallocate( rBuf_S)


!-----------------------------------------------------------------------
                        endsubroutine v02v

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocosH1                              &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies an (im,jm,Lm) array at generation one with halo made of     !
! the first and the last line.  Technically, here (im,jm)=(lon2,lat2)  !
!                                                                      !
!**********************************************************************!
(Warray,im,jm,Lm)
!-----------------------------------------------------------------------
use mpi
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_ne,Fitarg_se,Fitarg_sw,Fitarg_nw            &
                    ,Flcorner_sw,Flcorner_nw,Flcorner_se,Flcorner_ne   
implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: im,jm,Lm
real(r_kind),dimension(im,jm,Lm),intent(inout):: Warray
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:):: sBuf_N,sBuf_E,sBuf_S,sBuf_W &
                                           ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           
real(r_kind), dimension(LM):: sBuf_NE,sBuf_SE,sBuf_SW,sBuf_NW           &
                             ,rBuf_NE,rBuf_SE,rBuf_SW,rBuf_NW        

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e                         &
               ,itarg_ne,itarg_se,itarg_sw,itarg_nw                     &
               ,imax,jmax
logical:: lcorner_sw,lcorner_nw,lcorner_se,lcorner_ne                   &
         ,lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay,nbxy
integer(i_kind) g_ind
!-----------------------------------------------------------------------

      g_ind = 1

!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          itarg_ne = Fitarg_ne(g_ind)
          itarg_se = Fitarg_se(g_ind)
          itarg_sw = Fitarg_sw(g_ind)
          itarg_nw = Fitarg_nw(g_ind)

          lcorner_sw = Flcorner_sw(g_ind)
          lcorner_nw = Flcorner_nw(g_ind)
          lcorner_se = Flcorner_se(g_ind)
          lcorner_ne = Flcorner_ne(g_ind)   

     imax=im-2
     jmax=jm-2

!-----------------------------------------------------------------------
      ndatay = imax*Lm
      ndatax = jmax*Lm
      nbxy   =      Lm

!
!     SEND boundaries 
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:imax,1:Lm), stat = iaerr )

              do L=1,Lm
                do i=1,imax
                  sBuf_S(i,L) = Warray(i+1,1,L)
                enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:imax,1:Lm), stat = iaerr )

              do L=1,Lm
                do i=1,imax
                  sBuf_N(i,L)=Warray(i+1,jm,L)
                enddo
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:jmax,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,jmax
                  sBuf_W(j,L) = Warray(1,j+1,L)
                enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:jmax,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,jmax
                  sBuf_E(j,L) = Warray(im,j,L)
                enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
!     RECEIVE boundaries 
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:imax,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:imax,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if


! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:jmax,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:jmax,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
!                           DEALLOCATE sBufferes
!


      if( itarg_n >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if

!-----------------------------------------------------------------------
!
!                           SEND corners
!
! --- toward SOUTH-WEST ---

      if( itarg_sw >= 0 ) then
        nebpe = itarg_sw

          do L=1,Lm
            sBuf_SW(L)= Warray(1,1,L)
          enddo


        call MPI_ISEND( sBuf_SW, nbxy, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se

          do L=1,Lm
            sBuf_SE(L)=Warray(im,1,L)
          enddo

        call MPI_ISEND( sBuf_SE, nbxy, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(2), isend)
      end if

! --- toward NORTH-EAST ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne

          do L=1,Lm
            sBuf_NE(L) = Warray(im,jm,L)
          enddo

        call MPI_ISEND( sBuf_NE, nbxy, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(1), isend)
      end if


! --- toward NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw

          do L=1,Lm
            sBuf_NW(L) = Warray(1,jm,L)
          enddo

       call MPI_ISEND( sBuf_NW, nbxy, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(4), isend)
      end if

!
!                           RECEIVE corners
!
! --- from NORTH-EAST  ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne
        call MPI_IRECV( rBuf_NE, nbxy, dtype, nebpe, nebpe,  &
                        mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )
      end if

! --- from NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw
        call MPI_IRECV( rBuf_NW, nbxy, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )
      end if

! --- from SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se
        call MPI_IRECV( rBuf_SE, nbxy, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )
      end if

! --- from SOUTH-WEST ---

      if(  itarg_sw >= 0 ) then
        nebpe = itarg_sw
        call MPI_IRECV( rBuf_SW, nbxy, dtype, nebpe, nebpe,  &
                        mpi_comm_comp, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )
      end if


!
! Assign received values from NORTH, SOUTH, EAST and WEST
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do i=1,imax
       Warray(1+i,jm,L)=Warray(1+i,jm-1,L)
     enddo
     enddo

   else

     do L=1,Lm
     do i=1,imax
       Warray(1+i,jm,L)=rBuf_N(i,L)
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do i=1,imax
       Warray(1+i,1,L)=Warray(1+i,2,L)
     end do
     end do

   else

     do L=1,Lm
     do i=1,imax
       Warray(1+i,1,L)=rBuf_S(i,L)
     enddo
     enddo

   endif

! From west

   if(lwest) then

     do L=1,Lm
     do j=1,jmax
       Warray(1,1+j,L)= Warray(2,1+j,L)
     end do
     end do

   else 

     do L=1,Lm
     do j=1,jmax
       Warray(1,1+j,L)= rBuf_W(i,L)
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=1,jmax
       Warray(im,1+j,L)=Warray(im-1,1+j,L)
     end do
     end do

   else 

     do L=1,Lm
     do j=1,jmax
       Warray(im,1+j,L)=rBuf_E(j,L)
     enddo
     enddo

   endif

! From South-West

   if(lcorner_sw)then

     do L=1,Lm
       Warray(1,1,L) = Warray(2,2,L)
     end do

   else &
   if(lwest) then

     do L=1,Lm
       Warray(1,1,L) = Warray(2,1,L)
     end do

   else &
   if(lsouth) then

     do L=1,Lm
       Warray(1,1,L) = Warray(1,2,L)
     end do

   else 

     do L=1,Lm
       Warray(1,1,L) = rBuf_SW(L)
     end do

   end if


! From North-West

   if(lcorner_nw)then

     do L=1,Lm
       Warray(1,jm,L) = Warray(2,jm-1,L)
     end do
    
   else &
   if(lwest) then

     do L=1,Lm
        Warray(1,jm,L) = Warray(2,jm,L)
      end do

   else &
   if(lnorth) then

     do L=1,Lm
       Warray(1,jm,L) = Warray(1,jm-1,L)
     enddo

   else 

     do L=1,Lm
       Warray(1,jm,L)=rBuf_NW(L)
     end do
      
   end if


! From South-East

   if(lcorner_se)then

     do L=1,Lm
       Warray(im,1,L) = Warray(im-1,2,L)
     end do

   else &
   if(least) then

     do L=1,Lm
       Warray(im,1,L) = Warray(im-1,1,L)
     end do

   else &
   if(lsouth) then

     do L=1,Lm
       Warray(im,1,L) = Warray(im,2,L)
     end do

   else

     do L=1,Lm
       Warray(im,1,L) = rBuf_SE(L)
     end do

   end if

! From North-East


   if(lcorner_ne)then
    
     do L=1,Lm
       Warray(im,jm,L) = Warray(im-1,jm-1,L)
     end do

   else &
   if(least) then

     do L=1,Lm
       Warray(im,jm,L) = Warray(im-1,jm,L)
     end do

   else &
   if(lnorth) then

     do L=1,Lm
       Warray(im,jm,L) = Warray(im,jm-1,L)
     enddo

   else

     do L=1,Lm
       Warray(im,jm,L) = rBuf_NE(L)
     end do

   end if


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if





!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocosH1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_2d_g1                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions just for generation 1                                     !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby,Fimax,Fjmax) 
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                     

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: &
                                  Warray
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay,nbxy
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby
      ndatax = km*(jmax+1+2*nby)*nbx


!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j)=Warray(:,i,jmax-nby-1+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=Warray(:,i,jmax-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=Warray(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif


!----------------------------------------------------------------------
!
! SEND extended boundaries toward WEST and EAST 
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,-nby:jmax+nby), stat = iaerr )

                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,-nby:jmax+nby), stat = iaerr )

                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = Warray(:,imax-nbx-1+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST 
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if


!
! Assign received values from EAST and WEST
!

! From west

   if(lwest) then

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= Warray(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,imax+i,j)=Warray(:,imax-i,j)
     end do
     end do

   else 

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_2d_g1                          &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions just for generation 1                                     !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatax =km*(jmax+1+2*nby)*(nbx+1)
      ndatay =km*(imax+1)*(nby+1)

!
! SEND extended halos towiard WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )

              do j=-nby,jmax+nby
              do i=-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )

              do j=-nby,jmax+nby
              do i=0,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!

! From west

   if(lwest) then
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,i,j)= W(:,i,j)+W(:,-i,j)
     end do
     end do
   else
     do j=-nby,jmax+nby
     do i=0,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx-i,j)
     end do
     end do
   else 
     do j=-nby,jmax+nby
     do i=0,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND boundaries SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby), stat = iaerr )

              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_world, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby), stat = iaerr )

              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_world, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! ASSIGN received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+W(:,i,-j)
     end do
     end do
   else
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby-j)
     enddo
     enddo
   else
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!-----------------------------------------------------------------------
                        endsubroutine bocoT_2d_g1


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_2d_gh                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of domain by assuming mirror boundary conditions !
! for generations higher then one                                      !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: &
                                  Warray
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby
      ndatax = km*(jmax+1+2*nby)*nbx


!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j)=Warray(:,i,jmax-nby-1+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=Warray(:,i,jmax-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=Warray(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif

!
!  SEND extended boundaries to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,-nby:jmax+nby), stat = iaerr )

                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,-nby:jmax+nby), stat = iaerr )

                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = Warray(:,imax-nbx-1+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= Warray(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,imax+i,j)=Warray(:,imax-i,j)
     end do
     end do

   else 

     do j=-nby,jmax+nby
     do i=1,nbx
       Warray(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if




!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_2d_gh                          &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain. Developed for high grid generations.      *
!                                                                      *
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+1+2*nby)*(nbx+1)
      ndatay =km*(imax+1)*(nby+1)
!
! SEND extended halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )

              do j=-nby,jmax+nby
              do i=-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)
!T                              mpi_comm_world, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )

              do j=-nby,jmax+nby
              do i=0,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)
!T                              mpi_comm_world, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
!T                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,0:nbx,-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
!T                       mpi_comm_world, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,i,j)= W(:,i,j)+W(:,-i,j)
     end do
     end do
   else
     do j=-nby,jmax+nby
     do i=0,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx-i,j)
     end do
     end do
   else 
     do j=-nby,jmax+nby
     do i=0,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby), stat = iaerr )

              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
!T                              mpi_comm_world, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby), stat = iaerr )

              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)
!T                             mpi_comm_world, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
!T                      mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
!T                       mpi_comm_world, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+W(:,i,-j)
     end do
     end do
   else
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby-j)
     enddo
     enddo
   else
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoT_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_g1_old                    &
!***********************************************************************
!                                                                      *
!         Upsend data from generation one to generation two            *
!                                                                      *
!***********************************************************************
(Harray,Warray,Lm_all)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: Lm_all
real(r_kind), dimension(lm_all,0:imL,0:jmL),intent(in):: Harray
real(r_kind), dimension(lm_all,-hx:im+hx,-hy:jm+hy),intent(out):: Warray

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SE
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!
       g_ind=1

       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)


       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:)=0.
   endif

     ndata =lm_all*(imL+1)*(jmL+1)


      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=0,jmL
             do i=0,imL
                dBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      end if

!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
          
        else

        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,i,j)=Rbuf_SW(:,i,j)
             enddo
             enddo

        endif

      endif

      call barrierMPI

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
                dBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)

        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

        else

        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       

        endif

      endif

      call barrierMPI

!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
 
        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

        else

        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,i,jmL+j)=rBuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

        end if

      end if

      call barrierMPI
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(4), isend)

         call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

        else

        allocate( rBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,jmL+j)=rBuf_NE(:,i,j)
             enddo
             enddo

          deallocate( rBuf_NE, stat = iderr)

        endif
      endif

      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_g1_old

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_gh                        &
!***********************************************************************
!                                                                      *
!         Upsend data from one grid generation to another              *
!         (Just for high grid generations)                             *
!                                                                      *
!***********************************************************************
(Harray,Warray,Lm_all,mygen_dn,mygen_up)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: Lm_all
real(r_kind), dimension(lm_all,0:imL,0:jmL),intent(in):: Harray
real(r_kind), dimension(lm_all,-hx:im+hx,-hy:jm+hy),intent(out):: Warray
integer(i_kind),intent(in):: mygen_dn,mygen_up

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SE
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
!
! Define generational flags
!
 
       g_ind=2 

       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:)=0.
   endif

     ndata =lm_all*(imL+1)*(jmL+1)


      if(  lsendup_sw ) then

        nebpe = itarg_up
    

        allocate( sBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
!T                       mpi_comm_comp, sHandle(1), isend)
                       mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )


      end if

!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
!T                       mpi_comm_comp, rHandle(1), irecv)
                       mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,i,j)=Rbuf_SW(:,i,j)
             enddo
             enddo

      endif

!TT      call barrierMPI

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( sBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
!T                       mpi_comm_comp, sHandle(2), isend)
                       mpi_comm_work, sHandle(2), isend)

        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

      end if

!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se


        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
!T                       mpi_comm_comp, rHandle(2), irecv)
                       mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

      endif

!TT      call barrierMPI

!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        allocate( sBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
!T                        mpi_comm_comp, sHandle(3), isend)
                        mpi_comm_work, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )


    end if

!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
 

        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
!T                       mpi_comm_comp, rHandle(3), irecv)
                       mpi_comm_work, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,i,jmL+j)=rBuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

      end if

!TT      call barrierMPI
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( sBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
!T                       mpi_comm_comp, sHandle(4), isend)
                       mpi_comm_work, sHandle(4), isend)

         call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne

        allocate( rBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
!T                       mpi_comm_comp, rHandle(4), irecv)
                       mpi_comm_work, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,jmL+j)=rBuf_NE(:,i,j)
             enddo
             enddo

          deallocate( rBuf_NE, stat = iderr)

      endif

!TT      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_all_g2                      &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                       (MPI version)                                  *
!                                                                      *
!***********************************************************************
(Warray,Harray,lm_all)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne     &
                    ,Fitarg_up                                           &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw       
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: lm_all
real(r_kind), dimension(lm_all,0:im,0:jm),intent(in):: Warray
real(r_kind), dimension(lm_all,-1:imL+1,-1:jmL+1),intent(out):: Harray
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE             

real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SE
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer:: mygen_up,mygen_dn
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
!-----------------------------------------------------------------------
!
! Define generational flags
!
    mygen_up=2
    mygen_dn=1

         g_ind=1
       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)

       itarg_up=Fitarg_up(g_ind)


      ndata =lm_all*(imL+1)*(jmL+1)


!
! Send data down to generation 1
!
LSEND:  if(my_hgen==mygen_up) then
!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation
 
        nebpe = itargdn_sw

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )

        endif
!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

        nebpe = itargdn_se

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

        endif

! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

        nebpe = itargdn_nw

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
                dBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )

        endif

!
! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

        nebpe = itargdn_ne
        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
                dBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

        endif


    endif LSEND   

!
! --- Receive SW portion of data at lower generation
!

      if( lsendup_sw .and. mype /= itarg_up ) then

        nebpe = itarg_up


        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )


      else &

!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se .and. mype /= itarg_up) then

        nebpe = itarg_up

        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )


      else &


!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw .and. mype /= itarg_up) then

        nebpe = itarg_up

        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )


      else &


!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne .and. mype /= itarg_up) then
        nebpe = itarg_up

        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )


      end if
   
!
! Assign received and prescribed values
!     
      if( lsendup_sw ) then

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_se ) then

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_nw ) then

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_ne ) then

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo

       endif


!-----------------------------------------------------------------------
                        endsubroutine downsend_all_g2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_g1                        &
!***********************************************************************
!                                                                      *
!         Upsend data from generation one to generation two            *
!                                                                      *
!***********************************************************************
(Harray,Warray,Lm_all)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: Lm_all
real(r_kind), dimension(lm_all,0:imL,0:jmL),intent(in):: Harray
real(r_kind), dimension(lm_all,-hx:im+hx,-hy:jm+hy),intent(out):: Warray

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SE
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!
       g_ind=1

       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)


       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:)=0.
   endif

!
! --- Send data to SW portion of processors at higher generation
!
     ndata =lm_all*(imL+1)*(jmL+1)


      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=0,jmL
             do i=0,imL
                dBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      end if
!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
                dBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)

        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if
!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=0,jmL
             do i=0,imL
               dBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

!T        if(nebpe == mype) then
        if(nebpe /= mype) then

!T             do j=0,jmL
!T             do i=0,imL
!T               Warray(:,i,j)=dBuf_SW(:,i,j)
!T             enddo
!T             enddo
          
!T        else

!T        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

!T        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )
 
!TEST
          endif
!TEST


             do j=0,jmL
             do i=0,imL
!T               Warray(:,i,j)=rBuf_SW(:,i,j)
               Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

!T        deallocate( rBuf_SW, stat = iderr)

!T        endif

      endif

!      call barrierMPI


!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then

        nebpe = itargdn_se

!T        if(nebpe == mype) then
        if(nebpe /= mype) then

!T             do j=0,jmL
!T             do i=0,imL
!T               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
!T             enddo
!T             enddo

!T        else

!T        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

!T        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

!TEST
        endif
!TEST
             do j=0,jmL
             do i=0,imL
!T               Warray(:,imL+i,j)=rBuf_SE(:,i,j)
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

!T        deallocate( rBuf_SE, stat = iderr)

!T        endif

      endif

!      call barrierMPI


!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
!T        if(nebpe == mype) then
        if(nebpe /= mype) then

!T             do j=0,jmL
!T             do i=0,imL
!T               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
!T             enddo
!T             enddo

!T        else

!T        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

!T        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

!TEST
        endif
!TEST

             do j=0,jmL
             do i=0,imL
!T               Warray(:,i,jmL+j)=rBuf_NW(:,i,j)
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

!T        deallocate( rBuf_NW, stat = iderr)

!T        end if

      endif

!      call barrierMPI

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

!T        if(nebpe == mype) then
        if(nebpe /= mype) then

!T             do j=0,jmL
!T             do i=0,imL
!T               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
!T             enddo
!T             enddo

!T        else

!T        allocate( rBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

!T        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )
!TEST
        endif
!TEST

             do j=0,jmL
             do i=0,imL
!T               Warray(:,imL+i,jmL+j)=rBuf_NE(:,i,j)
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

!T          deallocate( rBuf_NE, stat = iderr)

!T        endif

      endif

!      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_all_gh                      &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                       (MPI version)                                  *
!                                                                      *
!***********************************************************************
(Warray,Harray,lm_all,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne     &
                    ,Fitarg_up                                           &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw       
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: lm_all
real(r_kind), dimension(lm_all,0:im,0:jm),intent(in):: Warray
real(r_kind), dimension(lm_all,-1:imL+1,-1:jmL+1),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_SE
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NW
real(r_kind),dimension(1:lm_all,0:imL,0:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
!-----------------------------------------------------------------------
     Harray = 0.
!
! Define generational flags
!

         g_ind=2
       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)

       ndata =lm_all*(imL+1)*(jmL+1)

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation

 
  if(my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw


        allocate( sBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
!T                        mpi_comm_comp, sHandle(1), isend)
                        mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )


  endif
!
! --- Receive SW portion of data at lower generation


      if( lsendup_sw ) then

        nebpe = itarg_up


        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
!T                        mpi_comm_comp, rHandle(1), irecv)
                        mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

      endif

!TT      call barrierMPI

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(my_hgen==mygen_up .and.  itargdn_se >= 0 ) then
        nebpe = itargdn_se

        allocate( sBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
               sBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
!T                       mpi_comm_comp, sHandle(2), isend)
                       mpi_comm_work, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

!        endif

  endif
!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
!T                        mpi_comm_comp, rHandle(2), irecv)
                        mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
     end if

!TT     call barrierMPI

!
! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

  if(my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw


        allocate( sBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
!T                        mpi_comm_comp, sHandle(3), isend)
                        mpi_comm_work, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )


  endif
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
!T                       mpi_comm_comp, rHandle(3), irecv)
                       mpi_comm_work, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)


      end if

!TT      call barrierMPI

! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne


        allocate( sBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

             do j=0,jmL
             do i=0,imL
                sBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
!T                        mpi_comm_comp, sHandle(4), isend)
                        mpi_comm_work, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )


  endif
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( rBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
!T                        mpi_comm_comp, rHandle(4), irecv)
                        mpi_comm_work, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=0,jmL
             do i=0,imL
               Harray(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)


      end if
   
!T      call barrierMPI
     

!-----------------------------------------------------------------------
                        endsubroutine downsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_3d_g1                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions                                                           !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby,nbz
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby,1-nbz:Lm+nbz)         &
                      ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay,nbxy
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
!         l_sidesend=.false.
!
!       if(mygen_min==1.and.mygen_max==1) then
!        g_ind=1
!        g = 1
!         l_sidesend=.true.
!       else &
!       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
!         g_ind=2 
!         g = my_hgen
!         l_sidesend=.true.
!       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

! FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind=1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im  
          jmax = jm

!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby*Lm
      ndatax = km*(jmax+1+2*nby)*nbx*Lm
      nbxy   = km*nbx*nby*Lm


!
! SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j,L)=W(:,i,jmax-nby-1+j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,jmax+j,L)=W(:,i,jmax-j,L)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,jmax+j,L)=rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,-nby-1+j,L)=W(:,i,nby+1-j,L)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,-nby-1+j,L)=rBuf_S(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!
! SEND extended boundaries toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j,L) = W(:,imax-nbx-1+i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries WEST and EAST
!

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

!
! Assign received values from  EAST and WEST
!
! From west

   if(lwest) then

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx-1+i,j,L)= W(:,nbx+1-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx-1+i,j,L)= rBuf_W(:,i,j,L)
     enddo
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=W(:,imax-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=rBuf_E(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_3d_gh                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions                                                           !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby,1-nbz:Lm+nbz)         &
                      ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby*Lm
      ndatax = km*(jmax+1+2*nby)*nbx*Lm

!
! SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j,L)=W(:,i,jmax-nby-1+j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from SOUTH and NORTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!TEST
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
!TEST

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,jmax+j,L)=W(:,i,jmax-j,L)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,jmax+j,L)=rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,-nby-1+j,L)=W(:,i,nby+1-j,L)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       W(:,i,-nby-1+j,L)=rBuf_S(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!TEST
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
!TEST

!
!                           DEALLOCATE sBufferes
!

!      if( itarg_n >= 0 ) then
!        call MPI_WAIT( sHandle(1), istat, ierr )
!        deallocate( sBuf_N, stat = ierr )
!      end if
!      if( itarg_s >= 0 ) then
!        call MPI_WAIT( sHandle(3), istat, ierr )
!        deallocate( sBuf_S, stat = ierr )
!      end if
!      if( itarg_e >= 0 ) then
!        call MPI_WAIT( sHandle(2), istat, ierr )
!        deallocate( sBuf_E, stat = ierr )
!      end if
!      if( itarg_w >= 0 ) then

!
! SEND extended boundaries to WEST and EAST   
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j,L) = W(:,imax-nbx-1+i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

!
! Deallocate send bufferes from EAST and WEST
!
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if

!
! Assign received values from WEST and EAST
!
! From west

   if(lwest) then

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx-1+i,j,L)= W(:,nbx+1-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx-1+i,j,L)= rBuf_W(:,i,j,L)
     enddo
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=W(:,imax-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=rBuf_E(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1-L )=W(:,:,:, 1+L)
          W(:,:,:,LM+L)=W(:,:,:,LM-L)
        end do


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_3d_g1                          &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain                                            *
!                                                                      *
!***********************************************************************
(W,km,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby,nbz
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby,1-nbz:Lm+nbz)         &
                       ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

         g_ind=1

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          imax = im
          jmax = jm

!----------------------------------------------------------------------
      ndatax =km*(jmax+1+2*nby)*(nbx+1) *Lm
      ndatay =km*(imax+1)*(nby+1) *Lm

!
! SEND extended halos toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,jmax+nby
              do i=-nbx,0
                sBuf_W(:,i+nbx,j,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,jmax+nby
              do i=0,nbx
                sBuf_E(:,i,j,L) = W(:,imax+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(2), isend)

      end if
!
! RECEIVE extended halos from EAST and WEST
!
! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


          allocate( rBuf_W(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if
!
! Assign received extended halos from WEST and EAST to interior of domains
!

! From west

   if(lwest) then
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,i,j,L)= W(:,i,j,L)+W(:,-i,j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
      W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
     end do
     end do
     end do
   endif

! From east

   if(least) then
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+nbx-i,j,L)
     end do
     end do
     end do
   else 
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx  
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
     end do
     end do
     end do
   endif

!
! Send halos SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_world, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j,L)=W(:,i,jmax+j,L)
              enddo
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_world, sHandle(1), isend)

      end if


!
! RECEIVE boundaries from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,j,L)= W(:,i,j,L)+W(:,i,-j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
     end do
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+nby-j,L)
     enddo
     enddo
     enddo
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo
   endif

!----------------------------------------------------------------------
!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1+L )=W(:,:,:, 1+L)+W(:,:,:, 1-L)
          W(:,:,:,LM-L)=W(:,:,:,LM-L)+W(:,:,:,LM+L)
        end do


!----------------------------------------------------------------------
!
!                           DEALLOCATE sBufferes
!


      if( itarg_w >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif 
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif 
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif 
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif


!-----------------------------------------------------------------------
                        endsubroutine bocoT_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_3d_gh                          &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain                                            *
!                                                                      *
!***********************************************************************
(W,km,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax,mygen_min,mygen_max)

!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi
implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby,1-nbz:Lm+nbz)         &
                       ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+1+2*nby)*(nbx+1) *Lm
      ndatay =km*(imax+1)*(nby+1) *Lm

!
! SEND extended halos toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,jmax+nby
              do i=-nbx,0
                sBuf_W(:,i+nbx,j,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,jmax+nby
              do i=0,nbx
                sBuf_E(:,i,j,L) = W(:,imax+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)
      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


          allocate( rBuf_W(1:km,0:nbx,-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if

!
! Assign received extended halos from WEST and EAST
!

! From west

   if(lwest) then
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,i,j,L)= W(:,i,j,L)+W(:,-i,j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
      W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
     end do
     end do
     end do
   endif

! From east

   if(least) then
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+nbx-i,j,L)
     end do
     end do
     end do
   else 
     do L=1,lm
     do j=-nby,jmax+nby
     do i=0,nbx  
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
     end do
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j,L)=W(:,i,jmax+j,L)
              enddo
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if


!-----------------------------------------------------------------------
!
! Assign received halos from SOUTH and NORTH
!

   if(lsouth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,j,L)= W(:,i,j,L)+W(:,i,-j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
     end do
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+nby-j,L)
     enddo
     enddo
     enddo
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo
   endif


!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1+L )=W(:,:,:, 1+L)+W(:,:,:, 1-L)
          W(:,:,:,LM-L)=W(:,:,:,LM-L)+W(:,:,:,LM+L)
        end do


!-----------------------------------------------------------------------
!
!                           DEALLOCATE sBufferes
!

      if( itarg_w >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoT_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocox_2d_g1                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for generation 1                              !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_w,Fitarg_e,Flwest,Fleast

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W              &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatax = km*(jmax+1)*nbx

!----------------------------------------------------------------------
!
! SEND extended boundaries toward WEST and EAST 
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,0:jmax), stat = iaerr )

                do j=0,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,0:jmax), stat = iaerr )

                do j=0,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = Warray(:,imax-nbx-1+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST 
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if


!
! Assign received values from EAST and WEST
!

! From west

   if(lwest) then

     do j=0,jmax
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= Warray(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=0,jmax
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=0,jmax
     do i=1,nbx
       Warray(:,imax+i,j)=Warray(:,imax-i,j)
     end do
     end do

   else 

     do j=0,jmax
     do i=1,nbx
       Warray(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocox_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_g1                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for generation 1                              !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flsouth,Flnorth                     

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby


!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j)=Warray(:,i,jmax-nby-1+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=Warray(:,i,jmax-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=Warray(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocox_2d_gh                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for high generations                          !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_w,Fitarg_e,Flwest,Fleast,Flsouth,Flnorth

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatax = km*(jmax+1)*nbx

!
!  SEND halos to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,0:jmax), stat = iaerr )

                do j=0,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,0:jmax), stat = iaerr )

                do j=0,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = Warray(:,imax-nbx-1+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=0,jmax
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= Warray(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=0,jmax
     do i=1,nbx
       Warray(:,-nbx-1+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=0,jmax
     do i=1,nbx
       Warray(:,imax+i,j)=Warray(:,imax-i,j)
     end do
     end do

   else 

     do j=0,jmax
     do i=1,nbx
       Warray(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocox_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_gh                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in u direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for high generations                          !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby

!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j)=Warray(:,i,jmax-nby-1+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=Warray(:,i,jmax-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=Warray(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_g1                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nbx lines close to East-West edges of the subdomains,       !
! including values at the edges from halos, assuming mirror boundary   !
! conditions. Just for generation 1                                    ! 
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatax =km*(jmax+1)*(nbx+1)

!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(1), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=0,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_world, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!

! From west

   if(lwest) then
     do j=0,jmax
     do i=0,nbx
       W(:,i,j)= W(:,i,j)+W(:,-i,j)
     end do
     end do
   else
     do j=0,jmax
     do i=0,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=0,jmax
     do i=0,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx-i,j)
     end do
     end do
   else 
     do j=0,jmax
     do i=0,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTy_2d_g1                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nby lines close to North and South edges of the subdomains, !
! including values at the edges, from halos of North and South         !
! neighbors assuming mirror boundary conditions at the boundaries of   !
! the whole domain. Version for generation 1                           ! 
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flsouth,Flnorth,Fitarg_n,Fitarg_s
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatay =km*(imax+1)*(nby+1)

!
! SEND SOUTH and NORTH halos
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby), stat = iaerr )

              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_world, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby), stat = iaerr )

              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_world, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )


      end if

!
! ASSIGN received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+W(:,i,-j)
     end do
     end do
   else
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby-j)
     enddo
     enddo
   else
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_gh                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nby lines close to West and East edges of the subdomains,   !
! including values at the edges, from halos of North and South         !
! neighbors assuming mirror boundary conditions at the boundaries of   !
! the whole domain. Version for high generations                       ! 
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flnorth,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   
integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+1)*(nbx+1)
!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=0,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=0,jmax
     do i=0,nbx
       W(:,i,j)= W(:,i,j)+W(:,-i,j)
     end do
     end do
   else
     do j=0,jmax
     do i=0,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=0,jmax
     do i=0,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx-i,j)
     end do
     end do
   else 
     do j=0,jmax
     do i=0,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTy_2d_gh                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nby lines close to North and South edges of the subdomains, !
! including values at the edges, from halos of North and South         !
! neighbors assuming mirror boundary conditions at the boundaries of   !
! the whole domain. Version for high generations                       !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fleast,Flsouth,Flnorth                             &
                    ,Fitarg_n,Fitarg_s
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::  sBuf_N,sBuf_S            &
                                              ,rBuf_N,rBuf_S
integer(i_kind) itarg_n,itarg_s,itarg_e,imax,jmax
logical least,lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------

      ndatay =km*(imax+1)*(nby+1)
!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby), stat = iaerr )

              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby), stat = iaerr )

              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+W(:,i,-j)
     end do
     end do
   else
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby-j)
     enddo
     enddo
   else
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_gh

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_bocos
