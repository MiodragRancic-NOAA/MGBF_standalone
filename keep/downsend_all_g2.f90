!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_all_g2_new                  &
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

      if( lsendup_sw .and. mype<>itarg_up ) then

        nebpe = itarg_up

!        if(nebpe == mype) then
!
!             do j=0,jmL
!             do i=0,imL
!               Harray(:,i,j)=dBuf_SW(:,i,j)
!             enddo
!             enddo
!
!        else

        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=0,jmL
             do i=0,imL
!T               Harray(:,i,j)=rBuf_SW(:,i,j)  
               dBuf_SW(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)


      else &

!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se .and. mype<>itarg_up) then

        nebpe = itarg_up

!        if(nebpe == mype) then
!
!             do j=0,jmL
!             do i=0,imL
!               Harray(:,i,j)=dBuf_SE(:,i,j)
!             enddo
!             enddo
!
!        else

        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=0,jmL
             do i=0,imL
!               Harray(:,i,j)=Rbuf_SE(:,i,j)
               dBuf_SE(:,i,j)=rBuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
!       endif

     else &


!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw .and. mype<>itarg_up) then

        nebpe = itarg_up

!        if(nebpe == mype) then
!
!             do j=0,jmL
!             do i=0,imL
!               Harray(:,i,j)=dBuf_NW(:,i,j)
!             enddo
!             enddo
!
!        else

        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=0,jmL
             do i=0,imL
!               Harray(:,i,j)=Rbuf_NW(:,i,j)
               dBuf_NW(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

!        endif

      else &


!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne .and. mype<>itarg_up) then
        nebpe = itarg_up

!        if(nebpe == mype) then
!
!             do j=0,jmL
!             do i=0,imL
!               Harray(:,i,j)=dBuf_NE(:,i,j)
!             enddo
!             enddo
!
!        else

        allocate( rBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=0,jmL
             do i=0,imL
!T               Harray(:,i,j)=rBuf_NE(:,i,j)
               dBuf_NE(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

!        endif

      end if
   
!
! Assign received values
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
                        endsubroutine downsend_all_g2_new
