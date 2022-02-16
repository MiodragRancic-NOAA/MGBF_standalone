!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_g1_new                    &
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

!   if(my_hgen==mygen_up) then
      Warray(:,:,:)=0.
!   endif

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

        if(nebpe /= mype) then

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,i,j)=dBuf_SW(:,i,j)
!             enddo
!             enddo
!          
!        else

!T        allocate( rBuf_SW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

!        endif
!
!             do j=0,jmL
!             do i=0,imL
!               Warray(:,i,j)=dBuf_SW(:,i,j)
!             enddo
!             enddo

      else &

!      call barrierMPI


!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then

        nebpe = itargdn_se

        if(nebpe /= mype) then

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
!             enddo
!             enddo

!        else

!        allocate( rBuf_SE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
!             enddo
!             enddo

        endif

      else &

!      call barrierMPI


!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
        if(nebpe /= mype) then

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
!             enddo
!             enddo

!        else

!        allocate( rBuf_NW(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
!             enddo
!             enddo

!        deallocate( rBuf_NW, stat = iderr)

        end if

      else &

!      call barrierMPI

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

        if(nebpe /= mype) then

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
!             enddo
!             enddo

!        else

!        allocate( dBuf_NE(1:lm_all,0:imL,0:jmL), stat = iaerr )

        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )

!             do j=0,jmL
!             do i=0,imL
!               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
!             enddo
!             enddo

!          deallocate( rBuf_NE, stat = iderr)

        endif

      endif
!
! Assign received values
!
      if( my_hgen==mygen_up) then

             do j=0,jmL
             do i=0,imL
               Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

             do j=0,jmL
             do i=0,imL
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

             do j=0,jmL
             do i=0,imL
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

      endif

!      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_g1_new
