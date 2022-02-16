!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_gh_group                  &
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
use mg_mppstuff, only: mpi_comm_work
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

!
! --- Send data to SW portion of processors at higher generation
!

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

      call barrierMPI


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

      call barrierMPI


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

      call barrierMPI
!

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

      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_gh_group
