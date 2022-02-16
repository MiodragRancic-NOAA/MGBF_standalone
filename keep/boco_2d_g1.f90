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
