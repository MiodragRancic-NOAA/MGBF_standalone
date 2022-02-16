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

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_g1
