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
      nbxy   =km*(nbx+1)*(nby+1)*Lm

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

