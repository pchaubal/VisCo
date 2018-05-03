Program test
  use readsnap
  use utils
  use fields
  use mpi
  implicit none
  Integer :: downsample_size = 20
  type(snapshotdata) :: snapdata
  Real, allocatable  :: xi(:,:)
  Integer :: ierr,i,rank


  



  snapdata%filename = './snapshots/snapshot_041'
  call readposvel(snapdata)
  call downsample_snap(snapdata,downsample_size)


  call MPI_INIT ( ierr )

  call MPI_comm_rank(MPI_comm_World, rank, ierr)

  call calculate_xi(snapdata,xi)

 
  ! allocate(xi1(downsample_size,3))

  ! if (rank == 0) then
  !   xi1=xi
  !   deallocate(xi)

  !   ! print *, ''
  !   ! do i = 1,downsample_size
  ! 	 !  print *,'xi', xi(i,:)
  !   ! end do
  !   ! print *, 'done'
  ! end if
  
  ! call MPI_BCAST(xi1, Size(xi1), Mpi_Real, 0, MPI_comm_World, IERR)

  print *, ''
  do i = 1,downsample_size
  	  print *,'xi',rank, xi(i,:)
  end do


  call MPI_FINALIZE ( ierr )



End program test