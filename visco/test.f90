Program test
  use readsnap
  use utils
  use fields
  implicit none
  Integer :: downsample_size = 20
  type(snapshotdata) :: snapdata
  Real, allocatable  :: xi(:,:)
  Integer :: ierr,i


  



  snapdata%filename = './snapshots/snapshot_041'
  call readposvel(snapdata)
  call downsample_snap(snapdata,downsample_size)


  call MPI_INIT ( ierr )

  call calculate_xi(snapdata,xi)

  call MPI_FINALIZE ( ierr )

  print *, ''
  do i = 1,downsample_size
  	print *, xi(i,:)
  end do

End program test