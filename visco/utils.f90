Module utils
implicit none

contains

  subroutine downsample_snap(snapdata,downsample_size)
      use readsnap
      implicit none
      Integer, intent(in) :: downsample_size
      type(snapshotdata)  :: snapdata
      Real*4, allocatable :: xds(:,:), vds(:,:)
  
  
      call downsample(snapdata%positions,downsample_size,xds)
      call downsample(snapdata%velocities,downsample_size,vds)
!       xds = xds/1000.0 ! this is only because the snapdata is in kpc and output data in Mpc
      deallocate(snapdata%positions)
      deallocate(snapdata%velocities)
      allocate(snapdata%positions(downsample_size,3))
      allocate(snapdata%velocities(downsample_size,3))
  
      snapdata%positions  = xds
      snapdata%velocities = vds
      snapdata%mass       = snapdata%mass*(snapdata%n_particles/downsample_size)
      snapdata%n_particles = downsample_size
  end subroutine downsample_snap


  subroutine downsample(array,downsample_size,downsampled)
    	! USE IFPORT
      implicit none
      real :: array(:,:)
      real, allocatable :: downsampled(:,:)
      integer, intent(in) :: downsample_size
      integer :: i,random_index
      ! integer :: j
      integer :: array_size
  
      allocate(downsampled(downsample_size,3))

      array_size = size(array,1)
      
      do i=1,downsample_size
        random_index =  int(rand(0)*(array_size+1-1))+1 
        downsampled(i,:) = array(random_index,:)
      end do
  
  end subroutine downsample
end Module utils