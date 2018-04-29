module readsnap 
implicit none 

contains      

subroutine readposvel(positions,velocities,a)
   integer, parameter :: SNAPSHOT = 0       ! number of dump
   integer, parameter :: FILES = 1          ! number of files per snapshot

   ! character*200, parameter :: path='..'
   character*200 filename
   integer*4 :: flag_sfr,flag_feedback,nstart=0,fn
   integer*4 npart(0:5), nall(0:5)
   real*8    massarr(0:5)
   real*8    a
   real*8    redshift
   integer*4 unused(34)
   integer*4 Ntot, N
   real*4,allocatable    :: PartPos(:,:), positions(:,:)
   real*4,allocatable    :: pos(:,:)
   real*4,allocatable    :: vel(:,:), velocities(:,:)


   filename= '/snapshots/snapshot_041'
   ! filename= '../snapshot_021'
   ! print *,'opening...  '//filename(1:len_trim(filename))
   
   
   open (1, file=filename, form='unformatted')
   read (1) npart, massarr ,a, redshift, flag_sfr,flag_feedback, nall, unused


   N=sum(npart)
   allocate(pos(1:3,1:N))
   allocate(vel(1:3,1:N))

   read (1) pos
   read (1) vel
   close (1)
   
   ! print *,'Redshift=',redshift, 'a=',a
   
   ! do i=1,6
   !    print *,'Type=',i,'    Particles=', nall(i), "mass=", massarr(i)
   ! end do
   
   ! Now we just read in the coordinates, and only keep those of type=1
   ! The array PartPos(1:3,1:Ntot) will contain all these particles
   Ntot= nall(1)
   allocate(PartPos(1:3, 1:Ntot))
   
   ! print *,'reading...  '//filename(1:len_trim(filename))
   
   
   
   
   PartPos(1:3,1+nstart:nstart+npart(1))=pos(1:3, 1 + npart(0):sum(npart(0:1)))
   deallocate(pos)
   
   !Transposing PartPos matrix and storing it in positions
   allocate(positions(1:Ntot,1:3))
   positions = Transpose(PartPos)
   deallocate(PartPos)

   !Trans the velocity matrix
   allocate(velocities(1:Ntot,1:3))
   velocities = Transpose(vel)
   deallocate(vel)

   
   ! print *,'Done with reading.'
 

   ! At this point, the coordinates of all particles of type 1 will
   ! be stored in the array PartPos(,)    
end subroutine readposvel
   
end module readsnap

