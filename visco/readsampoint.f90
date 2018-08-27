Module readsampoint

contains

subroutine read_spts(smp_pts)
    character*200       ::  filename
    Real*4              :: smp_pts(560000,3)

    filename = './sampling_point.dat'

    open (1, file=filename, status='old')
    read (1,*) smp_pts

    print*, smp_pts(100:105,:)

end subroutine read_spts



end module readsampoint