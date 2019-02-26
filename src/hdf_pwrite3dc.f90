! Program to use MPI_Cart and Parallel HDF5
!
program hdf_pwrite3dc

        use mpi
        use hdf5
        use kinds, only : r_dp

        implicit none

        ! Local array size with halo

        integer, parameter :: ndims = 3
        integer, parameter :: halo  = 0

        integer :: g_N   = 8000 ! Image slab size (in pixels)
        integer :: ch_N  = 1000 ! Chunk size (in pixels)
        integer :: argc         ! Number of command line arguments
        integer :: ierr         ! Error status
        integer :: id           ! My rank/ID
        integer :: np           ! Number of processors
        integer :: iunit        ! File descriptor
        integer :: i,j          ! Loop indexers
        integer :: n(ndims)     ! Local N for i and j directions
        integer :: total(ndims) ! Local total dimension size

        ! MPI IO/Lustre file striping
        integer :: lcount       ! Lustre count size
        integer :: lsize        ! Lustre stripe size
        character(len=1024) :: clcount, clsize ! Strings of LFS

        integer :: info                 ! MPI IO Info
        integer :: m_dims(ndims)        ! MPI cart dims
        integer :: coords(ndims)        ! Co-ords of procs in the grid
        logical :: is_periodic(ndims)   ! Periodic boundary conditions
        logical :: reorder              ! Reorder the MPI structure
        integer :: MPI_COMM_3D          ! New communicator

        character(len=1024) :: filename
        integer :: nsets, nsetsproc ! total nb of sets, nb of sets per process
        ! nb of hdf5 datasets, nb of sets per dataset and process
        integer :: ndsets, nsetsdsetproc 
        integer(kind=hid_t) :: p_id, f_id, x_id, d_id, c_id
        integer(kind=hid_t) :: memspace, filespace
        ! Chunk sizes
        integer(kind=hsize_t) :: c_size(ndims)
        ! Local hyper slab info
        integer(kind=hsize_t) :: d_size(ndims), s_size(ndims), h_size(ndims), &
                                 stride(ndims), block(ndims)
        ! Global hyper slab info
        integer(kind=hsize_t) :: g_size(ndims), g_start(ndims)

        ! Local data array
        real(kind=r_dp), allocatable :: ld(:,:)

        ! Micellenaous variables
        character(len=1024) :: hostname ! Host name
        character(len=1024) :: sbuffer ! buffer for cmd arg
        character(len=256)  :: romio_cb_write = "enable"  ! romio_cb_write option
        character(len=256)  :: romio_ds_write = "disable" ! romio_ds_write option
        integer :: dsetfilling = 0 ! Turn hdf5 dataset dsetfilling on (default off)
        integer :: buflen
        ! Clocks array and rate
        integer :: clock_counts(4)
        integer :: clock_rate
        ! Couters
        integer :: iset, idset
        ! (transfer) statistics 
        real :: stat_total, stat_time_create, stat_rate_create, stat_time_wr, stat_rate_wr
        
        argc = 0
        ierr = 0
        m_dims = (/ 0, 0, 0/)
        is_periodic = .false.      ! Non-periodic
        reorder     = .false.      ! Not allowed to reorder

        ! Init MPI
        call mpi_init(ierr)
        call mpi_comm_size(MPI_COMM_WORLD, np, ierr)

        ! Set up the MPI cartesian topology
        !call mpi_dims_create(np, ndims, m_dims, ierr)
        m_dims(1) = 1
        m_dims(2) = 1
        m_dims(3) = np
        call mpi_cart_create(MPI_COMM_WORLD, ndims, m_dims, is_periodic, &
                             reorder, MPI_COMM_3D, ierr)
        call mpi_comm_rank(MPI_COMM_3D, id, ierr)
        call mpi_cart_coords(MPI_COMM_3D, id, ndims, coords, ierr)

        call mpi_get_processor_name(hostname, buflen, ierr)

        ! Get program arguments        
        write(0,*) "Process:", id, ", chckpoint:", 1, ", hostname: ", trim(hostname(1:buflen))

        if (id .eq. 0) then
           argc = command_argument_count()
           ! get the filename
           if (argc .lt. 1 ) then
              write(0, *) 'Must supply a filename'
              call exit(1)
           endif
           call get_command_argument(1, filename)
           
           ! get the total number of sets
           if (argc .lt. 2 ) then
              write(0, *) 'Must supply a total nb of 2d datasets'
              call exit(1)
           endif
           call get_command_argument(2, sbuffer)
           read(sbuffer, *) nsets

           ! get the number of datasets
           if (argc .lt. 3 ) then
              write(0, *) 'Must supply a nb of hdf5 datasets'
              call exit(1)
           endif
           call get_command_argument(3, sbuffer)
           read(sbuffer, *) ndsets
           
           if (mod(nsets,np) .ne. 0) then
              write(0, *) 'Must use divisiable number of procs (', nsets, '/', np, ').'
              call mpi_abort(MPI_COMM_WORLD, 1, ierr)
           endif

           if (mod(nsets,ndsets) .ne. 0) then
              write(0, *) 'Must use divisiable number of hdf5 datasets (', nsets, '/', ndsets, ').'
              call mpi_abort(MPI_COMM_WORLD, 1, ierr)
           endif
           
           if (mod(nsets,ndsets*np) .ne. 0) then
              write(0, *) 'Must use divisiable number of proc and hdf5 datasets (', nsets, '/', ndsets, '/', np, ').'
              call mpi_abort(MPI_COMM_WORLD, 1, ierr)
           endif
           
           nsetsproc = nsets / np
           nsetsdsetproc = nsetsproc / ndsets

           ! read additional parameters from file
           if (argc .ge. 4 ) then
              call get_command_argument(4, sbuffer)
              open(unit=2, file=sbuffer)
              read(2, *) g_N
              read(2, *) ch_N
              read(2, *) romio_cb_write
              read(2, *) romio_ds_write
              read(2, *) dsetfilling
              close(2)
           endif
           
           print *, "np: ", np, ", nsets: ", nsets, ", ndsets: ", ndsets, ", g_N: ", g_N, ", ch_N: ", ch_N
           print *, "romio_cb_write: ", trim(romio_cb_write), ", romio_ds_write: ", trim(romio_ds_write), ", filling: ", dsetfilling 
        endif

        ! Broadcast the filename
        call mpi_bcast(filename, len(filename), MPI_CHAR, 0, &
                       MPI_COMM_WORLD, ierr)

        ! Broadcast the number of sets
        call mpi_bcast(nsets, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(ndsets, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(nsetsproc, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(nsetsdsetproc, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)


        ! Broadcast the slab size (g_N, ch_N), romio and dataset filling settings
        call mpi_bcast(g_N, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(ch_N, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)
        
        call mpi_bcast(romio_cb_write, len(romio_cb_write), MPI_CHAR, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(romio_ds_write, len(romio_ds_write), MPI_CHAR, 0, &
                       MPI_COMM_WORLD, ierr)

        call mpi_bcast(dsetfilling, 1, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ierr)
        
        ! Init the HDF5 library
        call h5open_f(ierr)

        ! Set a stripe count of 4 and a stripe size of 4MB
        ! lcount = 4
        ! lsize  = 4 * 1024 * 1024
        !lcount = 1
        !lsize = 8 * 1048576
        !write(clcount, '(I4)') lcount
        !write(clsize, '(I8)') lsize

        call mpi_info_create(info, ierr)
        ! Enable the  collective  buffering  optimisation
        !call mpi_info_set(info, "romio_print_hints", "enable", ierr)
        !call mpi_info_set(info, "cb_nodes", "15", ierr)
        !call mpi_info_set(info, "cb_buffer_size", "3944304", ierr)
        !call mpi_info_set(info, "ind_wr_buffer_size", "1972152", ierr)
        call mpi_info_set(info, "romio_cb_write", trim(romio_cb_write), ierr)
        call mpi_info_set(info, "romio_ds_write", trim(romio_ds_write), ierr)
        !call mpi_info_set(info, "cb_nodes", "2", ierr)
        !call mpi_info_set(info, "romio_no_indep_rw", "enable", ierr)
        ! set striping parameters
        !call mpi_info_set(info, "striping_factor", trim(clcount), ierr)
        !call mpi_info_set(info, "striping_unit", trim(clsize), ierr)

        ! Set up the access properties
        call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, ierr)
        call h5pset_fapl_mpio_f(p_id, MPI_COMM_3D, info, ierr)
        call h5pset_fapl_mpio_f(p_id, MPI_COMM_WORLD, info, ierr)

        ! Open the file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, ierr, &
                         access_prp = p_id)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to open: ', trim(filename), ': ', ierr
                call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        endif

        write(0,*) "Process:", id, ", chckpoint:", 2

        ! Generate our local matrix
        do i = 1, ndims-1
                n(i) = g_N / m_dims(i)
                total(i) = n(i) + (2 * halo)
        end do
        n(ndims) = 1
        total(ndims) = n(ndims) + (2 * halo)

        if (halo .ne. 0) then
                allocate(ld(0:total(1)-1, 0:total(2)-1), stat=ierr)
        else
                allocate(ld(total(1),total(2)), stat=ierr)
        end if
        if (ierr .ne. 0) then
                write(0,*) 'Unable to allocate local data array: ', ierr
                call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ld = -99.99
        ! init the local data
        do j = 1, n(2)
                do i = 1, n(1)
                        ld(i,j) = id
                enddo
        enddo

        write(0,*) "Process:", id, ", chckpoint:", 3

        ! Create the local memory space and hyperslab
        do i = 1, ndims
                d_size(i) = total(i)
                s_size(i) = n(i)
                h_size(i) = halo
                stride(i) = 1
                block(i)  = 1
        enddo

        call h5screate_simple_f(ndims, d_size, memspace, ierr)
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                   h_size, s_size, ierr,       &
                                   stride, block)

        write(0,*) "Process:", id, ", chckpoint:", 4

        ! Create the global file space and hyperslab
        g_size = (/ g_N, g_N, nsets/ndsets /)
        do i = 1, ndims
                g_start(i) = n(i) * coords(i)
        enddo

        call h5screate_simple_f(ndims, g_size, filespace, ierr)
        !call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
        !                           g_start, s_size, ierr,       &
        !                           stride, block)

        write(0,*) "Process:", id, ", chckpoint:", 5

        ! Create a data chunking property
        c_size = (/ ch_N, ch_N, 1/)
        call h5pcreate_f(H5P_DATASET_CREATE_F, c_id, ierr)
        write(0,*) "Process:", id, ", chckpoint:", 7
        call h5pset_chunk_f(c_id, ndims, c_size, ierr)
        ! Set dataset fill property
        if (dsetfilling .eq. 1) then
           call h5pset_fill_time_f(c_id, H5D_FILL_TIME_ALLOC_F, ierr)
        elseif (dsetfilling .eq. 2) then
           call h5pset_fill_time_f(c_id, H5D_FILL_TIME_ERROR_F, ierr)
        else
           call h5pset_fill_time_f(c_id, H5D_FILL_TIME_NEVER_F, ierr)
        endif

        write(0,*) "Process:", id, ", chckpoint:", 8

        ! ----------------------------------------------------------------------

        call mpi_barrier(MPI_COMM_3D, ierr)

        if (id .eq. 0) then
           call system_clock(clock_counts(1), clock_rate)
        endif
        
        ! Create the datasets
        do idset = 1, ndsets
           ! Dataset name
           write(sbuffer,"(A6,i0.4)") "/data_", idset-1
           
           call h5dcreate_f(f_id, trim(sbuffer), H5T_IEEE_F64LE, filespace, d_id, &
                            ierr, dcpl_id=c_id)

           call h5dclose_f(d_id, ierr)
        end do
        
        write(0,*) "Process:", id, ", chckpoint:", 9

        ! ----------------------------------------------------------------------

        call mpi_barrier(MPI_COMM_3D, ierr)

        if (id .eq. 0) then
           call system_clock(clock_counts(2), clock_rate)
        endif

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        !call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_INDEPENDENT_F, ierr)

        ! Do write for all datasets associated by process
        do idset = 1, ndsets
           ! Dataset name
           write(sbuffer,"(A6,i0.4)") "/data_", idset-1
           
           ! Reopen the dataset, get dataset id
           call h5dopen_f(f_id, trim(sbuffer), d_id, ierr)

           ! First process hyperslab
           g_start(3) = n(3) * coords(3)
           
           do i = 1, nsetsdsetproc
              ! Select hyperslab
              call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                                         g_start, s_size, ierr,       &
                                         stride, block)
              ! Write the data
              call h5dwrite_f(d_id, H5T_IEEE_F64LE, ld, s_size, ierr,         &
                              file_space_id=filespace, mem_space_id=memspace, &
                              xfer_prp=x_id)
              ! Next process hyperslab
              g_start(3) = g_start(3) + np
           enddo

           ! Close dataset
           call h5dclose_f(d_id, ierr)
        end do

        ! Write the data
        !call h5dwrite_f(d_id, H5T_IEEE_F64LE, ld, s_size, ierr,         &
        !                file_space_id=filespace, mem_space_id=memspace, &
        !                xfer_prp=x_id)

        ! ----------------------------------------------------------------------

        call mpi_barrier(MPI_COMM_3D, ierr)

        if (id .eq. 0) then
           call system_clock(clock_counts(3), clock_rate)
        endif
        
        if (allocated(ld)) then
                deallocate(ld)
        endif

        write(0,*) "Process:", id, ", chckpoint:", 10

        ! Close everything and exit
        !call h5dclose_f(d_id, ierr)
        call h5sclose_f(filespace, ierr)
        call h5sclose_f(memspace, ierr)
        call h5pclose_f(c_id, ierr)
        call h5pclose_f(x_id, ierr)
        call h5pclose_f(p_id, ierr)
        call h5fclose_f(f_id, ierr)
        call h5close_f(ierr)

        ! ----------------------------------------------------------------------

        call mpi_barrier(MPI_COMM_3D, ierr)

        if (id .eq. 0) then
           call system_clock(clock_counts(4), clock_rate)
        endif
 
        if (id .eq. 0) then
           print *, "elapsed time(create): ", real(clock_counts(2)-clock_counts(1))/real(clock_rate)
           print *, "elapsed time( write): ", real(clock_counts(3)-clock_counts(2))/real(clock_rate)
           ! Calculate statistics
           stat_total = real((real(g_N)/1.e3)**2 * nsets * 8) ! (MBytes)
           stat_time_create = real(clock_counts(2)-clock_counts(1))/real(clock_rate) ! (sec)
           stat_rate_create = stat_total/stat_time_create
           stat_time_wr = real(clock_counts(3)-clock_counts(2))/real(clock_rate) ! (sec)
           stat_rate_wr = stat_total/stat_time_wr
           print *, "saved ", stat_total, "(MB), rate(create): ", stat_rate_create, "(MB/s)", &
                ", rate(write): ", stat_rate_wr, "(MB/s)"
        endif
       
        call mpi_finalize(ierr)
end program hdf_pwrite3dc
