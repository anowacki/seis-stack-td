program stack_fk_prog
   ! Read in a list of files from stdin and save the fk stack to a NetCDF file
   ! if desired

   use f90sac
   use stack
   use, intrinsic :: iso_c_binding, only: C_FLOAT
   use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, &
      nf90_put_var, nf90_close, NF90_CLOBBER, NF90_FLOAT, nf90_noerr, nf90_strerror, &
      NF90_GLOBAL, nf90_put_att

   implicit none

   ! Maximum number of files
   integer, parameter :: nmax = 1000
   character(len=250) :: infile, file, type = 'linear'
   type(SACtrace) :: s(nmax)
   real(C_FLOAT) :: picks(nmax), pick_temp, t1, t2, sxmin, sxmax, symin, symax, ds
   real(C_FLOAT), allocatable :: fk(:,:), ux(:), uy(:)
   real(C_FLOAT) :: lon, lat, gcarc, baz, evdp
   integer :: i, j, n, iostat, n_stack, nx, ny
   logical :: use_picks = .false., write_ncf = .false.
#ifdef USE_XAPIIR
   real(C_FLOAT) :: corner1, corner2
   integer :: poles, passes
   character(len=4) :: filter = 'none'
#endif

   call get_args

   ! Read SAC files
   iostat = 0
   n = 0
   do while (iostat == 0)
      if (use_picks) then
         read(*,*,iostat=iostat) infile, pick_temp
      else
         read(*,*,iostat=iostat) infile
      endif
      if (iostat > 0) then
         write(0,'(a,i0.1)') 'stack_fk: Error: Problem getting file (and pick) from stdin, line ', &
            n + 1
         error stop
      endif
      if (iostat < 0) exit
      if (infile == "") cycle
      n = n + 1
      if (n > nmax) then
         write(0,'(2(a,i0.1))') 'stack_fk: Error: Number of traces (', n, &
            ') greater than precompiled limits (', nmax, ')'
         error stop
      endif
      call f90sac_readtrace(infile, s(n))
      call filter_trace(s(n))
      if (use_picks) picks(n) = pick_temp
   enddo

   ! Make vespagram
   if (use_picks) then
      call stack_fk(s(1:n), t1, t2, sxmin, sxmax, symin, symax, ds, fk, ux=ux, uy=uy, &
         type=type, n=n_stack, pick=picks(1:n))
   else
      call stack_fk(s(1:n), t1, t2, sxmin, sxmax, symin, symax, ds, fk, ux=ux, uy=uy, &
         type=type, n=n_stack)
   endif

   nx = size(ux)
   ny = size(uy)

   ! Write out stack
   if (write_ncf) then
      ! Get info for comments
      call stack_array_geography(s(1:n), lon, lat, gcarc, baz, evdp)
      call write_netcdf_file
   else
      do i = 1, ny
         do j = 1, nx
            write(*,*) ux(j), uy(i), fk(j, i)
         enddo
      enddo
   endif

contains
   subroutine usage(error)
      logical, optional, intent(in) :: error
      integer :: lu_usage
      lu_usage = 6 ! stdout
      if (present(error)) then
         if (error) lu_usage = 0
      endif
      write(lu_usage,'(a)') &
         'Usage: stack_fk (options) [t1] [t2] [sxmin] [sxmax] [symin] [symax] [ds] < (list of SAC files (pick times))', &
         '   Read a list of SAC files on stdin and write a frequency-wavenumber stack of', &
         '   (ux, uy, power) triplets to stdout', &
         'Arguments:', &
         '   t1, t2       : Time window start and end relative to O marker (s)', &
         '   sxmin, sxmax : Range of slowness to search in east direction (s/deg)', &
         '   symin, symax : Range of slowness to search in north direction (s/deg)', &
         '   ds           : Slowness increment (s/deg)', &
         'Options:', &
#ifdef USE_XAPIIR
         '   -bp [c1] [c2] [poles] [npasses] :', &
         '                  Perform band-pass filter between corners <c1> and <c2>,', &
         '                  specifying number of poles and passes, before stacking', &
         '   -hp [corner] [poles] [npasses]:', &
         '                  High-pass filter before stacking', &
         '   -lp [corner] [poles] [npasses]:', &
         '                  Low-pass filter before stacking', &
#endif
         '   -n [n]       : For nthroot or phaseweighted, use power n', &
         '   -o [file]    : Output a NetCDF file instead of writing to stdout', &
         '   -p           : Align on picks in column 2 of input', &
         '   -type [type] : Choose from the following stack types:', &
         '        [l]inear, [p]haseweighted, [n]throot', &
         '   -h, --help   : Print this usage message to stdout and exit'
         if (present(error)) then
            if (error) error stop
         else
            stop
         endif
   end subroutine usage

   subroutine get_args
      integer :: iarg, narg
      character(len=250) :: arg
      narg = command_argument_count()
      if (narg == 1) then
         call get_command_argument(1, arg)
         if (trim(arg) == '-h' .or. trim(arg) == '--help') then
            call usage()
         else
            call usage(error=.true.)
         endif
      elseif (narg < 7) then
         call usage(error=.true.)
      endif
      iarg = 1
      do while (iarg < narg - 6)
         call get_command_argument(iarg, arg)
         select case (arg)
#ifdef USE_XAPIIR
            case ('-bp')
               filter = 'bp'
               call get_command_argument(iarg + 1, arg)
               read(arg,*) corner1
               call get_command_argument(iarg + 2, arg)
               read(arg,*) corner2
               call get_command_argument(iarg + 3, arg)
               read(arg,*) poles
               call get_command_argument(iarg + 4, arg)
               read(arg,*) passes
               iarg = iarg + 5
            case ('-lp')
               filter = 'lp'
               call get_command_argument(iarg + 1, arg)
               read(arg,*) corner1
               call get_command_argument(iarg + 2, arg)
               read(arg,*) poles
               call get_command_argument(iarg + 3, arg)
               read(arg,*) passes
               iarg = iarg + 4
            case ('-hp')
               filter = 'hp'
               call get_command_argument(iarg + 1, arg)
               read(arg,*) corner1
               call get_command_argument(iarg + 2, arg)
               read(arg,*) poles
               call get_command_argument(iarg + 3, arg)
               read(arg,*) passes
               iarg = iarg + 4
#endif
            case ('-h', '--help')
               call usage(error=.false.)
            case ('-n')
               call get_command_argument(iarg + 1, arg)
               read(arg,*) n_stack
               iarg = iarg + 2
            case ('-o')
               call get_command_argument(iarg + 1, file)
               write_ncf = .true.
               iarg = iarg + 2
            case ('-p')
               use_picks = .true.
               iarg = iarg + 1
            case ('-type')
               call get_command_argument(iarg + 1, type)
               iarg = iarg + 2
            case default
               write(0,'(a)') 'stack_fk: Error: Unrecognised option "' // trim(arg) // '"'
               error stop
         end select
      enddo
      call get_command_argument(narg - 6, arg)
      read(arg,*) t1
      call get_command_argument(narg - 5, arg)
      read(arg,*) t2
      call get_command_argument(narg - 4, arg)
      read(arg,*) sxmin
      call get_command_argument(narg - 3, arg)
      read(arg,*) sxmax
      call get_command_argument(narg - 2, arg)
      read(arg,*) symin
      call get_command_argument(narg - 1, arg)
      read(arg,*) symax
      call get_command_argument(narg, arg)
      read(arg,*) ds
   end subroutine get_args

   subroutine check_ncf(ierror)
      integer, intent(in) :: ierror
      if (ierror /= nf90_noerr) then
         write(0,'(a)') 'stack_fk: Error: Problem creating NetCDF grid: ' &
            // trim(nf90_strerror(ierror))
         error stop
      endif
   end subroutine check_ncf

   subroutine write_netcdf_file
      integer :: ncid, ux_dimid, uy_dimid, ux_varid, uy_varid, a_varid
      character(len=250) :: comment

      call check_ncf(nf90_create(trim(file), NF90_CLOBBER, ncid))

      ! Define dimensions
      call check_ncf(nf90_def_dim(ncid, 'ux', nx, ux_dimid))
      call check_ncf(nf90_def_dim(ncid, 'uy', ny, uy_dimid))
      ! Set variables for coordinates
      call check_ncf(nf90_def_var(ncid, 'ux', NF90_FLOAT, ux_dimid, ux_varid))
      call check_ncf(nf90_def_var(ncid, 'uy', NF90_FLOAT, uy_dimid, uy_varid))
      call check_ncf(nf90_put_att(ncid, ux_varid, 'units', 's/deg'))
      call check_ncf(nf90_put_att(ncid, uy_varid, 'units', 's/deg'))
      call check_ncf(nf90_put_att(ncid, ux_varid, 'long_name', 'Slowness E'))
      call check_ncf(nf90_put_att(ncid, uy_varid, 'long_name', 'Slowness N'))
      call check_ncf(nf90_put_att(ncid, ux_varid, 'actual_range', [sxmin, sxmax]))
      call check_ncf(nf90_put_att(ncid, uy_varid, 'actual_range', [symin, symax]))
      call check_ncf(nf90_def_var(ncid, 'power', NF90_FLOAT, [ux_dimid, uy_dimid], &
         a_varid))
      ! Set variables for amplitude
      call check_ncf(nf90_put_att(ncid, a_varid, 'long_name', 'Normalised power'))
      call check_ncf(nf90_put_att(ncid, a_varid, 'actual_range', [minval(fk), maxval(fk)]))

      ! Add comments about data
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'title', &
         'FK response created by stack_fk'))
      write(comment, '(5(a,f0.6))') 'Stack_type: ' // trim(type) // ' Mean_lon: ', lon, &
         ' Mean_lat: ', lat, ' Mean_gcarc: ', gcarc, ' Mean_baz: ', baz, ' evdp: ', evdp
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(comment)))
      ! Finish data description
      call check_ncf(nf90_enddef(ncid))

      ! Fill structure
      call check_ncf(nf90_put_var(ncid, ux_varid, ux))
      call check_ncf(nf90_put_var(ncid, uy_varid, uy))
      call check_ncf(nf90_put_var(ncid, a_varid, fk))

      ! Finalise file
      call check_ncf(nf90_close(ncid))
   end subroutine write_netcdf_file

   subroutine filter_trace(trace)
      type(SACtrace), intent(inout) :: trace
#ifdef USE_XAPIIR
      select case (filter)
         case ('none')
            return
         case ('bp')
            call f90sac_bandpass_bu(trace, corner1, corner2, npoles=poles, npasses=passes)
         case ('hp')
            call f90sac_highpass_bu(trace, corner1, npoles=poles, npasses=passes)
         case ('lp')
            call f90sac_lowpass_bu(trace, corner1, npoles=poles, npasses=passes)
      end select
#endif
   end subroutine filter_trace

end program stack_fk_prog