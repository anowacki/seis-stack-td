program stack_vespa_prog
   ! Read in a list of files from stdin and save the linear stack to a SAC file

   use f90sac
   use stack
   use, intrinsic :: iso_c_binding, only: C_FLOAT

   implicit none

   ! Maximum number of files
   integer, parameter :: nmax = 1000
   character(len=250) :: infile, outfile, type = 'linear'
   type(SACtrace) :: s(nmax)
   real(C_FLOAT) :: picks(nmax), pick_temp, t1, t2, s1, s2, ds
   real(C_FLOAT), allocatable :: vesp(:,:), time(:), slowness(:)
   integer :: i, j, n, iostat, n_stack, nt, ns
   logical :: use_picks = .false.

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
         write(0,'(a,i0.1)') 'stack_sum: Error: Problem getting file (and pick) from stdin, line ', &
            n + 1
         error stop
      endif
      if (iostat < 0) exit
      if (infile == "") cycle
      n = n + 1
      if (n > nmax) then
         write(0,'(2(a,i0.1))') 'stack_sum: Error: Number of traces (', n, &
            ') greater than precompiled limits (', nmax, ')'
         error stop
      endif
      call f90sac_readtrace(infile, s(n))
      if (use_picks) picks(n) = pick_temp
   enddo

   ! Make vespagram
   if (use_picks) then
      call stack_vespa_slow(s(1:n), t1, t2, s1, s2, ds, vesp, time=time, &
         slowness=slowness, type=type, n=n_stack, pick=picks(1:n))
   else
      call stack_vespa_slow(s(1:n), t1, t2, s1, s2, ds, vesp, time=time, &
         slowness=slowness, type=type, n=n_stack)
   endif

   ! Write out stack
   do i = 1, size(slowness)
      do j = 1, size(time)
         write(*,*) time(j), slowness(i), vesp(j, i)
      enddo
   enddo

contains
   subroutine usage
      write(0,'(a)') &
         'Usage: stack_vespa [s1] [s2] [ds] [t1] [t2] (options) < (list of SAC files (pick times))', &
         '   Read a list of SAC files on stdin and write a vespagram', &
         '   of (t,s,amp) triplets to stdout', &
         'Arguments:', &
         '   s1, s2, ds   : Min and max slownesses and slowness increment (s/deg)', &
         '   t1, t2       : Time window start and end relative to O marker (s)', &
         'Options:', &
         '   -n [n]       : For nthroot or phaseweighted, use root n', &
         '   -p           : Use picks in column 2 of input', &
         '   -t [t1] [t2] : Use range of times (relative to picks if provided)', &
         '   -type [type] : Choose from the following stack types:', &
         '        [l]inear, [p]haseweighted, [n]throot'
      error stop
   end subroutine usage

   subroutine get_args
      integer :: iarg, narg
      character(len=250) :: arg
      narg = command_argument_count()
      if (narg < 5) call usage
      iarg = 1
      do while (iarg < narg - 4)
         call get_command_argument(iarg, arg)
         select case (arg)
            case ('-n')
               call get_command_argument(iarg + 1, arg)
               read(arg,*) n_stack
               iarg = iarg + 2
            case ('-p')
               use_picks = .true.
               iarg = iarg + 1
            case ('-type')
               call get_command_argument(iarg + 1, type)
               iarg = iarg + 2
            case default
               write(0,'(a)') 'stack_vespa: Error: Unrecognised option "' // trim(arg) // '"'
               error stop
         end select
      enddo
      call get_command_argument(narg - 4, arg)
      read(arg,*) s1
      call get_command_argument(narg - 3, arg)
      read(arg,*) s2
      call get_command_argument(narg - 2, arg)
      read(arg,*) ds
      call get_command_argument(narg - 1, arg)
      read(arg,*) t1
      call get_command_argument(narg, arg)
      read(arg,*) t2
   end subroutine get_args

end program stack_vespa_prog