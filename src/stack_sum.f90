program stack_sum_prog
   ! Read in a list of files from stdin and save the linear stack to a SAC file

   use f90sac
   use stack
   use, intrinsic :: iso_c_binding, only: C_FLOAT

   ! Maximum number of files
   integer, parameter :: nmax = 1000
   character(len=250) :: infile, outfile
   type(SACtrace) :: s(nmax), sum
   real(C_FLOAT) :: picks(nmax), pick_temp
   integer :: iostat, n
   logical :: window_set = .false., use_picks = .false.

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

   ! Stack over all the common data points by default
   if (.not.window_set) then
      t1 = maxval(s(1:n)%b)
      t2 = minval(s(1:n)%e)
   endif

   ! Perform stack
   if (use_picks) then
      call stack_sum(s(1:n), sum, t1=t1, t2=t2, pick=picks(1:n))
   else
      call stack_sum(s(1:n), sum, t1=t1, t2=t2)
   endif

   ! Write out stack
   if (outfile == '-') then
      do i = 1, sum%npts
         write(*,*) sum%b + real(i-1)*sum%delta, sum%trace(i)
      enddo
   else
      call f90sac_writetrace(outfile, sum)
   endif

contains
   subroutine usage
      write(0,'(a)') &
         'Usage: stack_sum (options) [outfile] < (list of SAC files (pick times))', &
         '   Read a list of SAC files on stdin and save the linear stack', &
         '   to <outfile>.', &
         '   If <outfile> is ''-'', write (t,amp) pairs to stdout', &
         'Options:', &
         '   -p           : Use picks in column 2 of input', &
         '   -t [t1] [t2] : Output stack time, relative to'
      error stop
   end subroutine usage

   subroutine get_args
      integer :: iarg, narg
      character(len=250) :: arg
      narg = command_argument_count()
      if (narg < 1) call usage
      iarg = 1
      do while (iarg < narg)
         call get_command_argument(iarg, arg)
         select case (arg)
            case ('-p')
               use_picks = .true.
               iarg = iarg + 1
            case ('-t')
               window_set = .true.
               call get_command_argument(iarg + 1, arg)
               read(arg,*) t1
               call get_command_argument(iarg + 2, arg)
               read(arg,*) t2
               iarg = iarg + 3
            case default
               write(0,'(a)') 'stack_sum: Error: Unrecognised option "' // trim(arg) // '"'
               error stop
         end select
      enddo
      call get_command_argument(narg, outfile)
   end subroutine get_args

end program stack_sum_prog