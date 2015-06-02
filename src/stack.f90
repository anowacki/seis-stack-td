!===============================================================================
module stack
!===============================================================================
!
! Routines for stacking seismograms stored in SAC format.
!
! Andy Nowacki
! a.nowacki@leeds.ac.uk
!

   use, intrinsic :: iso_c_binding, only: C_FLOAT, C_DOUBLE
   use f90sac
   ! Use the OpenMP module if needed
!$ use omp_lib

   implicit none

   private

   ! Precision selectors--based on C types because we deal with SAC data
   integer, parameter :: r4 = C_FLOAT, r8 = C_DOUBLE
   integer, parameter :: rs = r4

   ! I/O
   integer, parameter :: lu_stdin = 5, lu_stdout = 6, lu_stderr = 0, &
                         lu_in = 10, lu_out = 11

   ! Useful constants
   real(rs), parameter :: pi = real(4.d0*atan2(1.d0, 1.d0), rs)
   real(rs), parameter :: deg2km = 111.19_rs

   ! Length parameters
   integer, parameter :: STACK_CHAR_LEN = 30

   interface stack_allocate
      module procedure :: stack_allocate_r4_1d
   end interface stack_allocate

   public :: &
      stack_sum

contains

!===============================================================================
subroutine stack_sum(s, out, t1, t2, pick, type, n)
!===============================================================================
! Take an array of SAC traces and return the summed trace between times t1 and t2
! (if given), or just the whole trace.
!  INPUT:
!     s(:)   : SACtrace array of input traces
!  INPUT (OPTIONAL):
!     t1, t2 : Start and stop times for output stack.  Relative to O time, unless
!              pick array is supplied
!     pick(:): Array of pick times relative to which stack is performed.  Must be
!              an array of size(s)
!     type   : Character containing one of:
!                 '[l]inear'        : Simple linear stack [default]
!                 '[n]throot'       : N-th root stack
!                 '[p]haseweighted' : Phase-weighted stack
!     n      : If 'nthroot' or 'phaseweighted' are chosed, this option determined
!              the power for the stack.  [Default 2]
!  OUTPUT:
!     out(:) : SACtrace containing the stack
!
   type(SACtrace), intent(in) :: s(:)
   type(SACtrace), intent(inout) :: out
   real(rs), intent(in), optional :: t1, t2, pick(:)
   character(len=*), intent(in), optional :: type
   integer, intent(in), optional :: n
   real(rs) :: w1, w2
   integer :: npts
   real(rs), allocatable :: pick_in(:)
   character(len=STACK_CHAR_LEN) :: type_in
   integer :: n_in

   ! Basic checks on input
   call stack_check(s)

   ! Check input
   w1 = maxval(s%b)
   w2 = minval(s%e)
   if (present(t1)) w1 = t1
   if (present(t2)) w2 = t2
   if (w2 <= w1) call stack_error('stack_sum: Start time must be before end time')
   if (present(pick)) then
      if (size(pick) /= size(s)) &
         call stack_error('stack_sum: Pick array must be same length as number of traces')
      if (any(pick > s%e .or. pick < s%b)) &
         call stack_error('stack_sum: Picks must not be outside trace')
      if (any(pick + t2 > s%e .or. pick - t1 < s%b)) &
         call stack_error('stack_sum: Window around picks is outside data range')
   endif
   type_in = 'linear'
   if (present(type)) type_in = type
   n_in = 2
   if (present(n)) n_in = n

   ! Make space for stack
   npts = int((w2 - w1)/s(1)%delta) + 1
   call f90sac_newtrace(npts, s(1)%delta, out)
   out%b = w1
   out%e = w2

   ! Create internal pick array and populate
   allocate(pick_in(npts))
   if (present(pick)) then
      pick_in = pick
   else
      pick_in = 0._rs
   endif

   ! Peform the stack
   select case (type_in(1:1))
      case ('l', 'L')
         call stack_sum_linear(s, pick_in, w1, w2, npts, out%trace)
      case ('n', 'N')
         call stack_sum_nthroot(s, pick_in, w1, w2, n_in, npts, out%trace)
      case ('p', 'P')
         call stack_sum_phaseweighted(s, pick_in, w1, w2, n_in, npts, out%trace)
      case default
         call stack_error('stack_sum: Unimplemented stack type "'//trim(type_in)//'"')
   end select

   deallocate(pick_in)

end subroutine stack_sum
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_sum_linear(s, pick, w1, w2, npts, out)
!===============================================================================
! Perform a linear stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     w1, w2  : Start and stop times of stack, relative to pick(:)
!     npts    : Length of stacked trace
!  OUTPUT:
!     out(npts) : Array containing points of the stack
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: pick(:), w1, w2
   integer, intent(in) :: npts
   real(rs), intent(out) :: out(npts)
   integer :: i, j, iw1, iw2

   out = 0._rs
!$omp parallel do default(none) shared(s, npts, out, pick, w1, w2) &
!$omp    private(i, j, iw1, iw2)
   do i = 1, size(s)
      iw1 = nint((pick(i) + w1 - s(i)%b)/s(1)%delta) + 1
      iw2 = iw1 + npts - 1
      do j = iw1, iw2
         out(j-iw1+1) = out(j-iw1+1) + s(i)%trace(j)
      enddo
   enddo
!$omp end parallel do
   out = out/real(size(s))

end subroutine stack_sum_linear
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_sum_nthroot(s, pick, w1, w2, n, npts, out)
!===============================================================================
! Perform an N-th root stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     w1, w2  : Start and stop times of stack, relative to pick(:)
!     n       : Power to which to raise stack
!     npts    : Length of stacked trace
!  OUTPUT:
!     out(npts) : Array containing points of the stack
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: pick(:), w1, w2
   integer, intent(in) :: n, npts
   real(rs), intent(out) :: out(npts)
   integer :: i, j, iw1, iw2

   out = 0._rs
!$omp parallel do default(none) shared(s, npts, n, out, pick, w1, w2) &
!$omp    private(i, j, iw1, iw2)
   do i = 1, size(s)
      iw1 = nint((pick(i) + w1 - s(i)%b)/s(1)%delta) + 1
      iw2 = iw1 + npts - 1
      do j = iw1, iw2
         out(j-iw1+1) = out(j-iw1+1) + &
            sign(abs(s(i)%trace(j))**(1._rs/real(n)), s(i)%trace(j))
      enddo
   enddo
!$omp end parallel do
   out = sign(out**n, out)/real(size(s))
end subroutine stack_sum_nthroot
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_sum_phaseweighted(s, pick, w1, w2, nu, npts, out)
!===============================================================================
! Perform a phase-weighted stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     w1, w2  : Start and stop times of stack, relative to pick(:)
!     nu      : Power to which to raise coherency
!     npts    : Length of stacked trace
!  OUTPUT:
!     out(npts) : Array containing points of the stack
!
! The phase-weighted stack is given by:
!
!     v(t) = (1/N) \Sigma_{j=1}^N s_j(t) abs[(1/N) \Sigma_{k=1}^N exp(i\Phi_k(t))]^\nu
!
! (Eq. 4 of Schimmel & Paulssen, GJI, 1997)
! with N = number of traces, s_j(t) is the time-domain signal for the jth trace,
! i is sqrt(-1), \Phi_k(t) is the instantaneous phase of the kth trace, and \nu
! is the power parameter.  \Phi(t) is given by arg(s(t) + i*H(t)), where H(t)
! is the Hilbert transform of the trace s(t).
! We therefore must construct the Hilbert trace for each trace
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: pick(:), w1, w2
   integer, intent(in) :: nu, npts
   real(rs), intent(out) :: out(npts)
   real(rs), dimension(npts,size(s)) :: trace, hilbert, phase
   real(rs) :: weight
   integer :: n, i, j, iw1, iw2

   n = size(s)

   ! Create Hilbert traces
   do i = 1, n
      iw1 = nint((pick(i) + w1 - s(i)%b)/s(1)%delta) + 1
      iw2 = iw1 + npts - 1
      trace(:,i) = s(i)%trace(iw1:iw2)
   enddo
   hilbert = stack_hilbert(trace, n, npts)

   ! Create array of phase at each point for each trace
   phase = atan2(hilbert, trace)

   ! Do the stack
   out = 0._rs
!$omp parallel do default(none) shared(phase, n, out, trace, nu, npts) &
!$omp    private(i, j, weight)
   do i = 1, n
      do j = 1, npts
         weight = sum(abs(phase(j,:))/real(n))
         out(j) = out(j) + trace(j,i)*weight**nu
      enddo
   enddo
!$omp end parallel do
   out = out/real(n)

end subroutine stack_sum_phaseweighted
!-------------------------------------------------------------------------------

!===============================================================================
function stack_hilbert(y, n, npts) result(h)
!===============================================================================
! Return the hilbert transformed version of y(npts,n), which contains n traces
! of npts length.
! This routine uses FFTW3, and for various reasons it is impractical to assume
! that it has been compiled correctly to handle multithreading, so we do not use
! OpenMP here.
! Note that because we are using the real-to-complex transforms, we only have the
! positive frequencies in c, so we can multiply everything by -i.
!
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'
   real(C_FLOAT), intent(in) :: y(npts,n)
   integer, intent(in) :: n, npts
   real(rs) :: h(npts,n), h_temp(npts), y_temp(npts)
   complex(C_FLOAT) :: c(npts/2+1)
   type(C_PTR) :: plan_fwd, plan_rev
   integer(C_INT), parameter :: fftw_plan_type = FFTW_ESTIMATE
   integer :: i

   ! Create plans, which we can reuse if we use the same arrays
   call sfftw_plan_dft_r2c_1d(plan_fwd, int(npts, kind=C_INT), y_temp, c, fftw_plan_type)
   call sfftw_plan_dft_c2r_1d(plan_rev, int(npts, kind=C_INT), c, h_temp, fftw_plan_type)
   ! Make Hilbert traces by multiplying the frequency-domain complex traces by -i
   do i = 1, n
      y_temp = y(:,i)
      call sfftw_execute_dft_r2c(plan_fwd, y_temp, c)
      c = c*complex(0._rs, -1._rs)
      call sfftw_execute_dft_c2r(plan_rev, c, h_temp)
      h(:,i) = h_temp
   enddo
   ! Normalise values back to original amplitude
   h = h/real(npts)
   ! Deassociate pointers
   call sfftw_destroy_plan(plan_fwd)
   call sfftw_destroy_plan(plan_rev)
end function stack_hilbert
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_check(s)
!===============================================================================
! Check that an array of SAC traces conforms to the criteria necessary for use
! in routines in this module.
!
   type(SACtrace), intent(in) :: s(:)
   integer :: n

   n = size(s)
   if (n == 1) then
      call stack_warning('stack_check: Only one trace in stack')
      return
   endif
   if (any(s(2:n)%delta /= s(1)%delta)) &
      call stack_error('stack_check: All traces must have same delta')
end subroutine stack_check
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_error(string)
!===============================================================================
   character(len=*), intent(in) :: string
   write(lu_stderr,'(a)') 'Error: '// string
   error stop
end subroutine stack_error
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_warning(string)
!===============================================================================
   character(len=*), intent(in) :: string
   write(lu_stderr,'(a)') 'Warning: '// string
end subroutine stack_warning
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_allocate_r4_1d(a, n)
!===============================================================================
   real(r4), intent(inout), allocatable, dimension(:) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) == n) return
      deallocate(a)
      allocate(a(n))
      return
   else
      allocate(a(n))
   endif
end subroutine stack_allocate_r4_1d
!-------------------------------------------------------------------------------

end module stack
!-------------------------------------------------------------------------------
