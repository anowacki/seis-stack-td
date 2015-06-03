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
   real(rs), parameter :: todeg = 180._rs/pi, torad = pi/180._rs, deg2km = 111.19_rs
   real(rs), parameter :: R_EARTH_KM = 6371._rs
   real(rs), parameter :: s_deg2s_km = 180._rs/(pi*R_EARTH_KM)

   ! Length parameters
   integer, parameter :: STACK_CHAR_LEN = 30

   interface stack_allocate
      module procedure :: stack_allocate_r4_1d, stack_allocate_r4_2d
   end interface stack_allocate

   public :: &
      stack_fk, stack_sum, stack_vespa_slow

contains

!===============================================================================
subroutine stack_fk(s, t1, t2, smax, ds, out, pick, type, n, u)
!===============================================================================
! Perform f-k analysis on an array of SAC traces.  Search over a range of
! slowness vectors with maximum magnitude smax in ds increments.
!  INPUT:
!     s(:)    : SACtrace array of input traces
!     t1, t2  : Start and stop time window in s
!     smax    : Maximum slowness magnitude in s/deg
!     ds      : Slowness increment in s/deg
!  INPUT (OPTIONAL):
!     pick(:) : Array of pick times relative to which vespagram is made
!     type    : Character determining stack type.  See stack_sum() for details.
!     n       : Power of the stack; see stack_sum() for details.
!  OUTPUT:
!     out(:,:): Allocatable array of size(n,n), where n = 2*smax/ds + 1,
!               containing the normalised beam power at each point over the time
!               window of interest, which is the summed squared amplitudes.
!  OUTPUT (OPTIONAL):
!     u(:)    : Allocatable array containing the values of slowness in both
!               x and y directions.
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: t1, t2, smax, ds
   real(rs), allocatable, intent(inout) :: out(:,:)
   real(rs), intent(in), optional, dimension(size(s)) :: pick
   character(len=*), intent(in), optional :: type
   integer, intent(in), optional :: n
   real(rs), intent(out), allocatable, optional :: u(:)
   type(SACtrace) :: trace
   real(rs), dimension(size(s)) :: pick_in, x, y, r, phi, delay
   character(len=STACK_CHAR_LEN) :: type_in
   real(rs) :: lon, lat, ux, uy
   integer :: i, ix, iy, n_in, nu

   ! Basic checks on input
   call stack_check(s)

   if (t2 <= t1) call stack_error('stack_fk: Start time must be before end time')
   if (present(pick)) then
      if (size(pick) /= size(s)) &
         call stack_error('stack_fk: Pick array must be same length as number of traces')
      if (any(pick > s%e .or. pick < s%b)) &
         call stack_error('stack_fk: Picks must not be outside trace')
      if (any(pick + t2 > s%e .or. pick - t1 < s%b)) &
         call stack_error('stack_fk: Window around picks is outside data range')
   endif
   if (any(s%stlo == SAC_rnull .or. s%stla == SAC_rnull)) &
      call stack_error('stack_fk: All stations must headers STLO and STLA filled')

   ! Change defaults if necessary
   type_in = 'linear'
   if (present(type)) type_in = type
   n_in = 2
   if (present(n)) n_in = n
   pick_in = 0._rs
   if (present(pick)) pick_in = pick

   ! Calculate mean coordinates and get array of cartesian coordinates in km,
   ! relative to the mean point at (0,0)
   call stack_calc_station_coords(s%stlo, s%stla, lon, lat, x, y, r, phi)

   ! Calculate size of output and allocate
   nu = int(2.*smax/ds) + 1
   call stack_allocate(out, nu, nu)

   ! Create stack at each slowness point
   do iy = 1, nu
      uy = real(iy-1)*ds - smax
      uy = uy*s_deg2s_km
      do ix = 1, nu
         ux = real(ix-1)*ds - smax
         ux = ux*s_deg2s_km
         delay = pick_in - stack_station_delay(x, y, ux, uy)
         call stack_sum(s, trace, t1=t1, t2=t2, delay=delay, type=type_in, n=n_in)
         out(ix,iy) = sum(trace%trace**2)
      enddo
   enddo

   out = out/maxval(out)

   ! If requested, fill in the slowness axis
   if (present(u)) then
      call stack_allocate(u, nu)
      u = [(real(i-1)*ds - smax, i = 1, nu)]
   endif

end subroutine stack_fk
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_vespa_slow(s, t1, t2, s1, s2, ds, out, pick, type, n, time, &
                            slowness)
!===============================================================================
! Take an array of SAC traces and return the slowness vespagram.  This is done
! relative to the geographic mean point of the array, using the backazimuth from
! this point to the event which must be in at least the first trace's header.
!  INPUT:
!     s(:)    : SACtrace array of input traces
!     t1, t2  : Start and stop time window in s
!     s1, s2  : Start and end slownesses in s/deg
!     ds      : Slowness increment in s/deg
!  INPUT (OPTIONAL):
!     pick    : Array of pick times relative to which vespagram is made.
!               must be an array of size(s)
!     type    : Character determining stack type.  See stack_sum() for details.
!     n       : Power of the stack; see stack_sum() for details.
!  OUTPUT:
!     out(:,:): Alloctable array of size((t2 - t1)/s(1)%delta+1, (s2-s1)/ds+1)
!               containing the normalised beam power at each (s,t) point.
!  OUTPUT (OPTIONAL):
!     time(:) : Alloctable array of size((t2-t1)/s(1)%delta+1) containing the
!               values of the times in the vespa grid (s)
!     slow(:) : Alloctable array of size((s2-s1)/ds+1) containing the slownesses
!               in the vespa grid (s/deg)
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: t1, t2, s1, s2, ds
   real(rs), allocatable, intent(inout) :: out(:,:)
   real(rs), intent(in), optional :: pick(size(s))
   character(len=*), intent(in), optional :: type
   integer, intent(in), optional :: n
   real(rs), allocatable, intent(inout), optional :: time(:), slowness(:)
   type(SACtrace) :: trace
   real(rs) :: pick_in(size(s)), delay(size(s))
   character(len=STACK_CHAR_LEN) :: type_in
   integer :: n_in, ns, nt, i, j
   real(rs) :: lon, lat, x(size(s)), y(size(s)), r(size(s)), phi(size(s))
   real(rs) :: slow, ref_az, ref_baz, ref_delta, ux, uy

   ! Basic checks on input
   call stack_check(s)

   if (t2 <= t1) call stack_error('stack_vespa: Start time must be before end time')
   if (present(pick)) then
      if (size(pick) /= size(s)) &
         call stack_error('stack_vespa: Pick array must be same length as number of traces')
      if (any(pick > s%e .or. pick < s%b)) &
         call stack_error('stack_vespa: Picks must not be outside trace')
      if (any(pick + t2 > s%e .or. pick - t1 < s%b)) &
         call stack_error('stack_vespa: Window around picks is outside data range')
   endif
   if (s(1)%evlo == SAC_rnull .or. s(1)%evla == SAC_rnull) &
      call stack_error('stack_vespa: Event location headers must be in first trace')
   if (any(s%stlo == SAC_rnull .or. s%stla == SAC_rnull)) &
      call stack_error('stack_vespa: All stations must headers STLO and STLA filled')

   ! Change defaults if necessary
   type_in = 'linear'
   if (present(type)) type_in = type
   n_in = 2
   if (present(n)) n_in = n
   pick_in = 0._rs
   if (present(pick)) pick_in = pick

   ! Calculate mean coordinates and get array of cartesian coordinates in km,
   ! relative to the mean point at (0,0)
   call stack_calc_station_coords(s%stlo, s%stla, lon, lat, x, y, r, phi)
   ref_delta = stack_delta(s(1)%evlo, s(1)%evla, lon, lat)
   ref_baz = stack_azimuth(lon, lat, s(1)%evlo, s(1)%evla)
   ref_az = ref_baz + 180._rs

   ! Allocate vespagram array
   ns = int((s2 - s1)/ds) + 1
   nt = int((t2 - t1)/s(1)%delta) + 1
   call stack_allocate(out, nt, ns)

   ! Calculate stack at each slowness point, filling the vespagram
   do j = 1, ns
      slow = s1 + real(j-1)*ds ! s/deg
      slow = slow*s_deg2s_km ! s/km
      ux = slow*sin(torad*ref_az)
      uy = slow*cos(torad*ref_az)
      ! Create the delays at each station for this slowness
      delay = pick_in + stack_station_delay(x, y, ux, uy)
      call stack_sum(s, trace, t1=t1, t2=t2, delay=delay, type=type_in, n=n_in)
      out(:,j) = trace%trace
   enddo

   ! If desired, output time and slowness arrays
   if (present(time)) then
      call stack_allocate(time, nt)
      time = [(t1 + real(i-1)*s(1)%delta, i = 1, nt)]
   endif
   if (present(slowness)) then
      call stack_allocate(slowness, ns)
      slowness = [(s1 + real(i-1)*ds, i = 1, ns)]
   endif

   ! Normalise vespagram
   out = out/maxval(abs(out))

end subroutine stack_vespa_slow
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_sum(s, out, t1, t2, pick, delay, type, n)
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
!     delay(:):Array of delay times by which to shuffle traces.  Must be array of
!              size(s)
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
   real(rs), intent(in), optional :: t1, t2, pick(:), delay(:)
   character(len=*), intent(in), optional :: type
   integer, intent(in), optional :: n
   real(rs) :: w1, w2
   integer :: npts
   real(rs), dimension(size(s)) :: pick_in, delay_in
   character(len=STACK_CHAR_LEN) :: type_in
   integer :: n_in

   ! Basic checks on input
   call stack_check(s)

   ! Populate internal pick and delay arrays
   if (present(pick)) then
      pick_in = pick
   else
      pick_in = 0._rs
   endif
   if (present(delay)) then
      delay_in = delay
   else
      delay_in = 0._rs
   endif

   ! Check input
   w1 = maxval(s%b)
   w2 = minval(s%e)
   if (present(t1)) w1 = t1
   if (present(t2)) w2 = t2
   if (w2 <= w1) call stack_error('stack_sum: Start time must be before end time')
   if (present(pick)) then
      if (size(pick) /= size(s)) &
         call stack_error('stack_sum: Pick array must be same length as number of traces')
      if (any(pick - delay > s%e .or. pick - delay < s%b)) &
         call stack_error('stack_sum: Picks must not be outside trace')
      if (any(pick + t2 - delay > s%e .or. pick - t1 - delay < s%b)) &
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

   ! Peform the stack
   select case (type_in(1:1))
      case ('l', 'L')
         call stack_sum_linear(s, pick_in, delay_in, w1, w2, npts, out%trace)
      case ('n', 'N')
         call stack_sum_nthroot(s, pick_in, delay_in, w1, w2, n_in, npts, out%trace)
      case ('p', 'P')
         call stack_sum_phaseweighted(s, pick_in, delay_in, w1, w2, n_in, npts, out%trace)
      case default
         call stack_error('stack_sum: Unimplemented stack type "'//trim(type_in)//'"')
   end select

end subroutine stack_sum
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_sum_linear(s, pick, delay, w1, w2, npts, out)
!===============================================================================
! Perform a linear stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     delay(:): Delay by which to shift traces (one per SAC trace)
!     w1, w2  : Start and stop times of stack, relative to pick(:)
!     npts    : Length of stacked trace
!  OUTPUT:
!     out(npts) : Array containing points of the stack
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: pick(:), delay(:), w1, w2
   integer, intent(in) :: npts
   real(rs), intent(out) :: out(npts)
   integer :: i, j, iw1, iw2

   out = 0._rs
!$omp parallel do default(none) shared(s, npts, out, pick, delay, w1, w2) &
!$omp    private(i, j, iw1, iw2)
   do i = 1, size(s)
      iw1 = nint((pick(i) - delay(i) + w1 - s(i)%b)/s(1)%delta) + 1
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
subroutine stack_sum_nthroot(s, pick, delay, w1, w2, n, npts, out)
!===============================================================================
! Perform an N-th root stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     delay(:): Delay by which to shift traces (one per SAC trace)
!     w1, w2  : Start and stop times of stack, relative to pick(:)
!     n       : Power to which to raise stack
!     npts    : Length of stacked trace
!  OUTPUT:
!     out(npts) : Array containing points of the stack
!
   type(SACtrace), intent(in) :: s(:)
   real(rs), intent(in) :: pick(:), delay(:), w1, w2
   integer, intent(in) :: n, npts
   real(rs), intent(out) :: out(npts)
   integer :: i, j, iw1, iw2

   out = 0._rs
!$omp parallel do default(none) shared(s, npts, n, out, pick, delay, w1, w2) &
!$omp    private(i, j, iw1, iw2)
   do i = 1, size(s)
      iw1 = nint((pick(i) - delay(i) + w1 - s(i)%b)/s(1)%delta) + 1
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
subroutine stack_sum_phaseweighted(s, pick, delay, w1, w2, nu, npts, out)
!===============================================================================
! Perform a phase-weighted stack.
!  INPUT:
!     s(:)    : Array of SAC traces for stacking
!     pick(:) : Times relative to which w1 and w2 are made (one per SAC trace)
!     delay(:): Delay by which to shift traces (one per SAC trace)
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
   real(rs), intent(in) :: pick(:), delay(:), w1, w2
   integer, intent(in) :: nu, npts
   real(rs), intent(out) :: out(npts)
   real(rs), dimension(npts,size(s)) :: trace, hilbert, phase
   real(rs) :: weight
   integer :: n, i, j, iw1, iw2

   n = size(s)

   ! Create Hilbert traces
   do i = 1, n
      iw1 = nint((pick(i) - delay(i) + w1 - s(i)%b)/s(1)%delta) + 1
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
function stack_station_delay(x, y, ux, uy) result(t)
!===============================================================================
! Return an array containing the relative station delay times for a wave with
! horizontal slowness components ux and uy, at stations with relative coordinates
! in x and y.
!  INPUT:
!     x(:), y(:) : Station coordinates relative to the array centre (units
!                  determine output)
!     ux, uy     : Horizontal slowness vector components: ux is east, uy is north
!                  (units should be s/distance, i.e., s/km)
!  OUTPUT:
!     t(:)       : Relative delay times for stations in s
!
   real(rs), intent(in) :: x(:), y(:), ux, uy
   real(rs) :: t(size(x))

   if (size(y) /= size(x)) &
      call stack_error('stack_station_delay: x and y locations must have same length')
   t = -(x*ux + y*uy)
end function stack_station_delay
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
subroutine stack_geog_mean(lon, lat, mean_lon, mean_lat)
!===============================================================================
! Calculate  mean geographic coordinates, given arrays of longitude and
! latitude in degrees.  Assumes points are on the same sphere.
!
   real(rs), intent(in), dimension(:) :: lon, lat
   real(rs), intent(out) :: mean_lon, mean_lat
   real(rs) :: x(3), xs(3), r
   integer :: i

   if (size(lon) /= size(lat)) &
      call stack_error('stack_geog_mean: Input arrays must be same length')
   xs = 0.
   do i = 1, size(lon)
      call stack_geog2cart(lon(i), lat(i), 1., x(1), x(2), x(3))
      xs = xs + x
   enddo
   mean_lon = atan2(xs(2), xs(1))*todeg
   r = sqrt(sum(xs**2))
   mean_lat = asin(xs(3)/r)*todeg
end subroutine stack_geog_mean
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_calc_station_coords(lon, lat, mean_lon, mean_lat, xs, ys, dist, phi)
!===============================================================================
! Calculate the cartesian coordinates of the stations, relative to the mean
! longitude and latitude.  This is done by projecting the points on a sphere
! to the tangent plane at the centre of the array (lon_mean, lat_mean).
!  INPUT:
!     lon, lat : Arrays of geographic points (degrees)
!  OUTPUT:
!     mean_lon, mean_lat : Spherical mean of input points
!     xs, ys   : Arrays of coordinates of array plane, with y parallel to local
!                north at the mean point, and y parallel to local east (km)
!     dist     : Array of distances from mean point to station (km)
!     phi      : Array of azimuths from mean point to station (deg)
   real(rs), intent(in), dimension(:) :: lon, lat
   real(rs), intent(out) :: mean_lon, mean_lat
   real(rs), intent(out), dimension(size(lon)) :: xs, ys, dist, phi
   real(rs) :: x, y, z, d
   real(rs) :: m(3), n(3), r(3), s(3), t(3), north(3), p(3), norm
   integer :: i

   if (.not.all([size(lat), size(xs), size(ys), size(dist), size(phi)] &
                == size(lon))) &
      call stack_error('stack_calc_station_coords: Arrays must be the same length')
   call stack_geog_mean(lon, lat, mean_lon, mean_lat)
   if (abs(mean_lat) >= 89.999999) &
      call stack_error('stack_calc_station_coords: Array mean is on the pole')

   ! Find mean point and normal to plane
   call stack_geog2cart(mean_lon, mean_lat, R_EARTH_KM, m(1), m(2), m(3))
   n = m/sqrt(sum(m**2))
   ! Local north unit vector at mean point
   call stack_lonlat2local_north(mean_lon, mean_lat, north(1), north(2), north(3))
!    write(*,*) 0,0,0,0
!    write(*,*) n,0
!    write(*,*)
!    write(*,*)
!    write(*,*) 0,0,0,1
!    write(*,*) north,1
!    if (abs(dot_product(north, n)) > 1.e30*tiny(1._rs)) & ! Should be orthogonal
!       call stack_error('stack_calc_station_coords: Error in calculating normal or local north')
   do i = 1, size(lon)
      ! Calculate cartesian point of station on surface in km, r
      call stack_geog2cart(lon(i), lat(i), R_EARTH_KM, r(1), r(2), r(3))
      ! Vector to point of station projected to plane, s
      d = dot_product(n, m - r)
      s = r + d*n
      ! Vector t from mean point to station has length dist km
      t = s - m
      dist(i) = sqrt(sum(t**2))
      ! Angle
      phi(i) = acos(dot_product(t, north)/dist(i))
      ! p is cross-product between local north and vector from mean point
      ! to station
      p = cross_product(t, north)/dist(i)
      ! The angle is in the opposite direction if the cross product is antiparallel
      ! to the plane normal (i.e., p ^ n < 0)
      if (dot_product(p, n) < 0.) phi(i) = -phi(i)
   enddo

   ! And rexpress things onto the plane
   xs = dist*sin(phi)
   ys = dist*cos(phi)
   phi = todeg*phi
   phi = modulo(phi, 360._rs)

end subroutine stack_calc_station_coords
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_geog2cart(lon, lat, r, x, y, z)
!===============================================================================
! Return the cartesian coordinates from geographic coordinates
   real(rs), intent(in) :: lon, lat, r
   real(rs), intent(out) :: x, y, z
   real :: lo, la

   lo = lon*torad
   la = lat*torad
   x = r*cos(lo)*cos(la)
   y = r*sin(lo)*cos(la)
   z = r*sin(la)
end subroutine stack_geog2cart
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_lonlat2local_north(lon, lat, x, y, z)
!===============================================================================
! Find the unit vector which points locally north on the unit sphere, for
! geographic coodinates in degrees
   real(rs), intent(in) :: lon, lat
   real(rs), intent(out) :: x, y, z
   real(rs) :: lo, la

   lo = lon*torad
   la = lat*torad
   x = -cos(lo)*sin(la)
   y = -sin(lo)*sin(la)
   z = cos(la)
end subroutine stack_lonlat2local_north
!-------------------------------------------------------------------------------

!===============================================================================
function stack_azimuth(lon1, lat1, lon2, lat2) result(azimuth)
!===============================================================================
! Return the azimuth from point (lon1,lat1) to (lon2,lat2).  I/O in degrees.

   real(rs), intent(in) :: lon1, lat1, lon2, lat2
   real(rs) :: azimuth, rlon1, rlat1, rlon2, rlat2

   rlon1 = torad*lon1;  rlon2 = torad*lon2
   rlat1 = torad*lat1;  rlat2 = torad*lat2

   azimuth = atan2(sin(rlon2-rlon1)*cos(rlat2) , &
                   cos(rlat1)*sin(rlat2) - sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1))
   azimuth = modulo(azimuth*todeg, 360._rs)

end function stack_azimuth

!===============================================================================
function stack_delta(lon1_in, lat1_in, lon2_in, lat2_in) result(delta)
!===============================================================================
! Returns the angular distance between two points on a sphere given
! the lat and lon of each using the Haversine formula.  I/O is in degrees.
!
   real(rs), intent(in) :: lat1_in, lon1_in, lat2_in, lon2_in
   real(rs) :: delta, lat1, lon1, lat2, lon2

   lat1 = torad*lat1_in;  lon1 = torad*lon1_in
   lat2 = torad*lat2_in;  lon2 = torad*lon2_in

   delta = atan2(sqrt((cos(lat2)*sin(lon2-lon1))**2 + (cos(lat1)*sin(lat2) - &
                       sin(lat1)*cos(lat2)*cos(lon2-lon1))**2) , &
                 sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))
   delta = delta*todeg

end function stack_delta
!-------------------------------------------------------------------------------

!===============================================================================
function cross_product(a, b) result(c)
!===============================================================================
   real(rs), intent(in) :: a(3), b(3)
   real(rs) :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function cross_product
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

!===============================================================================
subroutine stack_allocate_r4_2d(a, m, n)
!===============================================================================
   real(r4), intent(inout), allocatable, dimension(:, :) :: a
   integer, intent(in) :: m, n
   if (allocated(a)) then
      if (size(a,1) == m .and. size(a,2) == n) return
      deallocate(a)
      allocate(a(m,n))
      return
   else
      allocate(a(m,n))
   endif
end subroutine stack_allocate_r4_2d
!-------------------------------------------------------------------------------

end module stack
!-------------------------------------------------------------------------------
