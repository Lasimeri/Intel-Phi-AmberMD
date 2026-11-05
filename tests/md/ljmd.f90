! ljmd.f90: test simulation for the card's monitoring path. Liquid argon (Lennard-Jones),
! velocity Verlet, cell list, Berendsen thermostat, all forces and updates run on the thread
! team (src/card/team.c) through bind(C) callbacks. Each step it updates the MD stats block
! (src/card/mdstat.c) that phitop reads, and every few steps a frame of positions
! (src/card/mdframe.c) for the host viewer amview. It is a load generator with real physics, not Amber.
!
! Usage: ljmd [lattice cells per side (23 -> 48,668 atoms)] [steps] [threads] [print interval]
!             [frame interval in steps, 0 = no frames]
! Units: Angstrom, fs, amu, kcal/mol, K.
module ljmd_mod
  use iso_c_binding
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: eps = 0.2381_dp, sig = 3.405_dp, mass = 39.948_dp
  real(dp), parameter :: rc = 8.5_dp, kb = 0.0019872041_dp, t0 = 94.4_dp
  real(dp), parameter :: acc = 4.184e-4_dp            ! (kcal/mol/A)/amu -> A/fs^2
  integer :: n, ncell, nc3, nthr
  real(dp) :: boxl, cellw, dt, lambda
  real(dp), allocatable :: x(:,:), v(:,:), f(:,:), epart(:), kpart(:)
  integer, allocatable :: head(:), next(:)
  real(c_float), pointer :: fbuf(:,:) => null()        ! the frame slot being filled, in the shared region
  real(c_float), allocatable :: wpart(:)

  interface
    function mdframe_open(natoms) bind(C, name="mdframe_open") result(r)
      import :: c_int64_t, c_int
      integer(c_int64_t), value :: natoms
      integer(c_int) :: r
    end function mdframe_open
    function mdframe_buffer() bind(C, name="mdframe_buffer") result(p)
      import :: c_ptr
      type(c_ptr) :: p
    end function mdframe_buffer
    subroutine mdframe_publish(step, bx, by, bz, wmax) bind(C, name="mdframe_publish")
      import :: c_int64_t, c_float
      integer(c_int64_t), value :: step
      real(c_float), value :: bx, by, bz, wmax
    end subroutine mdframe_publish
    subroutine mdframe_close() bind(C, name="mdframe_close")
    end subroutine mdframe_close
    function team_init(k) bind(C, name="team_init") result(r)
      import :: c_int
      integer(c_int), value :: k
      integer(c_int) :: r
    end function team_init
    subroutine team_run(fn, ctx) bind(C, name="team_run")
      import :: c_funptr, c_ptr
      type(c_funptr), value :: fn
      type(c_ptr), value :: ctx
    end subroutine team_run
    function team_busy_ns() bind(C, name="team_busy_ns") result(r)
      import :: c_int64_t
      integer(c_int64_t) :: r
    end function team_busy_ns
    function team_ns() bind(C, name="team_ns") result(r)
      import :: c_int64_t
      integer(c_int64_t) :: r
    end function team_ns
    subroutine team_end() bind(C, name="team_end")
    end subroutine team_end
    function mdstat_open(name, natoms, nstlim, dt_mfs, threads) bind(C, name="mdstat_open") result(r)
      import :: c_char, c_int64_t, c_int
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int64_t), value :: natoms, nstlim, dt_mfs, threads
      integer(c_int) :: r
    end function mdstat_open
    subroutine mdstat_step(nstep, etot, temp, busy, stepns) bind(C, name="mdstat_step")
      import :: c_int64_t
      integer(c_int64_t), value :: nstep, etot, temp, busy, stepns
    end subroutine mdstat_step
    function mdstat_stop_requested() bind(C, name="mdstat_stop_requested") result(r)
      import :: c_int
      integer(c_int) :: r
    end function mdstat_stop_requested
    subroutine mdstat_close(state) bind(C, name="mdstat_close")
      import :: c_int32_t
      integer(c_int32_t), value :: state
    end subroutine mdstat_close
  end interface

contains

  pure subroutine span(total, tid, k, lo, hi)
    integer, intent(in) :: total, tid, k
    integer, intent(out) :: lo, hi
    lo = int((int(total, 8) * tid) / k) + 1
    hi = int((int(total, 8) * (tid + 1)) / k)
  end subroutine span

  ! Forces on the atoms of this thread's cells (full neighbour loop, no write to other threads' atoms).
  subroutine job_force(ctx, tid, k) bind(C)
    type(c_ptr), value :: ctx
    integer(c_int), value :: tid, k
    integer :: lo, hi, c, cx, cy, cz, dx, dy, dz, nb, i, j
    real(dp) :: xi(3), d(3), r2, ir2, ir6, fr, fi(3), e, rc2, s6
    if (c_associated(ctx)) continue
    rc2 = rc * rc
    s6 = sig**6
    e = 0
    call span(nc3, tid, k, lo, hi)
    do c = lo, hi
      cz = (c - 1) / (ncell * ncell)
      cy = mod((c - 1) / ncell, ncell)
      cx = mod(c - 1, ncell)
      i = head(c)
      do while (i > 0)
        xi = x(:, i)
        fi = 0
        do dz = -1, 1
          do dy = -1, 1
            do dx = -1, 1
              nb = 1 + modulo(cx + dx, ncell) + ncell * (modulo(cy + dy, ncell) + ncell * modulo(cz + dz, ncell))
              j = head(nb)
              do while (j > 0)
                if (j /= i) then
                  d = x(:, j) - xi
                  d = d - boxl * anint(d / boxl)
                  r2 = d(1)**2 + d(2)**2 + d(3)**2
                  if (r2 < rc2) then
                    ir2 = 1 / r2
                    ir6 = s6 * ir2 * ir2 * ir2
                    fr = 24 * eps * ir6 * (2 * ir6 - 1) * ir2
                    fi = fi - fr * d
                    e = e + 2 * eps * ir6 * (ir6 - 1)
                  end if
                end if
                j = next(j)
              end do
            end do
          end do
        end do
        f(:, i) = fi
        i = next(i)
      end do
    end do
    epart(tid) = e
  end subroutine job_force

  ! First half kick and drift, with the thermostat factor of the previous step.
  subroutine job_drift(ctx, tid, k) bind(C)
    type(c_ptr), value :: ctx
    integer(c_int), value :: tid, k
    integer :: lo, hi, i
    if (c_associated(ctx)) continue
    call span(n, tid, k, lo, hi)
    do i = lo, hi
      v(:, i) = lambda * v(:, i) + 0.5_dp * dt * acc / mass * f(:, i)
      x(:, i) = modulo(x(:, i) + dt * v(:, i), boxl)
    end do
  end subroutine job_drift

  ! Second half kick and the kinetic energy of this thread's atoms.
  subroutine job_kick(ctx, tid, k) bind(C)
    type(c_ptr), value :: ctx
    integer(c_int), value :: tid, k
    integer :: lo, hi, i
    real(dp) :: s
    s = 0
    if (c_associated(ctx)) continue
    call span(n, tid, k, lo, hi)
    do i = lo, hi
      v(:, i) = v(:, i) + 0.5_dp * dt * acc / mass * f(:, i)
      s = s + v(1, i)**2 + v(2, i)**2 + v(3, i)**2
    end do
    kpart(tid) = 0.5_dp * mass * s / acc
  end subroutine job_kick

  ! This thread's atoms into the frame slot: x, y, z in Angstrom and the speed in A/ps as w.
  subroutine job_frame(ctx, tid, k) bind(C)
    type(c_ptr), value :: ctx
    integer(c_int), value :: tid, k
    integer :: lo, hi, i
    real(c_float) :: sp, wmax
    wmax = 0
    if (c_associated(ctx)) continue
    call span(n, tid, k, lo, hi)
    do i = lo, hi
      sp = real(1000 * sqrt(v(1, i)**2 + v(2, i)**2 + v(3, i)**2), c_float)
      fbuf(1, i) = real(x(1, i), c_float)
      fbuf(2, i) = real(x(2, i), c_float)
      fbuf(3, i) = real(x(3, i), c_float)
      fbuf(4, i) = sp
      wmax = max(wmax, sp)
    end do
    wpart(tid) = wmax
  end subroutine job_frame

  ! Cell list (serial): head of each cell, next atom in the same cell.
  subroutine build_cells()
    integer :: i, c, ix(3)
    head = 0
    do i = n, 1, -1
      ix = min(int(x(:, i) / cellw), ncell - 1)
      c = 1 + ix(1) + ncell * (ix(2) + ncell * ix(3))
      next(i) = head(c)
      head(c) = i
    end do
  end subroutine build_cells

end module ljmd_mod

program ljmd
  use ljmd_mod
  implicit none
  integer :: m, nsteps, nprint, nframe, nframes, step, i, a, b, c, q, ios
  character(len=32) :: arg
  real(dp) :: lat, ek, ep, temp, u1, u2
  integer(c_int64_t) :: tstart, t1, tstep0, stepns, tf0, framens
  logical :: frames
  real(dp), parameter :: basis(3,4) = reshape([0d0,0d0,0d0, .5d0,.5d0,0d0, .5d0,0d0,.5d0, 0d0,.5d0,.5d0], [3,4])
  integer(8) :: seed
  logical :: stopped

  m = 23; nsteps = 2000; nthr = 120; nprint = 100; nframe = 1
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*,iostat=ios) m; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*,iostat=ios) nsteps; end if
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read(arg,*,iostat=ios) nthr; end if
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read(arg,*,iostat=ios) nprint; end if
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read(arg,*,iostat=ios) nframe; end if

  n = 4 * m**3
  lat = (4 / 0.0211_dp)**(1.0_dp / 3)          ! liquid argon, 0.0211 atoms per A^3
  boxl = m * lat
  ncell = max(3, int(boxl / rc))
  cellw = boxl / ncell
  nc3 = ncell**3
  dt = 4.0_dp
  lambda = 1
  allocate(x(3, n), v(3, n), f(3, n), head(nc3), next(n))
  allocate(epart(0:nthr-1), kpart(0:nthr-1), wpart(0:nthr-1))

  q = 0
  do c = 0, m - 1; do b = 0, m - 1; do a = 0, m - 1; do i = 1, 4
    q = q + 1
    x(:, q) = lat * ([real(a, dp), real(b, dp), real(c, dp)] + basis(:, i))
  end do; end do; end do; end do
  seed = 88172645463325252_8
  do i = 1, n                                 ! Box-Muller from an xorshift stream
    do a = 1, 3
      seed = ieor(seed, ishft(seed, 13)); seed = ieor(seed, ishft(seed, -7)); seed = ieor(seed, ishft(seed, 17))
      u1 = (real(ibits(seed, 11, 52), dp) + 0.5_dp) * 2.0_dp**(-52)
      seed = ieor(seed, ishft(seed, 13)); seed = ieor(seed, ishft(seed, -7)); seed = ieor(seed, ishft(seed, 17))
      u2 = real(ibits(seed, 11, 52), dp) * 2.0_dp**(-52)
      v(a, i) = sqrt(-2 * log(u1)) * cos(8 * atan(1.0_dp) * u2) * sqrt(kb * t0 * acc / mass)
    end do
  end do
  do a = 1, 3
    v(a, :) = v(a, :) - sum(v(a, :)) / n
  end do

  if (team_init(int(nthr, c_int)) /= nthr) stop 'team_init failed'
  if (mdstat_open('ljmd argon' // c_null_char, int(n, c_int64_t), int(nsteps, c_int64_t), &
                  int(dt * 1000, c_int64_t), int(nthr, c_int64_t)) /= 0) &
    write(6, '(a)') 'mdstat: no stats block (not on the card?), running anyway'
  frames = nframe > 0
  if (frames) frames = mdframe_open(int(n, c_int64_t)) == 0
  if (nframe > 0 .and. .not. frames) write(6, '(a)') 'mdframe: no frame ring, running without frames'
  nframes = 0
  framens = 0
  write(6, '(a,i0,a,i0,a,f8.3,a,i0,a,i0,a)') 'ljmd: ', n, ' atoms, ', nthr, ' threads, box ', boxl, &
       ' A, ', ncell, '^3 cells, ', nsteps, ' steps of 4 fs'

  call build_cells()
  call team_run(c_funloc(job_force), c_null_ptr)
  tstart = team_ns()
  stepns = 0
  stopped = .false.
  do step = 1, nsteps
    tstep0 = team_ns()
    call team_run(c_funloc(job_drift), c_null_ptr)
    call build_cells()
    call team_run(c_funloc(job_force), c_null_ptr)
    call team_run(c_funloc(job_kick), c_null_ptr)
    ep = sum(epart); ek = sum(kpart)
    temp = 2 * ek / (3 * (n - 1) * kb)
    lambda = sqrt(1 + (dt / 100.0_dp) * (t0 / temp - 1))   ! Berendsen, tau 100 fs
    t1 = team_ns()
    stepns = stepns + (t1 - tstep0)
    call mdstat_step(int(step, c_int64_t), nint((ep + ek) * 1000, c_int64_t), nint(temp * 1000, c_int64_t), &
                     team_busy_ns(), stepns)
    if (frames) then
      if (mod(step, nframe) == 0) then
        tf0 = team_ns()
        call c_f_pointer(mdframe_buffer(), fbuf, [4, n])
        call team_run(c_funloc(job_frame), c_null_ptr)
        call mdframe_publish(int(step, c_int64_t), real(boxl, c_float), real(boxl, c_float), real(boxl, c_float), &
                             maxval(wpart))
        framens = framens + (team_ns() - tf0)
        nframes = nframes + 1
      end if
    end if
    if (mod(step, nprint) == 0 .or. step == nsteps) &
      write(6, '(a,i7,a,f14.3,a,f14.3,a,f8.3,a,f9.3,a,f6.1,a)') 'step', step, '  Epot', ep, '  Etot', ep + ek, &
        '  T', temp, '  steps/s', step / (real(t1 - tstart, dp) * 1e-9_dp), '  team busy', &
        100 * real(team_busy_ns(), dp) / (real(t1 - tstart, dp) * nthr), '%'
    if (mdstat_stop_requested() /= 0) then
      stopped = .true.
      exit
    end if
  end do
  if (stopped) then
    write(6, '(a,i0)') 'ljmd: stopped by the host at step ', step
    call mdstat_close(4_c_int32_t)
  else
    call mdstat_close(2_c_int32_t)
  end if
  if (frames) then
    call mdframe_close()
    write(6, '(a,i0,a,f10.3,a)') 'ljmd: ', nframes, ' frames published, ', &
      real(framens, dp) * 1e-6_dp / max(nframes, 1), ' ms per frame (fill + publish)'
  end if
  call team_end()
  write(6, '(a,f10.3,a)') 'ljmd: done, ', real(team_ns() - tstart, dp) * 1e-9_dp, ' s'
end program ljmd
