! conf.f90: Fortran features sander depends on, checked on the card (G1 gate).
! Each check prints "PASS name" or "FAIL name detail". Values whose bits matter are also
! printed as "HEX name z" so the card can be compared with a host build bit for bit.
! The C side is conf_c.c (compiled by the KNC clang for the card, by gcc for the host).
module conf_mod
  use iso_c_binding
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  type :: atom_t
    real(dp) :: x(3)
    integer :: id
    character(len=4) :: name
  end type atom_t
  real(dp), allocatable :: frc(:,:)
  integer :: nfail = 0

  interface
    function c_cmul(a, b) bind(C, name="c_cmul") result(r)
      import :: c_double_complex
      complex(c_double_complex), value :: a, b
      complex(c_double_complex) :: r
    end function c_cmul
    function c_dot3(p, q) bind(C, name="c_dot3") result(r)
      import :: c_double
      real(c_double), intent(in) :: p(3), q(3)
      real(c_double) :: r
    end function c_dot3
    subroutine c_calls_fortran(r) bind(C, name="c_calls_fortran")
      import :: c_double_complex
      complex(c_double_complex), intent(out) :: r
    end subroutine c_calls_fortran
  end interface

contains

  subroutine check(ok, name, detail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name, detail
    if (ok) then
      write(6,'(a,1x,a)') 'PASS', name
    else
      write(6,'(a,1x,a,1x,a)') 'FAIL', name, detail
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine hex(name, v)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: v
    write(6,'(a,1x,a,1x,z16.16)') 'HEX', name, v
  end subroutine hex

  ! Called from C: the C side passes two complex numbers back through this function.
  function f_cadd(a, b) bind(C, name="f_cadd") result(r)
    complex(c_double_complex), value :: a, b
    complex(c_double_complex) :: r
    r = a + b
  end function f_cadd

end module conf_mod

program conf
  use conf_mod
  implicit none
  integer :: imin, nstlim, ntc, ios, i, j, u
  real(dp) :: dt, cut, temp0
  character(len=64) :: restraintmask
  namelist /cntrl/ imin, nstlim, dt, cut, ntc, temp0, restraintmask
  character(len=256) :: line
  real(dp) :: a, b, s, p(3), q(3)
  complex(dp) :: z, w
  type(atom_t) :: at(2)
  integer(8) :: big
  real(dp), allocatable :: buf(:)

  ! 1. NAMELIST input, the way sander reads mdin
  open(newunit=u, file='conf.mdin', status='replace', action='write')
  write(u,'(a)') 'test input'
  write(u,'(a)') ' &cntrl'
  write(u,'(a)') '   imin=0, nstlim=1000, dt=0.002,'
  write(u,'(a)') '   cut=8.0, ntc=2, temp0=310.0,'
  write(u,'(a)') "   restraintmask='!:WAT & !:Na+',"
  write(u,'(a)') ' /'
  close(u)
  imin = -1; nstlim = -1; dt = -1; cut = -1; ntc = -1; temp0 = -1; restraintmask = ''
  open(newunit=u, file='conf.mdin', status='old', action='read')
  read(u,'(a)') line
  read(u, nml=cntrl, iostat=ios)
  close(u)
  call check(ios == 0 .and. imin == 0 .and. nstlim == 1000 .and. abs(dt - 0.002_dp) < 1e-15_dp &
       .and. cut == 8.0_dp .and. ntc == 2 .and. temp0 == 310.0_dp &
       .and. trim(restraintmask) == '!:WAT & !:Na+', 'namelist', 'mdin values differ')

  ! 2. Formatted output and read-back, the mdout and restart formats
  write(line,'(f12.7,e15.7,es20.12,i8)') 3.25_dp, -1.5e-7_dp, 6.02214076e23_dp, 42
  call check(trim(line) == '   3.2500000 -0.1500000E-06  6.022140760000E+23      42', 'format-write', trim(line))
  read(line,*) a, b, s, i
  call check(a == 3.25_dp .and. abs(b + 1.5e-7_dp) < 1e-20_dp .and. abs(s/6.02214076e23_dp - 1) < 1e-12_dp &
       .and. i == 42, 'format-read', trim(line))
  write(line,'(6f12.7)') (0.1_dp * j, j = 1, 6)
  call check(line(1:12) == '   0.1000000' .and. line(61:72) == '   0.6000000', 'inpcrd-format', trim(line))

  ! 3. Unformatted stream round trip
  allocate(buf(1000))
  buf = [(sqrt(real(j, dp)), j = 1, 1000)]
  open(newunit=u, file='conf.bin', access='stream', form='unformatted', status='replace')
  write(u) buf
  close(u)
  buf = 0
  open(newunit=u, file='conf.bin', access='stream', form='unformatted', status='old')
  read(u) buf
  close(u)
  call check(all(buf == [(sqrt(real(j, dp)), j = 1, 1000)]), 'stream-io', 'values differ')

  ! 4. Allocatable module arrays, derived types, array syntax
  allocate(frc(3, 49910))
  frc = 0
  do j = 1, 49910
    frc(:, j) = [real(j, dp), -real(j, dp), 0.5_dp]
  end do
  call check(abs(sum(frc(1,:)) - 49910.0_dp*49911.0_dp/2) == 0 .and. sum(frc(3,:)) == 24955.0_dp, &
       'allocatable', 'sum differs')
  at(1) = atom_t([1.0_dp, 2.0_dp, 3.0_dp], 7, 'CA  ')
  at(2) = at(1); at(2)%id = 8
  call check(at(2)%name == 'CA  ' .and. at(2)%id == 8 .and. at(1)%id == 7 .and. at(2)%x(3) == 3, 'derived-type', '')

  ! 5. COMPLEX and arrays across bind(C), both directions
  z = c_cmul((1.5_dp, -2.0_dp), (0.25_dp, 4.0_dp))
  call check(z == (8.375_dp, 5.5_dp), 'complex-to-c', '')
  call c_calls_fortran(w)
  call check(w == (4.0_dp, -1.0_dp), 'complex-from-c', '')
  p = [1.0_dp, 2.0_dp, 3.0_dp]; q = [4.0_dp, -5.0_dp, 6.0_dp]
  call check(c_dot3(p, q) == 12.0_dp, 'array-to-c', '')

  ! 6. Integer kinds and intrinsics used in sander
  big = 2_8**40 + 3
  call check(big == 1099511627779_8 .and. mod(big, 1000_8) == 779, 'integer8', '')
  call check(index('!:WAT & !:Na+', 'Na') == 11 .and. len_trim('ab  ') == 2 .and. adjustl('  x') == 'x  ', &
       'strings', '')

  ! 7. Math library: printed in hex for host and card comparison, plus sanity bounds
  call hex('sqrt2', sqrt(2.0_dp))
  call hex('exp1', exp(1.0_dp))
  call hex('log10', log(10.0_dp))
  call hex('erfc05', erfc(0.5_dp))
  call hex('erfc23', erfc(2.3_dp))
  call hex('sin1', sin(1.0_dp))
  call hex('cos1', cos(1.0_dp))
  call hex('atan2', atan2(1.0_dp, -2.0_dp))
  call hex('pow', 1.7_dp**2.3_dp)
  call hex('dsum', sum([(1.0_dp / j, j = 1, 100000)]))
  s = 0
  do j = 1, 100000
    s = s + 1.0_dp / real(j, dp)
  end do
  call hex('loopsum', s)
  call check(abs(erfc(0.5_dp) - 0.4795001221869535_dp) < 1e-15_dp .and. abs(sqrt(2.0_dp)**2 - 2) < 1e-15_dp, &
       'math', '')

  if (nfail == 0) then
    write(6,'(a)') 'RESULT all passed'
  else
    write(6,'(a,i3,a)') 'RESULT ', nfail, ' failed'
    stop 1
  end if
end program conf
