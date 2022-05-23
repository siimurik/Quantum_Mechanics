program fin_pot_well
  implicit none

  real(8) :: h, hbar, U0, a, m
  real(8) :: E, Evana, limit, f, df, vahed
  integer(8) :: count
  real(8), parameter :: pi = 4.0*atan(1.0)

  E = 1.312E-40
  Evana = 0.0
  vahed = abs(E-Evana)
  count = 0

  h   = 6.626068e-34 ! m^2 kg / s
  hbar = h/(2*pi)
  E   = 5*1.60217663E-19
  U0  = 10*1.60217663E-19
  m   = 1.67262192E-27
  a   = 0.5

  limit = 1e-44

  do while ( abs(E-Evana) .gt. limit )
    print*, 'E_', count,'=', E, '|E - Evana| =', vahed
    count = count + 1
    Evana = E
    f  = tan(a*sqrt(2*m*E)/hbar) - 2*sqrt(E*(U0 - E))/(2*E - U0)
    df = -((-2*E + U0)/(sqrt(E)*(2*E - U0)*sqrt(-E + U0))) + &
          (4*sqrt(E)*sqrt(-E + U0))/(2*E - U0)**2 + &
          a*sqrt(m)*1/cos((a*sqrt(2*E)*sqrt(m))/hbar**2)/(sqrt(2*E)*hbar)
    E  = E - f/df
    vahed = abs(E-Evana)
  end do

end program fin_pot_well
