subroutine Friction_coefficient(NI, NJ, dy, Cf, U, MU, R0, U0)
implicit none
integer :: NI, NJ, i
real :: dy, mu, r0, u0
real, dimension(0 : NI, 0 : NJ) :: U
real, dimension(0 : NI) :: Cf

do i = 0, NI
Cf(i) = (2 * mu * (U(i, 1) - U(i, 0)) / (dy)) / (r0 * u0**2)
end do

end subroutine
