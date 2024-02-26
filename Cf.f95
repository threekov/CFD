subroutine Cf(NI, NJ, dy, C_f, U, mu, r0, u0)
implicit none
integer :: NI, NJ, i
real :: dy, mu, r0, u0
real, dimension(0 : NI, 0 : NJ) :: U
real, dimension(0 : NI) :: Cf

do i = 0, NI
	Cf(i) = -(mu * (3 * U(i, 1) - 4 * U(i, 2) + U(i, 3)) / (2 * dy)) / (0.5 * r0 * u0**2)
end do

end subroutine