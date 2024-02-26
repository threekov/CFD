subroutine NavierStokes(NI, NJ, dx, dy, U_c, V_c, P_c, R_c, &
CFL, u_ref, MAXITER, eps, u0, mu, r0, p0, gamma, k, X_cell)
implicit none
integer :: NI, NJ, MAXITER
real :: dx, dy, CFL, u_ref, eps, u0, mu, r0, p0, gamma, k
real, dimension(0 : NI, 0 : NJ) :: X_cell, U_c, V_c, P_c, R_c
integer :: i, j, counter
real :: A, dt, u12, v12, p12, r_u_c12, r_v_c12
real :: Mass, Imp_x, Imp_y, Imp_x2, Imp_y2
real, dimension(0 : NI, 0 : NJ) :: Res_P, Res_U, Res_V

A = 1.0 / (u_ref * u_ref)

counter = 0

if (dx .gt. dy) then
	dt = CFL * dy / u0
else
	dt = CFL * dx / u0
end if

open(1, file = 'Residuals.plt')
write(1, *) "Variables = Iteration, Residual P, Residual U, Residual V"
write(1, *) 'Zone i =', MAXITER

do while ((counter .eq. 0) .or. (((maxval(abs(Res_U)) .gt. eps) .or. &
(maxval(abs(Res_V)) .gt. eps) .or. (maxval(abs(Res_P)) .gt. eps)) &
.and. (counter .lt. MAXITER)))

	Res_U(:,:) = 0.0
	Res_V(:,:) = 0.0
	Res_P(:,:) = 0.0

	! inlet boundary condition
	do j = 1, NJ - 1
		U_c(0, j) = 2 * u0 - U_c(1, j) !u0
		V_c(0, j) = - V_c(1, j) !0.0
		P_c(0, j) = P_c(1, j)
	end do

	! outlet boundary condition
	do j = 1, NJ - 1
		U_c(NI, j) = U_c(NI - 1, j)
		V_c(NI, j) = V_c(NI - 1, j)
		P_c(NI, j) = 2 * p0 - P_c(NI - 1, j) !p0
	end do

	! wall boundary condition
	do i = 1, NI - 1
		U_c(i, 0) = - U_c(i, 1)
		V_c(i, 0) = - V_c(i, 1) + 2.0 * k * X_cell(i, 0) * u0
		P_c(i, 0) = P_c(i, 1)
	end do


	! external boundary condition
	do i = 1, NI - 1
		V_c(i, NJ) = V_c(i, NJ - 1)
		if (V_c(i, NJ - 1) .gt. 0.0) then
			U_c(i, NJ) = U_c(i, NJ - 1)
			P_c(i, NJ) = 2 * p0 - P_c(i, NJ - 1)
		else
			U_c(i, NJ) = 2 * u0 - U_c(i, NJ - 1)
			P_c(i, NJ) = P_c(i, NJ - 1)
		end if
	end do

	! edge flux calculation (x)
	do j = 1, NJ - 1
		do i = 0, NI - 1

			r_u_c12 = (R_c(i + 1, j) * U_c(i + 1, j) + R_c(i, j) * U_c(i, j)) / 2.0

			if (r_u_c12 .ge. 0.0) then
				u12 = U_c(i, j)
				v12 = V_c(i, j)
				p12 = P_c(i + 1, j)
			else
				u12 = U_c(i + 1, j)
				v12 = V_c(i + 1, j)
				p12 = P_c(i, j)
			end if

			Mass = (R_c(i, j) + R_c(i + 1, j)) / 2.0 * u12
			Imp_x = r_u_c12 * u12 + p12 - 4.0 / 3.0 * mu * (U_c(i + 1, j) - U_c(i, j)) / dx
			Imp_y = r_u_c12 * v12 - mu * (V_c(i + 1, j) - V_c(i, j)) / dx
			Imp_x2 = 1.0 / 3.0 * mu * (V_c(i + 1, j + 1) - V_c(i + 1, j - 1)) / (2.0 * dy)
			Imp_y2 = - mu * (U_c(i + 1, j + 1) - U_c(i + 1, j - 1)) / (4.0 * dy)

			if (i .gt. 0) then
				Imp_x2 = Imp_x2 - 1.0 / 3.0 * mu * (V_c(i - 1, j + 1) - V_c(i - 1,j - 1)) / (2.0 * dy)
				Imp_y2 = Imp_y2 + mu * (U_c(i - 1,j + 1) - U_c(i - 1,j - 1)) / (4.0 * dy)
			end if

			! i, j
			Res_P(i, j) = Res_P(i, j) + Mass / dx
			Res_U(i, j) = Res_U(i, j) + (Imp_x + Imp_x2) / dx
			Res_V(i, j) = Res_V(i, j) + (Imp_y + Imp_y2) / dx

			! i + 1, j
			Res_P(i + 1, j) = Res_P(i + 1, j) - Mass / dx
			Res_u(i + 1, j) = Res_U(i + 1, j) - Imp_x / dx
			Res_v(i + 1, j) = Res_V(i + 1, j) - Imp_y / dx

		end do

		Res_P(0, j) = 0.0
		Res_U(0, j) = 0.0
		Res_V(0, j) = 0.0
		Res_P(NI, j) = 0.0
		Res_U(NI, j) = 0.0
		Res_V(NI, j) = 0.0

	end do

	! edge flux calculation (y)
	do i = 1, NI - 1
		do j = 0, NJ - 1

			r_v_c12 = (R_c(i, j + 1) * V_c(i, j + 1) + R_c(i, j) * V_c(i, j)) / 2.0

			if (r_v_c12 .ge. 0.0) then
				u12 = U_c(i, j)
				v12 = V_c(i, j)
				p12 = P_c(i, j + 1)
			else
				u12 = U_c(i, j + 1)
				v12 = V_c(i, j + 1)
				p12 = P_c(i, j)
			end if

			Mass = (R_c(i, j) + R_c(i, j + 1)) / 2.0 * v12
			Imp_x = r_v_c12 * u12 - mu * (U_c(i, j + 1) - U_c(i, j)) / dy
			Imp_y = r_v_c12 * v12 + p12 - 4.0/3.0 * mu * (V_c(i, j + 1) - V_c(i, j)) / dy
			Imp_x2 = - mu * (V_c(i + 1, j + 1) - V_c(i - 1, j + 1)) / (4.0 * dx)
			Imp_y2 = 2.0 / 3.0 * MU * (U_c(i + 1, j + 1) - U_c(i - 1, j + 1)) / (4.0 * dx)

			if (j .gt. 0) then
				Imp_x2 = Imp_x2 + mu * (V_c(i + 1, j - 1) - V_c(i - 1, j - 1)) / (4.0 * dx)
				Imp_y2 = Imp_x2 - 2.0 / 3.0 * mu * (U_c(i + 1, j - 1) - U_c(i - 1, j - 1)) / (4.0 * dx)
			end if

			! correction at the required boundaries
			if (j .eq. 0) then
				Mass = (R_c(i, j) + R_c(i, j + 1)) / 2.0 * r_v_c12
			end if

			! i, j
			Res_P(i, j) = Res_P(i, j) + Mass / dy
			Res_U(i, j) = Res_U(i, j) + (Imp_x + Imp_x2) / dy
			Res_V(i, j) = Res_V(i, j) + (Imp_y + Imp_y2) / dy

			! i, j + 1
			Res_P(i, j + 1) = Res_p(i, j + 1) - Mass / dy
			Res_U(i, j + 1) = Res_u(i, j + 1) - Imp_x / dy
			Res_V(i, j + 1) = Res_v(i, j + 1) - Imp_y / dy

		end do

		Res_P(i, 0) = 0.0
		Res_U(i, 0) = 0.0
		Res_V(i, 0) = 0.0
		Res_P(i, NJ) = 0.0
		Res_U(i, NJ) = 0.0
		Res_V(i, NJ) = 0.0

	end do

	do i = 0, NI
		do j = 0, NJ
			U_c(i, j) = U_c(i, j) - dt * Res_U(i, j) / R_c(i, j)
			V_c(i, j) = V_c(i, j) - dt * Res_V(i, j) / R_c(i, j)
			P_c(i, j) = P_c(i, j) - dt / A * Res_P(i, j)
		end do
	end do
	
	! cycle progress
	if (mod(counter + 1, 1000) .eq. 0) then
		print *, 'Iteration =', counter + 1
	end if

	do i = 0, NI
		do j = 0, NJ
			R_c(i,j) = R0 * (P_c(i, j) / p0)**(1.0 / gamma)
		end do
	end do

	counter = counter + 1

	write(1, *) counter, maxval(abs(Res_P)) / A, maxval(abs(Res_U)), maxval(abs(Res_V))

end do

close(1)

end subroutine