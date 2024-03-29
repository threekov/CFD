program var5
implicit none
integer, parameter :: IO = 12 ! input-output unit
integer ::  i, j, NI, NJ, NITER
real :: L, H, u0, mu, nu, r0, p0, gamma, k
real :: dx, dy, CFL, u_ref, eps
real, allocatable :: X_node(:,:), Y_node(:,:)
real, allocatable :: X_cell(:,:), Y_cell(:,:)
real, allocatable :: R_c(:,:), P_c(:,:), U_c(:,:), V_c(:,:)
real, allocatable :: Cf(:), U_max(:)

write(*, *) 'Read input file'
open(IO, file = 'Input.txt')
read(IO, *) L
read(IO, *) H
read(IO, *) NI
read(IO, *) NJ
read(IO, *) NITER
read(IO, *) CFL
read(IO, *) u_ref
read(IO, *) eps
read(IO, *) u0
read(IO, *) mu
read(IO, *) r0
read(IO, *) p0
read(IO, *) gamma
read(IO, *) k
close(IO)

allocate(X_node(0 : NI + 1, 0 : NJ + 1)) ! mesh nodes x-coordinates
allocate(Y_node(0 : NI + 1, 0 : NJ + 1)) ! mesh nodes y-coordinates
allocate(X_cell(0 : NI, 0 : NJ)) ! cell centers x-coordinates
allocate(Y_cell(0 : NI, 0 : NJ)) ! cell centers y-coordinates
 
!—— Cell-centered variables
allocate(U_c(0 : NI, 0 : NJ)) ! velocity u
allocate(V_c(0 : NI, 0 : NJ)) ! velocity v
allocate(P_c(0 : NI, 0 : NJ)) ! pressure
allocate(R_c(0 : NI, 0 : NJ)) ! density
allocate(Cf(0 : NI)) ! friction coefficient
allocate(U_max(0 : NI)) ! max value of velocity u
dx = L / (NI - 1)
dy = H / (NJ - 1)
 
!—— Coordinate of nodes
do i = 0, NI + 1
  do j = 0, NJ + 1
    X_node(i, j) = (i - 1) * dx
    Y_node(i, j) = (j - 1) * dy
  end do
end do

!—— Coordinate of cell centers
do i = 0, NI
  do j = 0, NJ
    X_cell(i, j) = X_node(i, j) + dx / 2
    Y_cell(i, j) = Y_node(i, j) + dy / 2
  end do
end do

!—————– Parameters ————————
nu = mu / r0 ! kinematic viscosity
write(*, *) 'L = ', L, 'NI = ', NI, 'dx = ', dx
write(*, *) 'H = ', H, 'NJ = ', NJ, 'dy = ', dy
write(*,*) 'Re = ', u0 * H / nu

!—————– Initial fields —————————–
do i = 0, NI
  do j = 0, NJ
    U_c(i, j) = u0
    V_c(i, j) = 1.0e-5
    P_c(i, j) = p0
    R_c(i,j) = r0
  end do
end do

!—————- Solve Navier-Stokes equations ———————
write(*, *) 'Solve Navier-Stokes equations'
call NavierStokes(NI, NJ, dx, dy, U_c, V_c, P_c, R_c, &
CFL, u_ref, NITER, eps, u0, mu, r0, p0, gamma, k, X_cell, U_max)

call Friction_coefficient(NI, NJ, dy, Cf, U_c, mu, r0, u0)

!—————– Output data ——————————
write(*, *) 'Output data cell (Navier-Stokes)'
open(IO, file = 'data.plt')
call OutputFields(IO, NI, NJ, X_node, Y_node, &
U_c, V_c, P_c, R_c)
close(IO)

open(IO, file = 'Cf.plt')
call Output_Cf(IO, NI, X_node(:,1), Cf, U_max)
close(IO)

end program
