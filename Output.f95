subroutine OutputFields(IO, NI, NJ, x, y, U, V, P, R)
implicit none
integer :: NI, NJ, IO
real, dimension(0 : NI + 1, 0 : NJ + 1) :: X, Y
real, dimension(0 : NI, 0 : NJ) :: U, V, P, R

write(IO, *) &
'Variables = "X", "Y", "U", "V", "P", "Rho"'
write(IO, *) 'Zone i =', NI + 2, ', j =', NJ + 2, &
',DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
write(IO,'(100e25.16)') X(0 : NI + 1, 0 : NJ + 1)
write(IO,'(100e25.16)') Y(0 : NI + 1, 0 : NJ + 1)
write(IO,'(100e25.16)') U(0 : NI, 0 : NJ)
write(IO,'(100e25.16)') V(0 : NI, 0 : NJ)
write(IO,'(100e25.16)') P(0 : NI, 0 : NJ)
write(IO,'(100e25.16)') R(0 : NI, 0 : NJ)

end subroutine


subroutine Output_Cf(IO, NI, X, Cf, U_max)
implicit none
integer :: IO, NI
real, dimension(0 : NI + 1) :: X
real, dimension(0 : NI) :: Cf, U_max

write(IO, *) 'Variables = "X", "Cf", "U_max"'
write(IO, *) 'Zone i =', NI + 2, &
',DATAPACKING=BLOCK, VARLOCATION=([2-20]=CELLCENTERED)'
write(IO,'(100E25.16)') X(0 : NI + 1)
write(IO,'(100E25.16)') Cf(0 : NI)
write(IO,'(100E25.16)') U_max(0 : NI)

end subroutine
