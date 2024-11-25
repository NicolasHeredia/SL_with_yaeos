program calc_presion
   use yaeos__constants, only: pr, R
   use pressure__SL, only: pressure_SL  ! Importa la subrutina
   implicit none

   ! Declaración de variables
   real(pr) :: T, V, P_SL, T_SL, rho_SL, r_SL, P
   integer :: i

   ! Inicialización de variables (puedes ajustar estos valores)
   T = 175.0_pr          ! Temperatura en K
   V = 0.1_pr            ! Volumen molar en L/mol
   P_SL = 2480.0_pr      ! Presión característica en bar
   T_SL = 224.0_pr       ! Temperatura característica en K
   rho_SL = 31.17_pr     ! Densidad característica en mol/L
   r_SL = 4.26_pr        ! Segment

   ! Llamada a la subrutina
   call pressure_SL(T, V, P_SL, T_SL, rho_SL, r_SL, P)

   ! Mostrar el resultado
   print *, "Presión calculada: ", P

   ! Bucle para calcular presión
   ! Abrir archivo directamente
   open(unit=10, file="presion.txt", status="replace")

   do i=1,1000
      V = real(i, pr)/1000
      call pressure_SL(T, V, P_SL, T_SL, rho_SL, r_SL, P)
      write(10, '(E20.12, 2X, E20.12)') V, P
    end do

   ! Cerrar el archivo
   close(10)

end program calc_presion

