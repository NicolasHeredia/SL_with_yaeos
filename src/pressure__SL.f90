module pressure__SL
   use yaeos__constants, only: pr, R
   implicit none
contains
   subroutine pressure_SL(T, V, P_SL, T_SL, rho_SL, r_SL, P)
       ! Argumentos 
       real(pr), intent(in) :: T
       real(pr), intent(in) :: V

       ! Parametros de Sanchez-Lacombe
       real(pr), intent(in) :: P_SL
       real(pr), intent(in) :: T_SL
       real(pr), intent(in) :: rho_SL
       real(pr), intent(in) :: r_SL

       ! Parametros internos
       real(pr) :: T_red ! reducida
       real(pr) :: rho_red
       real(pr) :: P_red

       ! Presión (salida)
       real(pr), intent(out) :: P
   

       ! Cuerpo de la subrutina
       T_red = T / T_SL
       rho_red = ((1 / V) / rho_SL)
       
       P_red = - rho_red**2 - T_red * (log(1 - rho_red) + (1 - 1 / r_SL) * rho_red)

       P = P_red * P_SL

   end subroutine pressure_SL
end module pressure__SL
