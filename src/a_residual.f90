module a_residual
   use yaeos__constants, only: pr, R
   implicit none
contains
   subroutine helmholtz_res(n, V, T, P_SL, T_SL, rho_SL, r_SL, a_res)
       ! Argumentos
       real(pr), intent(in) :: n
       real(pr), intent(in) :: V
       real(pr), intent(in) :: T

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
       real(pr), intent(out) :: a_res
   

       ! Cuerpo de la subrutina
       T_red = T / T_SL
       rho_red = ((1 / V) / rho_SL)
       
       a_res = r_SL * n *  R * T * (&
               (- rho_red / T_red) &
               + (1.0_pr / rho_SL) - 1.0_pr &
               * log(1.0_pr - rho_SL) &
               + 1.0_pr)

   end subroutine helmholtz_res
end module a_residual