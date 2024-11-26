program flash_SL
   use autodiff_hyperdual_SanchezLacombe
   use yaeos__constants, only: pr
   use yaeos__adiff_hyperdual_ar_api
   use yaeos__models_solvers, only: volume_michelsen

   use yaeos, only: EquilibriumState, flash
   implicit none

   type(hd_SL) :: model_SL
   type(hyperdual) :: Ar

   type(EquilibriumState) :: flash_result !! Result of Flash calculation

   real(pr) :: n(2), t, p, k0(2)
   real(pr) :: r_SL(2), vol_SL(2), epsilon_SL(2)
   real(pr) :: kij(2, 2), lij(2, 2)
   integer :: iter

   ! Code starts here:
   print *, "FLASH Sanchez-Lacombe test:"
   print *, "Methane/Butane mixture"

   ! Methane/ Butane mixture

   n = [0.4, 0.6]                      ! Composition
   r_SL = [4.26, 7.59]                 ! segment SL
   vol_SL = [0.00752, 0.0104]          ! volume SL
   epsilon_SL = [18.6244, 33.5073]     ! epsilon SL


   ! The mixrule
   kij = 0.0_pr
   lij = 0.0_pr

   ! Conditions
   T = 294.0_pr ! K
   P = 60.0_pr ! bar

   ! Sanchez Lacombe model with autodiff

   model_SL = setup(r_SL, vol_SL, epsilon_SL, kij, lij)

   ! K-wilson factors for initialization
   k0 = [5.1875, 3.5659E-002]
   flash_result = flash(model_SL, n, t=t, p_spec=p, k0=k0, iters=iter)
   
   ! Print results with format statements
   print "(A,5x, *(F6.4,1x))", "X:", flash_result%x
   print "(A,5x, *(F6.4,1x))", "Y:", flash_result%y
   print "(A,1x, *(E15.7))", "Vx: ", flash_result%Vx
   print "(A,1x, *(E15.7))", "Vy: ", flash_result%Vy
   print "(A,1x, F10.3)", "P: ", flash_result%p
   print "(A,1x, F10.3)", "T: ", flash_result%T

end program
