program flash_SL_polymer
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
   real(pr) :: mw(2), P_SL(2), T_SL(2), rho_SL(2), wi(2)
   real(pr) :: r_SL(2), vol_SL(2), epsilon_SL(2)
   real(pr) :: kij(2, 2), lij(2, 2)

   real(pr) :: phase_1_mass(2), phase_2_mass(2), phase_1_mol(2), phase_2_mol(2)
   integer :: iter

   ! Code starts here:
   print *, "FLASH Sanchez-Lacombe for a polymer:"
   print *, "pentane / polyethylene(16400) mixture"

   ! carbon dioxide / PMMA mixture

   wi = [0.9, 0.1]                     ! mass fraction
   mw = [72.15, 16400.0]               ! molecular weight
   P_SL = [3101.0, 3590.0]
   T_SL = [441.0, 650.0]
   rho_SL = [755.0, 895.0] / mw        ! units: mol/L


   ! variables i must to use
   n = (wi / mw) / sum(wi / mw)        ! Composition
   epsilon_SL = T_SL * R               ! epsilon SL
   vol_SL = epsilon_SL / P_SL          ! volume SL
   r_SL = 1.0_pr / (rho_SL * vol_SL)   ! segment SL


   ! The mixrule
   kij = 0.0_pr
   lij = 0.0_pr

   ! Conditions
   T = 460.0_pr ! K
   P = 300.0_pr ! bar

   ! Sanchez Lacombe model with autodiff

   model_SL = setup(r_SL, vol_SL, epsilon_SL, kij, lij)

   ! initialization
   phase_1_mass = [1.0 - 7.44e-03, 7.443-03]      ! mass fraction
   phase_2_mass = [1.0 - 1.977e-01, 1.977e-01]    ! mass fraction

   phase_1_mol = (phase_1_mass * mw) / sum(phase_1_mass * mw)
   phase_2_mol = (phase_2_mass * mw) / sum(phase_2_mass * mw)


   k0 = phase_2_mol / phase_1_mol
   print *, 'k0 = ', k0

   flash_result = flash(model_SL, n, t=t, p_spec=p, k0=k0, iters=iter)
   
   ! Print results with format statements
   print "(A,5x, *(F6.4,1x))", "X:", flash_result%x
   print "(A,5x, *(F6.4,1x))", "Y:", flash_result%y
   print "(A,1x, *(E15.7))", "Vx: ", flash_result%Vx
   print "(A,1x, *(E15.7))", "Vy: ", flash_result%Vy
   print "(A,1x, F10.3)", "P: ", flash_result%p
   print "(A,1x, F10.3)", "T: ", flash_result%T

end program