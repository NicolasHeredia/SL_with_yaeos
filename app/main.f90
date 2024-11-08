program main
  use autodiff_hyperdual_SanchezLacombe
  use yaeos__constants, only: pr
  use yaeos__adiff_hyperdual_ar_api
  use yaeos__models_solvers, only: volume_michelsen

  implicit none

  integer :: n_components
  real(pr) :: n(2), segments(2), volumes(2), epsilons(2)
  real(pr) :: kij(2, 2), lij(2, 2)
  real(pr) :: dPdV, dPdT, dPdn(2)
  real(pr) :: m(2), phi_v(2)

  real(pr) :: V, T, P
  type(hd_SL) :: model
  type(hyperdual) :: Ar

  integer :: i

  ! Parametros
  n = [0.4_pr, 0.6_pr]
  segments = [5.30_pr, 11.49_pr]
  volumes = [4.41_pr, 7.69_pr]
  epsilons = [2709.0_pr, 3850.1_pr]

  ! The mixrule
  kij = 0.0_pr
  lij = 0.0_pr

  ! Conditions
  T = 300.0_pr
  P = 10.0_pr

  ! Crear el modelo y calcular Ar
  model = setup(segments, volumes, epsilons, kij, lij)

  call model%lnphi_pt(n, P, T, lnPhi=phi_v, root_type="stable")
  print *, phi_v

  V = 5
  ! do i=1,5000
  !  V = real(i, pr)/1000
  !  call model%pressure(n, V, T, P=P)
  !  write(1, *) V, P
  !end do

  call model%volume(n, P=P, T=T, V=V, root_type="vapor")
  print *, V, P


end program