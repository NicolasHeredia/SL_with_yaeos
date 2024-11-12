program main
  use autodiff_hyperdual_SanchezLacombe
  use yaeos__constants, only: pr
  use yaeos__adiff_hyperdual_ar_api
  use yaeos__models_solvers, only: volume_michelsen

  implicit none

  integer :: n_components
  real(pr) :: n(1), segments(1), volumes(1), epsilons(1)
  real(pr) :: kij(1, 1), lij(1, 1)
  real(pr) :: dPdV, dPdT, dPdn(1)


  real(pr) :: V, T, P
  type(hd_SL) :: model
  type(hyperdual) :: Ar

  integer :: i

  ! Parametros
  n = [1.0_pr]
  segments = [4.26_pr]
  volumes = [7.52_pr]
  epsilons = [1864.96_pr]

  ! The mixrule
  kij = 0.0_pr
  lij = 0.0_pr

  ! Conditions
  T = 110.5_pr
  P = 46.04_pr

  ! Crear el modelo y calcular Ar
  model = setup(segments, volumes, epsilons, kij, lij)

  !call model%lnphi_pt(n, P, T, lnPhi=phi_v, root_type="stable")
  !print *, phi_v

  V = 1.0
  ! do i=1,5000
  !  V = real(i, pr)/1000
  !  call model%pressure(n, V, T, P=P)
  !  write(1, *) V, P
  !end do

  call model%volume(n, P=P, T=T, V=V, root_type="stable")
  print *, V, P

  call model%pressure(n, V=V, T=T, P=P)
  print *, P, T

end program