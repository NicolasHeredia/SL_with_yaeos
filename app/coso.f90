program main
  use autodiff_hyperdual_SanchezLacombe
  use yaeos__constants, only: pr
  use yaeos__adiff_hyperdual_ar_api
  use yaeos__models_solvers, only: volume_michelsen

  implicit none

  real(pr) :: n(1), segments(1), volumes(1), epsilons(1)
  real(pr) :: kij(1, 1), lij(1, 1)

  real(pr) :: V, T, P
  real(pr) :: phi_v(1)
  type(hd_SL) :: model
  type(hyperdual) :: Ar

  integer :: i

  ! Abrir archivo directamente
  open(unit=10, file="output.txt", status="replace")

  ! Parametros del metano en las unidades correspondientes
  n = [1.0_pr]
  segments = [4.26_pr]
  volumes = [0.00752_pr] ! L * bar
  epsilons = [18.6244_pr] ! L*bar / mol

  ! The mixrule
  kij = 0.0_pr
  lij = 0.0_pr

  ! Conditions
  T = 175.0_pr ! K
  ! P = 60.0_pr ! bar

  ! Crear el modelo y calcular Ar
  model = setup(segments, volumes, epsilons, kij, lij)

  call model%lnphi_pt(n, P, T, lnPhi=phi_v, root_type="stable")
  print *, phi_v

  do i=1,1000
    V = real(i, pr)/1000
    call model%pressure(n, V, T, P=P)
    write(10, '(E20.12, 2X, E20.12)') V, P
  end do

  ! Cerrar el archivo
  close(10)

  P = 60.0_pr ! bar

  call model%volume(n, P=P, T=T, V=V, root_type="stable")
  print *, V, P, 1.0_pr / V

  call model%pressure(n, V=V, T=T, P=P)
  print *, P, T

end program