module autodiff_hyperdual_SanchezLacombe
   use yaeos__constants, only: pr, R
   use yaeos__adiff_hyperdual_ar_api
   use yaeos__models_solvers, only: volume_michelsen
   implicit none

   type, extends(ArModelAdiff) :: hd_SL
      real(pr), allocatable :: kij(:, :), lij(:, :), segment(:), vol(:), eppsilon(:)
   contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
      procedure :: volume => volume
   end type

contains

   type(hd_SL) function setup(segment, vol, eppsilon, kij, lij) result(self) ! parametros de cada compuesto, kij y lij
      !! Seup an Autodiff_Sanchez-Lacombe model
      real(pr) :: segment(:)
      real(pr) :: vol(:)
      real(pr) :: eppsilon(:)
      real(pr) :: kij(:, :)
      real(pr) :: lij(:, :)

      self%segment = segment
      self%vol = vol
      self%eppsilon = eppsilon

      self%kij = kij
      self%lij = lij
   end function

   function arfun(self, n, V, T) result(ar)
      class(hd_SL) :: self
      type(hyperdual), intent(in) :: n(:), v, t
      type(hyperdual) :: ar

      type(hyperdual) :: segment_mix
      type(hyperdual) :: vol_mix
      type(hyperdual) :: eppsilon_mix
      type(hyperdual) :: frac_segment(size(n))
      type(hyperdual) :: x(size(n))
      type(hyperdual) :: frac_segment_ij
      type(hyperdual) :: rho_SL
      type(hyperdual) :: T_SL
      type(hyperdual) :: T1, T2, T3, T4

      integer :: i, j

      associate(segment => self%segment, vol => self%vol, &
         eppsilon => self%eppsilon, &
         kij => self%kij, lij => self%lij &
         )
         ! ==================================================================
         ! Fraction segment and mixture segment
         ! ------------------------------------------------------------------
         x = n / sum(n)
         segment_mix = sum(segment * x)

         frac_segment = (x * segment) / segment_mix

         ! ==================================================================
         ! Regla de mezclado
         ! ------------------------------------------------------------------
         vol_mix = 0.0_pr
         eppsilon_mix = 0.0_pr
         do i=1,size(n)-1
            do j=i+1,size(n)
               frac_segment_ij = frac_segment(i) * frac_segment(j)
               vol_mix = vol_mix + (1.0_pr / 2.0_pr) * frac_segment_ij * (vol(i) + vol(j)) * (1 - lij(i, j))
               eppsilon_mix = eppsilon_mix + frac_segment_ij * sqrt(eppsilon(i) * eppsilon(j)) * (1 - kij(i, j))
            end do
         end do
         ! Añadido de la diagonal de la regla de mezclado
         vol_mix = vol_mix + sum(frac_segment**2 * vol)
         eppsilon_mix = eppsilon_mix + sum(frac_segment**2 * eppsilon)


         ! ==================================================================
         ! Funcion Ar
         ! ------------------------------------------------------------------
         rho_SL = segment_mix * vol_mix / v ! densidad reducida
         T_SL = R * t / eppsilon_mix ! temperatura reducida

         T1 = segment_mix * sum(n) *  R * t 
         T2 = (- rho_SL / T_SL)
         T3 = (1.0_pr / rho_SL) - 1.0_pr
         T4 = log(1.0_pr - rho_SL)

         ar = T1 * (T2 + T3 * T4 + 1.0_pr) 

      end associate
   end function

   function v0(self, n, p, t)
      !! Initialization of volume
      class(hd_SL), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: vol_mix
      real(pr) :: frac_segment_ij
      real(pr) :: x(size(n))
      real(pr) :: segment_mix
      real(pr) :: frac_segment(size(n))
      real(pr) :: v0

      integer :: i, j

      x = n / sum(n)
      segment_mix = sum(self%segment * x)

      frac_segment = (x * self%segment) / segment_mix

      ! ==================================================================
      ! Regla de mezclado
      ! ------------------------------------------------------------------
      vol_mix = 0.0_pr

      do i=1,size(n)-1
         do j=i+1,size(n)
            frac_segment_ij = frac_segment(i) * frac_segment(j)
            vol_mix = vol_mix + (1.0_pr / 2.0_pr) * frac_segment_ij * (self%vol(i) + self%vol(j)) * (1 - self%lij(i, j))
         end do
      end do
      ! Añadido de la diagonal de la regla de mezclado
      vol_mix = vol_mix + sum(frac_segment**2 * self%vol)

      !------!
      !--v0--!
      !------!
      v0 = vol_mix * segment_mix
   end function

   subroutine volume(eos, n, P, T, V, root_type)
      class(hd_SL), intent(in) :: eos
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), intent(in) :: root_type

      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]

      call volume_michelsen(eos, n, P, T, V, root_type)

   end subroutine
end module
