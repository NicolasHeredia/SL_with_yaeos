program flasher
   !! Example program to make flash calculations with a PengRobinson76 
   !! EoS model.

   ! Import the relevant 
   use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel, k_wilson
   implicit none

   ! Variables definition:
   class(ArModel), allocatable :: model !! Model to use
   type(EquilibriumState) :: flash_result !! Result of Flash calculation

   real(pr) :: tc(2), pc(2), w(2)
   real(pr) :: n(2), t, p, k0(2)
   integer :: iter

   ! Code starts here:
   print *, "FLASH EXAMPLE:"
   print *, "Methane/Butane mixture"
   
   ! Methane/ Butane mixture
   n = [0.4, 0.6]                      ! Composition
   tc = [190.564, 425.12]              ! Critical temperatures
   pc = [45.99, 37.96]                 ! Critical pressures
   w = [0.0115478, 0.200164]           ! Acentric factors

   ! Use the PengRobinson76 model
   model = PengRobinson76(tc, pc, w)

   ! Set pressure and temperatures
   P = 60
   T = 294

   ! Calculate flashes
   !flash_result = flash(model, n, t=t, p_spec=p, iters=iter)
   !write(*, *) flash_result
   
   ! K-wilson factors for initialization
   k0 = k_wilson(model, T=T, P=P)
   !flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iter)
   ! Print results with format statements

   
   ! Specify volume
   flash_result = flash(model, n, t=t, V_spec=1.0_pr, iters=iter)
   !write(*, *) flash_result

   print "(A,5x, *(F6.4,1x))", "X:", flash_result%x
   print "(A,5x, *(F6.4,1x))", "Y:", flash_result%y
   print "(A,1x, *(E15.7))", "Vx: ", flash_result%Vx
   print "(A,1x, *(E15.7))", "Vy: ", flash_result%Vy
   print "(A,1x, F10.3)", "P: ", flash_result%p
   print "(A,1x, F10.3)", "T: ", flash_result%T

   print *, 'k-wilson:', k0
end program