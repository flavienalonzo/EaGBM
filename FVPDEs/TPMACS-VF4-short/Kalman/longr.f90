module longr
  implicit none
  integer, parameter     :: long = 8
  character (len=6)      :: prefix
  integer                :: Ndim,Nsigpt,Npar,Nsigpt_par,Nobs
  real (kind = long)     :: Lscale, deltat,rho,alpha,beta,gamma,rho_perturb,alpha_perturb,beta_perturb,gamma_perturb
  real(kind=long)        :: kappa, upsilon
  logical, dimension(:), allocatable :: isbomba
  real(kind=long), dimension(:),allocatable :: x_0
  real(kind=long),dimension(:,:), allocatable :: covXX
  CHARACTER (len = 5)   :: sigma_type, ChoixPb

end module longr
