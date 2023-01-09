subroutine init 
    use longr
    implicit none 
    Ndim = 2
    Nobs = 2
    Nsigpt = 5
    Npar = 1
    sigma_type = "usual"
    deltat=1.D-2
    rho=2.D-1
    rho_perturb = 1.D-1
    alpha = 1.D-1
    alpha_perturb = 3.D-1
    beta = 5.D-2
    beta_perturb = 3.D-2
    gamma = 2.D-2
    gamma_perturb = 1.D-2
    kappa = 0.D0
    upsilon = 1.D-3
end subroutine init