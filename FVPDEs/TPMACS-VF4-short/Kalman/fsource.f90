module fsource
    implicit none 
    contains 

    function state_transition(x)
        use longr
        implicit none 
        real(kind=long),dimension(Ndim), intent(in) :: x 
        real(kind=long), dimension(Ndim) :: state_transition 

        state_transition(1) = x(1) + deltat*rho_perturb*x(1)
        state_transition(2) = x(2) + deltat*alpha*x(1) -deltat*beta*x(2) &
        & - deltat*gamma*x(1)*x(2)
    end function state_transition

    subroutine sigma_points(weights,points,m,cov)
        use longr
        implicit none 
        real(kind=long), intent(in) :: m(:) 
        real(kind=long), dimension(size(m),size(m)), intent(in) :: cov(:,:)
        real(kind=long), dimension(2*size(m)+1), intent(out) :: weights
        real(kind=long), dimension(2*size(m)+1,size(m)), intent(out) :: points
        real(kind=long), dimension(size(m),size(m)) :: root_cov
        integer :: q
        select case(sigma_type)
        case("usual")
            Nsigpt = 2*(size(m))+1
            Lscale = upsilon-size(m)
            root_cov = sqrt_matrix(cov)
            weights = 0.5D0/(size(m)+Lscale)
            weights(Nsigpt) = Lscale/(size(m)+Lscale)
            do q=1,size(m)
                points(q,:) = m + sqrt(upsilon)*root_cov(q,:)
                points(size(m)+q,:) = m - sqrt(upsilon)*root_cov(q,:)
            end do
            points(Nsigpt,:) = m
        end select
    end subroutine sigma_points 

    subroutine unscented_transform_obser(mY,cYY,cXY,mX,cXX)
        use longr
        implicit none 
        real(kind=long), intent(in) ::mX(:)
        real(kind=long),dimension(size(mX,1),size(mX,1)), intent(in) :: cXX
        real(kind=long),dimension(Nobs), intent(out) :: mY
        real(kind=long),dimension(size(mX,1),Nobs), intent(out) :: cXY
        real(kind=long),dimension(Nobs,Nobs), intent(out) :: cYY
        real(kind=long),dimension(2*(size(mX,1))+1) :: weights
        real(kind=long), dimension(2*(size(mX,1))+1,size(mX,1)) :: points 
        integer :: q
        real(kind=long), dimension(size(mX,1)) :: vec_x
        real(kind=long), dimension(Nobs) :: vec_y
        call sigma_points(weights,points,mX,cXX)
            !print*,'weights = ',weights
            !print*,'points = ',points
        mY = 0.D0
        cYY=0.D0
        cXY=0.D0
        ! y est l observation
        do q=1,2*(size(mX,1))+1
            mY = mY + weights(q)*phi(points(q,:))
        end do
        weights(2*(size(mX,1))+1) = weights(2*(size(mX,1))+1) + 1 - upsilon**2 + 2
        do q=1,2*(size(mX,1))+1 
            vec_y = phi(points(q,:))-mY
                !print*,'vec_y = ',vec_y
            vec_x = points(q,:)-mX
                !print*,'vec_x = ',vec_x
            cYY = cYY + weights(q)*prod_vec_to_mat(vec_y,vec_y)
            cXY = cXY + weights(q)*prod_vec_to_mat(vec_x,vec_y)
                !print*,q,weights(q)*prod_vec_to_mat(vec_x,vec_y)
        end do
                !print*,'cXY',cXY
                !print*,'cYY',cYY
                !print*,'cXY.cYY-1',matmul(cXY,transpose(cYY))
    end subroutine unscented_transform_obser 

    subroutine unscented_transform_state(mY,cYY,mX,cXX,mTheta)
        use longr
        implicit none 
        real(kind=long), intent(in) ::  mX(:)
        real(kind=long), dimension(Npar), intent(in) :: mTheta
        real(kind=long),dimension(size(mX,1),size(mX,1)), intent(in) :: cXX
        real(kind=long),dimension(size(mX,1)), intent(out) :: mY
        real(kind=long),dimension(size(mX,1),size(mX,1)), intent(out) :: cYY
        real(kind=long),dimension(2*(size(mX,1))+1) :: weights
        real(kind=long), dimension(2*(size(mX,1))+1,size(mX,1)) :: points 
        integer :: q
        real(kind=long), dimension(size(mX,1)) :: vec_y
        call sigma_points(weights,points,mX,cXX)
            !print*,'weights = ',weights
            !print*,'points = ',points
        mY = 0.D0
        cYY=0.D0
        !y est l ensemble des variables 
        do q=1,2*(size(mX,1))+1
            mY = mY + weights(q)*GOO(mTheta,points(q,:))
        end do
            !print*,'x_suiv = ',mY
        weights(2*(size(mX,1))+1) = weights(2*(size(mX,1))+1) + 1 - upsilon**2 + 2
        do q=1,2*(size(mX,1))+1
            vec_y = GOO(mTheta,points(q,:))-mY
            cYY = cYY + weights(q)*prod_vec_to_mat(vec_y,vec_y)
        end do
    end subroutine unscented_transform_state

    function phi(point)
        use longr 
        implicit none 
        real(kind=long),dimension(Ndim),intent(in) :: point
        real(kind=long), dimension(Nobs) :: phi 
        phi = point(1:Nobs)
    end function phi

    function exact_solution(t)
        use longr
        implicit none
        real(kind=long), intent(in) :: t 
        real(kind=long), dimension(Ndim) :: exact_solution 
        real(kind=long) :: F0
        integer :: I
        real(kind=long), dimension(10000) :: x,y
        F0 = F_primi(0.D0)
        x = (/ (I, I = 1, 10000) /)*t/10000
        y = (/ (g(I*t/10000)*exp(-F_primi(I*t/10000)),I = 1,10000)/)
        exact_solution(1) = x_0(1)+rho*t!x_0(1)*exp(rho*t)
        exact_solution(2) = x_0(2)*exp(F_primi(t)-F0) +exp(F_primi(t))*integrate(x,y)
    end function exact_solution

    function f(t)
        use longr
        implicit none 
        real(kind=long), intent(in) :: t 
        real(kind=long) :: f 
        f = -beta - gamma *exp(x_0(1))*exp(rho*t)
        !-beta - gamma *x_0(1)*exp(rho*t)
    end function f

    function F_primi(t)
        use longr
        implicit none 
        real(kind=long), intent(in) :: t 
        real(kind=long) :: F_primi
        F_primi = -beta*t-gamma*exp(x_0(1))/rho*exp(rho*t)+gamma*exp(x_0(1))/rho
        !-beta*t-gamma*x_0(1)/rho*exp(rho*t)+gamma*x_0(1)/rho
    end function F_primi

    function g(t) 
        use longr
        implicit none 
        real(kind=long), intent(in) :: t 
        real(kind=long) :: g
        g = alpha*exp(x_0(1))*exp(rho*t)
        !alpha*x_0(1)*exp(rho*t)
    end function g

    function integrate(x, y) result(r)
        use longr
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(kind=long),dimension(10000), intent(in)  :: x      !! Variable x
    real(kind=long),dimension(10000), intent(in)  :: y   !! Function y(x)
    real(kind=long)              :: r            !! Integral ∫y(x)·dx
    integer :: n
    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
    end function

    function inv(A) result(Ainv)
        use longr
        implicit none
        real(kind=long),intent(in) :: A(:,:)
        real(kind=long)            :: Ainv(size(A,1),size(A,2))
        real(kind=long)            :: work(size(A,1))            ! work array for LAPACK
        integer         :: n,info,ipiv(size(A,1))     ! pivot indices
    
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
        call dpotrf	('L',n,Ainv,n,info)	
        if (info.ne.0) stop 'Matrix is numerically singular!'
        call dpotri('L',n,Ainv,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function inv

    function MatVecProduct(A,x)
        use longr 
        implicit none 
        real(kind=long), intent(in) :: A(:,:) 
        real(kind=long), dimension(size(A,2)), intent(in) :: x 
        real(kind=long), dimension(size(A,1)) :: MatVecProduct 
        integer :: i
        do i=1,size(A,1)
            MatVecProduct(i) = sum(A(i,:)*x)
        end do
    end function MatVecProduct

    function sqrt_matrix(Y)
        use longr 
        implicit none 
        real(kind=long), intent(in) :: Y(:,:)
        real(kind=long) :: sqrt_matrix(size(Y,1),size(Y,2))
        integer :: info,n
        n = size(Y,1)
        sqrt_matrix = Y
        call dpotrf	('L',n,sqrt_matrix,n,info)	

    end function sqrt_matrix

    function prod_vec_to_mat(u,v)
        use longr 
        implicit none 
        real(kind=long), intent(in) :: u(:)
        real(kind=long), intent(in) :: v(:)
        real(kind=long), dimension(size(u),size(v)) :: prod_vec_to_mat
        integer :: i,j 
        do i=1,size(u)
            do j=1,size(v)
            prod_vec_to_mat(i,j) = u(i)*v(j)
            end do
        end do
    end function prod_vec_to_mat

    function normale(N) !generate N values of N(0,1)
        use longr 
        implicit none 
        integer, intent(in) :: N
        real(kind=long), dimension(N) :: normale 
        integer :: i 
        real(kind=long) :: u,v 
        real(kind=long), parameter :: pi=4.D0*DATAN(1.D0)
        do i=1,N 
            call random_number(u)
            call random_number(v)
            normale(i) = sqrt(-2*log(u))*cos(2*pi*v)
        end do
    end function normale

    function GOO(Param,State)
        use longr 
        implicit none 
        real(kind=long), dimension(Npar), intent(in) :: Param 
        real(kind=long), dimension(Ndim), intent(in) :: State 
        real(kind=long), dimension(Nobs) :: GOO
        GOO(1) = State(1) + deltat*exp(Param(1))!*State(1)
        GOO(2) = State(2) + alpha*deltat*exp(State(1))-beta*deltat*State(2) &
        &- gamma*deltat*exp(State(1))*State(2)
        !Param(2)*deltat*exp(State(1))-Param(3)*deltat*State(2) - Param(4)*deltat*exp(State(1))*State(2)
        !alpha*deltat*exp(State(1))-beta*deltat*State(2) - gamma*deltat*exp(State(1))*State(2)
        !alpha*deltat*State(1)-beta*deltat*State(2) - gamma*deltat*State(1)*State(2)
    end function GOO

    subroutine unscented_transform_param(cTT,cThetaY,cYY,Theta,cThetaTheta,Xprev)
        use longr
        implicit none 
        real(kind=long), dimension(Npar), intent(in) ::  Theta
        real(kind=long),dimension(Npar,Npar), intent(in) :: cThetaTheta
        real(kind=long), dimension(Ndim), intent(in) :: Xprev
        real(kind=long),dimension(Nobs,Nobs), intent(out) :: cYY
        real(kind=long), dimension(Npar,Nobs), intent(out) :: cThetaY
        real(kind=long), dimension(Npar,Npar), intent(out) :: cTT
        real(kind=long), dimension(Npar) ::  mTheta
        real(kind=long), dimension(Nobs) :: mY
        real(kind=long),dimension(2*Npar+1) :: weights
        real(kind=long), dimension(2*Npar+1,Npar) :: points 
        integer :: q
        real(kind=long), dimension(Npar) :: vec_theta
        real(kind=long), dimension(Nobs) :: vec_y
        call sigma_points(weights,points,Theta,cThetaTheta)
            !print*,'weights = ',weights
            !print*,'points = ',points
        mY = 0.D0!GOO(points(2*Npar+1,:),Xprev)!
        cTT = 0.D0
        cThetaY = 0.D0
        cYY = 0.D0
        mTheta = Theta
        !y est l ensemble des variables 
        do q=1,2*Npar+1
            mY = mY + weights(q)*GOO(points(q,:),Xprev)
        end do
            !print*,'x_suiv = ',mY
        weights(2*Npar+1) = weights(2*Npar+1) + 1 + 2 - upsilon/(Npar+kappa)
        do q=1,2*Npar+1
            vec_theta = points(q,:)-mTheta
                !print*,'vec_theta',vec_theta
            vec_y = GOO(points(q,:),Xprev)-mY
                !print*,'vec_y',vec_y
            cTT = cTT + weights(q)*prod_vec_to_mat(vec_theta,vec_theta)
            cThetaY = cThetaY + weights(q)*prod_vec_to_mat(vec_theta,vec_y)
            cYY = cYY + weights(q)*prod_vec_to_mat(vec_y,vec_y)
        end do
    end subroutine unscented_transform_param


    subroutine unscented_transform_param_forecast(Theta_f,Var_f,Err,Theta_p,Var_p)
        use longr
        implicit none 
        real(kind=long), dimension(Npar), intent(in) ::  Theta_p
        real(kind=long),dimension(Npar,Npar), intent(in) :: Var_p,Err
        real(kind=long), dimension(Npar,Npar), intent(out) :: Var_f
        real(kind=long), dimension(Npar), intent(out) ::  Theta_f
        real(kind=long),dimension(2*Npar+1) :: weights
        real(kind=long), dimension(2*Npar+1,Npar) :: points 
        integer :: q
        real(kind=long), dimension(Npar) :: vec_theta
        call sigma_points(weights,points,Theta_p,Var_p)
            !print*,'weights = ',weights
            !print*,'points = ',points
        Theta_f = 0.D0!Theta_p!
        Var_f = 0.D0
        !y est l ensemble des variables 
        do q=1,2*Npar+1
            Theta_f = Theta_f + weights(q)*points(q,:)
        end do
            !print*,'Theta_f = ',Theta_f
        weights(2*Npar+1) = weights(2*Npar+1) + 1 + 2 - upsilon/(Npar+kappa)
        do q=1,2*Npar+1
            vec_theta = points(q,:)-Theta_f
            Var_f = Var_f + weights(q)*prod_vec_to_mat(vec_theta,vec_theta)
        end do
        Var_f = Var_f + Err
    end subroutine unscented_transform_param_forecast

end module fsource 