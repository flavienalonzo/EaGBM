program kalman
    use longr
    use fsource
    implicit none 
    real(kind=long), dimension(:,:), allocatable :: R,Q,err_pred,err_prev,cYY_pred,cXY,cYY,Var_pred,err_suiv
    real(kind=long), dimension(:), allocatable :: x_pred,y_pred,x_prev,y_obs,x_suiv, y_interm, x_interm,x_noK
    real(kind=long), dimension(:,:), allocatable :: cInterm, err_Interm,Cov_pred,cfact,V_f,V_a,V_prev,Ctt,Cty,Kt
    real(kind=long), dimension(:), allocatable :: theta_0,theta_f,theta_a,theta_prev,y_true,d_t
    real(kind=long), dimension(:,:), allocatable :: Mat1,Mat2,Q_state,cYY_state
    integer :: n
    call init
    open(unit=666,file='data_Kalman.txt')

    ChoixPb = "UD_KF"
    if (ChoixPb=="state") then
        !!! Modèle : 
        ! Y(n) = h_n(x(n)) + w(n)               !Observations
        ! X(n+1) = f_n(x(n)) + v(n)             !State
        !! Filtre de Kalman unscented (UKF) with additive noise
        allocate(x_0(Ndim),x_pred(Ndim),y_pred(Nobs),x_prev(Ndim),y_obs(Nobs),x_suiv(Ndim))
        allocate(Q(Ndim,Ndim),R(Nobs,Nobs),err_pred(Ndim,Ndim),theta_0(Npar))
        allocate(err_prev(Ndim,Ndim),cYY_pred(Nobs,Nobs),cXY(Ndim,Nobs),cYY(Nobs,Nobs),y_interm(Nobs))
        allocate(Var_pred(Ndim,Ndim),err_suiv(Ndim,Ndim),cInterm(Ndim,Nobs))
        allocate(x_interm(Ndim),err_Interm(Ndim,Ndim),x_noK(Ndim),cfact(Ndim,Ndim))

        !Initialisation 
        x_0 = (/log(1.D-3),1.D0/)
        x_noK = x_0
        theta_0 = (/rho_perturb/)
        R = 0.D0; R(1,1)=deltat**2; R(2,2)=deltat**2
        Q = 0.D0; Q(1,1) = 1.D-2; Q(2,2) = 1.D-2
        cfact = 0.D0; cfact(1,1)=1.D0; cfact(2,2)=1.D0
        write(666,*) 0.D0, x_0(1), x_0(2), exact_solution(0.D0), x_noK(1), x_noK(2)
        !prediction du premier etat
        x_pred = GOO(theta_0,x_0)
        err_pred = Q
        !call unscented_transform_state(x_suiv,Var_pred,x_pred,err_pred)
        x_prev = x_pred 
        err_prev = err_pred
        do n=1,5000
            !Prediction observation
            call unscented_transform_obser(y_pred,cYY_pred,cXY,x_prev,err_prev)
                !print*,'x_prev = ',x_prev
                !print*,'err_prev',err_prev
                !print*,'y_pred = ',y_pred
                !print*,'cYY_pred = ',cYY_pred
                !print*,'cXY = ',cXY
            cYY = cYY_pred+ R
                !print*,'cYY = ',cYY
            !Observation update
            y_obs = exact_solution(n*deltat)
                !print*,'y_obs = ',y_obs
            !Estimation state
            if (modulo(n,100)==0)  then
                cInterm = matmul(cXY,inv(cYY))
                y_interm = y_obs-y_pred
                x_interm = MatVecProduct(cInterm,y_interm)
                    !print*,'C_interm = ',cInterm(1:Ndim,1)
                    !print*,'C_interm = ',cInterm(1:Ndim,2)
                    !print*,'y_interm = ',y_interm
                    !print*,'x_interm = ',x_interm
                print*,n*deltat
            else 
                x_interm = 0.D0
                cInterm = matmul(cXY,inv(cYY))
            end if
            x_pred = x_prev + x_interm
            err_Interm = matmul(cInterm,transpose(cXY))
                print*,'err_interm = ',err_Interm
            err_pred = err_prev
            err_pred = err_pred - err_Interm
                print*,'x_pred = ',x_pred
                print*,'err_pred = ',err_pred
            !Prediction of future state
            call unscented_transform_state(x_suiv,Var_pred,x_pred,err_pred,theta_0)
            err_suiv = Var_pred + Q 
                print*,'x_suiv = ',x_suiv
                print*,'err_suiv = ',err_suiv
            !Actualisation
            x_noK = GOO(theta_0,x_noK)
            x_prev = x_suiv
            err_prev = err_suiv 
            if (AINT(n*deltat+deltat)-AINT(n*deltat)>0) then
                write(666,*) n*deltat, x_suiv(1), x_suiv(2), y_obs(1), y_obs(2), x_noK(1), x_noK(2)
            end if
        end do
        print*,'done'
        close(666)
        deallocate(x_0,x_pred,y_pred,x_prev,y_obs,x_suiv,Q,R,err_pred,err_prev,cYY_pred,cXY,cYY,y_interm)
        deallocate(Var_pred,err_suiv,cInterm,x_interm,err_Interm,x_noK,cfact)

    else if (ChoixPb=="param") then
        !!!Modèle : 
        ! Y(n) = h_n(x(n)) + w(n)               !Observations
        ! X(n+1) = f_n(x(n),theta)              !State depending on a parameter
        allocate(x_0(Ndim),theta_0(Npar),x_prev(Ndim),theta_prev(Npar),Q(Npar,Npar),V_f(Npar,Npar),x_pred(Ndim))
        allocate(theta_f(Npar),V_prev(Npar,Npar),Ctt(Npar,Npar),Cty(Npar,Nobs),Cyy(Nobs,Nobs),V_a(Npar,Npar))
        allocate(R(Nobs,Nobs),y_true(Nobs),y_pred(Nobs),d_t(Nobs),Kt(Npar,Nobs),theta_a(Npar),Mat1(Npar,Nobs),Mat2(Npar,Npar))
        x_0 = (/log(1.D-3),1.D0/)
        theta_0 = (/rho_perturb/)
        !theta_0 = (/rho_perturb,alpha_perturb,beta_perturb,gamma_perturb/)
        Q = 1.D-3
        !Q = 0.D0; Q(1,1)=1.D-3; Q(2,2)=1.D-3; Q(3,3)=1.D-3; Q(4,4)=1.D-3;
        R = 0.D0; R(1,1)=deltat**2; R(2,2)=deltat**2
        !write(666,*) 'Time , ', 'Simulated data , ',' ', 'Observations , ',' ', 'Parameters calibration','Variance'
        write(666,*) 0.D0, x_0(1), x_0(2), exact_solution(0.D0),theta_0,Q 
        x_prev = x_0
        theta_prev = theta_0
        theta_a = theta_0
        V_prev = Q
        do n=1,5000
            if (modulo(n,500)==0) then
                !((AINT(9*n*deltat)-AINT(9*n*deltat-deltat)>0)) then!.and.(modulo(int(n*deltat),5)==0)) then
            !Prediction step (forecast)
                print*,'time',n*deltat
                print*,'theta_prev',theta_prev
            call unscented_transform_param_forecast(Theta_f,V_f,Q,theta_prev,V_prev)
            !theta_f = theta_prev !+ MatVecProduct(Q,normale(Npar))
                print*,'theta_f',theta_f
            !V_f = V_prev + Q
                print*,'V_f',V_f

            !Observation
            y_true = exact_solution(n*deltat)
                print*,'y_true',y_true
            y_pred = GOO(theta_f,x_prev)
                print*,'y_pred',y_pred
            d_t = y_true - y_pred

            !Update step (analysis)
            call unscented_transform_param(Ctt,Cty,Cyy,theta_f,V_f,x_prev)
                print*,'Ctt',Ctt
                print*,'Cty',Cty
                print*,'Cyy',Cyy
            Kt = matmul(Cty,inv(Cyy+R))
                print*,'Kt',Kt
            theta_a = theta_f + MatVecProduct(Kt,d_t)
                print*,'theta_a',theta_a
            Mat1 = matmul(Kt,Ctt)
            Mat2 = matmul(Mat1,transpose(Kt))
            V_a = V_f - Mat2!Ctt + Q
                print*,'V_a',V_a
            x_pred = GOO(theta_a,x_prev)
            else
                !No measurement available
            theta_a = theta_prev
            x_pred = GOO(theta_a,x_prev)
            V_a = V_prev
            end if
            !Actualisation algo
            x_prev = x_pred
            theta_prev = theta_a 
            V_prev = V_a 

            if (AINT(100*n*deltat)-AINT(100*n*deltat-deltat)>0) then
                print*,deltat*n
                write(666,*) n*deltat, x_prev(1), x_prev(2), exact_solution(n*deltat),theta_a,V_a
            end if
        end do
        print*,'done'
        close(666)
        deallocate(x_0,theta_0,x_prev,theta_prev,Q,V_f,x_pred)
        deallocate(theta_f,V_prev,Ctt,Cty,Cyy,V_a)
        deallocate(R,y_true,y_pred,d_t,Kt,theta_a,Mat1,Mat2)

    else if (ChoixPb=="UD_KF") then
        !!!Modèle : 
        ! Y(n) = h_n(x(n)) + w(n)               !Observations
        ! X(n+1) = f_n(x(n),theta) + v(n)           !State depending on a parameter
        allocate(x_0(Ndim),theta_0(Npar),x_prev(Ndim),theta_prev(Npar),Q(Npar,Npar),V_f(Npar,Npar),x_pred(Ndim))
        allocate(theta_f(Npar),V_prev(Npar,Npar),Ctt(Npar,Npar),Cty(Npar,Nobs),Cyy(Nobs,Nobs),V_a(Npar,Npar),Q_state(Ndim,Ndim))
        allocate(R(Nobs,Nobs),y_true(Nobs),y_pred(Nobs),d_t(Nobs),Kt(Npar,Nobs),theta_a(Npar),Mat1(Npar,Nobs),Mat2(Npar,Npar))
        allocate(err_prev(Ndim,Ndim),cYY_pred(Nobs,Nobs),cXY(Ndim,Nobs),y_interm(Nobs),x_suiv(Ndim),x_interm(Ndim))
        allocate(cInterm(Ndim,Nobs),cYY_state(Nobs,Nobs),err_pred(Ndim,Ndim),err_interm(Ndim,Ndim))
        allocate(err_suiv(Ndim,Ndim),Var_pred(Ndim,Ndim))
        x_0 = (/log(1.D-3),1.D0/)
        theta_0 = (/log(rho_perturb)/)
        !theta_0 = (/rho_perturb,alpha_perturb,beta_perturb,gamma_perturb/)
        Q = 1.D-3
        !Q = 0.D0; Q(1,1)=1.D-3; Q(2,2)=1.D-3; Q(3,3)=1.D-3; Q(4,4)=1.D-3;
        Q_state = 0.D0; Q_state(1,1) = 1.D-3; Q_state(2,2) = 1.D-3;
        R = 0.D0; R(1,1)=deltat**2; R(2,2)=deltat**2
        !write(666,*) 'Time , ', 'Simulated data , ',' ', 'Observations , ',' ', 'Parameters calibration','Variance'
        write(666,*) 0.D0, x_0(1), x_0(2), exact_solution(0.D0),exp(theta_0)!,Q 
        x_prev = x_0
        theta_prev = theta_0
        theta_a = theta_0
        V_prev = Q
        err_prev = Q_state
        do n=1,5000
            if (modulo(n,250)==0) then
                !((AINT(9*n*deltat)-AINT(9*n*deltat-deltat)>0)) then!.and.(modulo(int(n*deltat),5)==0)) then
            !Prediction step (forecast)
                print*,'time',n*deltat
                print*,'theta_prev',theta_prev
            call unscented_transform_param_forecast(Theta_f,V_f,Q,theta_prev,V_prev)
            !theta_f = theta_prev !+ MatVecProduct(Q,normale(Npar))
                print*,'theta_f',theta_f
            !V_f = V_prev + Q
                print*,'V_f',V_f
            call unscented_transform_obser(y_pred,cYY_pred,cXY,x_prev,err_prev)
            cYY_state = cYY_pred + Q_state
            cInterm = matmul(cXY,inv(cYY_state))
                    
            !Observation
            y_true = exact_solution(n*deltat)
                print*,'y_true',y_true
            !y_pred = GOO(theta_f,x_prev)
                print*,'y_pred',y_pred
            d_t = y_true - y_pred


            x_interm = MatVecProduct(cInterm,d_t)
            x_pred = x_prev + x_interm
            !Update step (analysis)
            call unscented_transform_param(Ctt,Cty,Cyy,theta_f,V_f,x_pred)
                print*,'Ctt',Ctt
                print*,'Cty',Cty
                print*,'Cyy',Cyy
            Kt = matmul(Cty,inv(Cyy+R))
                print*,'Kt',Kt
            theta_a = theta_f + MatVecProduct(Kt,d_t)
                print*,'theta_a',exp(theta_a)
            Mat1 = matmul(Kt,Ctt)
            Mat2 = matmul(Mat1,transpose(Kt))
            V_a = V_f - Mat2!Ctt + Q
                print*,'V_a',V_a

            
            err_Interm = matmul(cInterm,transpose(cXY))
                !print*,'err_interm = ',err_Interm
            err_pred = err_prev
            err_pred = err_pred - err_Interm
            call unscented_transform_state(x_suiv,Var_pred,x_pred,err_pred,theta_a)
            err_suiv = Var_pred + Q_state
            !x_pred = GOO(theta_a,x_prev)
            else
                !No measurement available
            theta_a = theta_prev
            x_suiv= GOO(theta_a,x_prev)
            V_a = V_prev 
            err_suiv = err_prev
            end if
            !Actualisation algo
            x_prev = x_suiv
            theta_prev = theta_a 
            V_prev = V_a 
            err_prev = err_suiv

            if (AINT(100*n*deltat)-AINT(100*n*deltat-deltat)>0) then
                print*,deltat*n
                write(666,*) n*deltat, x_suiv, exact_solution(n*deltat),exp(theta_a) !,V_a
            end if
        end do
        print*,'done'
        close(666)
        deallocate(x_0,theta_0,x_prev,theta_prev,Q,V_f,x_pred)
        deallocate(theta_f,V_prev,Ctt,Cty,Cyy,V_a)
        deallocate(R,y_true,y_pred,d_t,Kt,theta_a,Mat1,Mat2,Q_state,err_prev,cYY_pred,cXY,y_interm,x_suiv,x_interm)
        deallocate(cInterm,cYY_state,err_pred,err_interm,err_suiv,Var_pred)
    end if

end program kalman