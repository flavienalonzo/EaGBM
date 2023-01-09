SUBROUTINE NewtonAngioP1sommets(Atum,Anut,Aend,Aveg,TumOld,NutrimOld,EndothOld,VasegfOld,&
   &Tum,Nut,Endo,Vegf,dim,choixf,temps,dX)
    !--------
    ! Modules
    !--------
    USE longr
    USE imprime
    USE parmmage
    USE intmatvec
    USE intbigradc
    USE fsourcemod
    use fsourcebreast
  
    IMPLICIT NONE
  
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux)                                   :: Atum,Anut,Aend,Aveg
    TYPE(MatCreux)                                   :: Atumnut,Atumend,Atumveg
    TYPE(MatCreux)                                   :: Anuttum,Anutend,Anutveg
    TYPE(MatCreux)                                   :: Aendtum,Aendnut,Aendveg
    TYPE(MatCreux)                                   :: Avegtum,Avegnut,Avegend
    type(MatCreux),dimension(4,4)                    :: Tab_A
    Integer, intent(in)                              :: dim, choixf
    REAL(kind=long), intent(in)                      :: temps
    REAL(kind = long), dimension(dim), intent(inout) :: Tum,Nut,Endo,Vegf 
    REAL(kind=long), dimension (dim), intent(in)    :: TumOld,NutrimOld,EndothOld,VasegfOld
    REAL(kind=long), dimension (4*dim), intent(out)  :: dX
  
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: jt, ik, j, k, jL, iseg, ii, kiter, kitermax, is , js ,i
    REAL(kind=long)       :: CoefAjout, x1, y1
    REAL(kind = long)     :: Adegenbord
    REAL(kind = long)     :: dVKLmoins, dVKLplus
    REAL(kind = long)     :: dVKL, dVLK, Ubord, Fbord, temp
    REAL(kind = long), dimension(4*dim)            :: FF
    REAL(kind = long), dimension(:), ALLOCATABLE     :: T
    REAL(kind = long), dimension(3,3)                :: Aloc
    REAL(kind = long), dimension(2)                  :: DerAKL,DerEtaKL
    !--------------------------------------------------------
    !--------------------------------------------------------
    INTERFACE
       FUNCTION amatloc(jt)
         USE longr
         IMPLICIT NONE
         INTEGER, INTENT(in)               :: jt       ! numero du triangle
         REAL(kind=long), DIMENSION(3,3)   :: amatloc
       END FUNCTION amatloc
    END INTERFACE
  
    !-------------------
    ! Debut du programme
    !-------------------
    kitermax = 10000
    Tum=TumOld; Nut=NutrimOld; Endo=EndothOld; Vegf=VasegfOld;
    call matrixinitP1sommets(Atumnut)
    call matrixinitP1sommets(Atumend)
    call matrixinitP1sommets(Atumveg)
    call matrixinitP1sommets(Anuttum)
    call matrixinitP1sommets(Anutend)
    call matrixinitP1sommets(Anutveg)
    call matrixinitP1sommets(Aendtum)
    call matrixinitP1sommets(Aendnut)
    call matrixinitP1sommets(Aendveg)
    call matrixinitP1sommets(Avegtum)
    call matrixinitP1sommets(Avegnut)
    call matrixinitP1sommets(Avegend)

    do kiter = 1, kitermax

        !Initialisation des matrices
        Atum%TMAT=0.D0;Anut%TMAT=0.D0;Aend%TMAT=0.D0;Aveg%TMAT=0.D0
        Atumnut%TMAT=0.D0;Atumend%TMAT=0.D0;Atumveg%TMAT=0.D0
        Anuttum%TMAT=0.D0;Anutend%TMAT=0.D0;Anutveg%TMAT=0.D0
        Aendtum%TMAT=0.D0;Aendnut%TMAT=0.D0;Aendveg%TMAT=0.D0
        Avegtum%TMAT=0.D0;Avegnut%TMAT=0.D0;Avegend%TMAT=0.D0
        FF = 0.D0
        !Termes sans derivees spatiales:

        DO i = 1, dim
            !Equation en u:
            FF(i) = FF(i) + AireDsommet(i)*( Tum(i)-TumOld(i) )/dt -AireDsommet(i)*(- apop*Tum(i)&
            & -chemotherapy(temps)*Tum(i))- AireDsommet(i)*rate*Tum(i)*(1-Tum(i)-Endo(i))

            call ajout(i,i, AireDSommet(i)/dt + AireDSommet(i)*apop+chemotherapy(temps)&
            &- AireDsommet(i)*rate*(1-2*Tum(i)-Endo(i)),Atum)
            call ajout(i,i,AireDsommet(i)*rate*Tum(i),Atumend)

            if (isbomba(i).eqv..true.) then 
                FF(i) = FF(i) + AireDsommet(i)*radiotherapy(temps)*Tum(i)
                call ajout(i,i,AireDsommet(i)*radiotherapy(temps),Atum)
            end if          

            !Equation en c:
            FF(dim+i) = FF(dim+i) + AireDsommet(i)* (Nut(i)-NutrimOld(i))/dt &
            & - AireDsommet(i)*(Coef_prod*Endo(i) - nut_degra*Nut(i)-Coef_cons*Tum(i)*Nut(i) )

            call ajout(i,i,AireDsommet(i)/dt + AireDsommet(i)*nut_degra+AireDsommet(i)*Coef_cons*Tum(i),Anut)
            call ajout(i,i,-AireDsommet(i)*Coef_prod,Anutend)
            call ajout(i,i,AireDsommet(i)*Coef_cons*Nut(i),Anuttum)

            !Equation en ue:
            FF(2*dim+i) = FF(2*dim+i) + AireDsommet(i)*( Endo(i)-EndothOld(i) )/dt &
            & -AireDsommet(i)*( rate_endo*Endo(i)*(1-Endo(i)-Tum(i)) - degr_endo*Endo(i) )

            call ajout(i,i,AireDsommet(i)/dt -AireDsommet(i)*rate_endo*(1-Tum(i)-2*Endo(i)) &
            & +AireDsommet(i)*degr_endo,Aend)
            call ajout(i,i,AireDsommet(i)*rate_endo*Endo(i),Aendtum)

            !Equation en V:
            FF(3*dim+i) = FF(3*dim+i) + AireDsommet(i)* ( Vegf(i)-VasegfOld(i) )/dt - AireDsommet(i)*&
            & (- VEGF_degr*Vegf(i)-VEGF_cons*Endo(i)*Vegf(i) )
            call ajout(i,i, AireDsommet(i)/dt+ AireDsommet(i)*VEGF_degr+ AireDsommet(i)*VEGF_cons*Endo(i),Aveg)
            call ajout(i,i, AireDsommet(i)*VEGF_cons*Vegf(i),Avegend)
            if (Nut(i)<=seuil_hypo.and.Nut(i)>seuil_necro) then
                FF(3*dim+i) = FF(3*dim+i) - AireDsommet(i)*VEGF_prod*Tum(i)
                call ajout(i,i,- AireDsommet(i)*VEGF_prod,Avegtum)
            end if
        end do

        !Termes de diffusion:
        !equation en u:
        CoefTranspChi = chi_u; CoefDiffuAdeg = Diff_u; 
        call tenseur(choixanisu,temps)
        DO jt =1,Nbt
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
    
            Aloc = - amatloc(jt) 
    
            ! contribution dans la ligne i
            !-----------------------------
           
            IF ( Ntyps(i) /= 1) THEN
               IF (Ntyps(j) /= 1 ) THEN 
                  temp = AKL(Aloc(1,2), Tum(i), Tum(j))
                  FF(i) = FF(i) + Aloc(1,2)*temp*(p(Tum(i)) - p(Tum(j)))
                  DerAKL = DerivAKL(Aloc(1,2), Tum(i), Tum(j))
                  !! derivee de FF(i) par rapport a Tum(i)
                  CALL Ajout (i, i, Aloc(1,2)*( DerAKL(1)*(p(Tum(i)) - p(Tum(j))) + temp*Derivp(Tum(i)) ), Atum )
                  !! derivee de FF(i) par rapport a Tum(j)
                  CALL Ajout (i, j, Aloc(1,2)*( DerAKL(2)*(p(Tum(i)) - p(Tum(j))) - temp*Derivp(Tum(j)) ), Atum )
               ENDIF  ! fin (i,j) 
               !
               IF (Ntyps(k) /= 1 ) THEN
                  temp = AKL(Aloc(1,3), Tum(i), Tum(k))
                  FF(i) = FF(i) + Aloc(1,3)*temp*(p(Tum(i)) - p(Tum(k)))
                  DerAKL = DerivAKL(Aloc(1,3), Tum(i), Tum(k))
                  !! derivee de FF(i) par rapport a Tum(i)
                  CALL Ajout (i, i, Aloc(1,3)*( DerAKL(1)*(p(Tum(i)) - p(Tum(k))) + temp*Derivp(Tum(i)) ), Atum )
                  !! derivee de FF(i) par rapport a Tum(k)
                  CALL Ajout (i, k, Aloc(1,3)*( DerAKL(2)*(p(Tum(i)) - p(Tum(k))) - temp*Derivp(Tum(k)) ), Atum )
               ENDIF ! fin (i,k) interieur
            END IF  ! fin i 
            !print*,'max FF apres i',maxval(FF)
            ! contribution dans la ligne j 
            !-----------------------------
            
            IF ( Ntyps(j) /= 1) THEN
               IF (Ntyps(i) /= 1 ) THEN
                  temp = AKL(Aloc(2,1), Tum(j), Tum(i)) 
                  FF(j) = FF(j) + Aloc(2,1)*temp*(p(Tum(j)) - p(Tum(i)))
                  DerAKL = DerivAKL(Aloc(2,1), Tum(j), Tum(i))
                  !! derivee de FF(j) par rapport a Tum(j)
                  CALL Ajout (j, j, Aloc(2,1)*( DerAKL(1)*(p(Tum(j)) - p(Tum(i))) + temp*Derivp(Tum(j)) ), Atum )
                  !! derivee de FF(j) par rapport a Tum(i)
                  CALL Ajout (j, i, Aloc(2,1)*( DerAKL(2)*(p(Tum(j)) - p(Tum(i))) - temp*Derivp(Tum(i)) ), Atum )
               ENDIF  ! fin (j,i) 
               !
               IF (Ntyps(k) /= 1 ) THEN
                  temp = AKL(Aloc(2,3), Tum(j), Tum(k))
                  FF(j) = FF(j) + Aloc(2,3)*temp*(p(Tum(j)) - p(Tum(k)))
                  DerAKL = DerivAKL(Aloc(2,3), Tum(j), Tum(k))
                  !! derivee de FF(j) par rapport a Tum(j)
                  CALL Ajout (j, j, Aloc(2,3)*( DerAKL(1)*(p(Tum(j)) - p(Tum(k))) + temp*Derivp(Tum(j)) ), Atum )
                  !! derivee de FF(j) par rapport a Tum(k)
                  CALL Ajout (j, k, Aloc(2,3)*( DerAKL(2)*(p(Tum(j)) - p(Tum(k))) - temp*Derivp(Tum(k)) ), Atum )
               ENDIF ! fin (j,k) 
            END IF  ! fin j
            !print*,'max FF apres j',maxval(FF)
            ! contribution dans la ligne k
            !-----------------------------
            
            IF ( Ntyps(k) /= 1) THEN
               IF (Ntyps(i) /= 1 ) THEN
                  temp = AKL(Aloc(3,1), Tum(k), Tum(i))
                  FF(k) = FF(k) + Aloc(3,1)*temp*(p(Tum(k)) - p(Tum(i)))
                  DerAKL = DerivAKL(Aloc(3,1), Tum(k), Tum(i))
                  !! derivee de FF(k) par rapport a Tum(k)
                  CALL Ajout (k, k, Aloc(3,1)*( DerAKL(1)*(p(Tum(k)) - p(Tum(i))) + temp*Derivp(Tum(k)) ), Atum )
                  !! derivee de FF(k) par rapport a Tum(i)
                  CALL Ajout (k, i, Aloc(3,1)*( DerAKL(2)*(p(Tum(k)) - p(Tum(i))) - temp*Derivp(Tum(i)) ), Atum )
               ENDIF  ! fin (k,i) 
               !
               IF (Ntyps(j) /= 1 ) THEN
                  temp = AKL(Aloc(3,2), Tum(k), Tum(j))
                  FF(k) = FF(k) + Aloc(3,2)*temp*(p(Tum(k)) - p(Tum(j)))
                  DerAKL = DerivAKL(Aloc(3,2), Tum(k), Tum(j))
                  !! derivee de FF(k) par rapport a Tum(k)
                  CALL Ajout (k, k, Aloc(3,2)*( DerAKL(1)*(p(Tum(k)) - p(Tum(j))) + temp*Derivp(Tum(k)) ), Atum )
                  !! derivee de FF(k) par rapport a Tum(j)
                  CALL Ajout (k, j, Aloc(3,2)*( DerAKL(2)*(p(Tum(k)) - p(Tum(j))) - temp*Derivp(Tum(j)) ), Atum )
               ENDIF ! fin (k,j) 
            END IF  !! fin k 
         END DO
         deallocate(Sxxk,SxyK,SyyK)

         !equation en c:
         CoefDiffV = Coef_diffusion;
        
         call tenseur(choixanis,temps)

         DO jt =1,Nbt
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
    
            Aloc = - Amatloc(jt)*CoefDiffV !! We multiply by CoefDiffV
            
            !-----------------------------
            ! contribution dans la ligne i
            !-----------------------------
            
            IF ( NtypS(i) /= 1) THEN
               IF (NtypS(j) /= 1 ) THEN 
                  temp = EtaKL(Aloc(1,2), Nut(i), Nut(j))
                  FF(dim+i) = FF(dim+i) + Aloc(1,2)*temp*(p(Nut(i)) - p(Nut(j)))
                  DerEtaKL = DerivEtaKL(Aloc(1,2), Nut(i), Nut(j))
                  !! derivee de FF(i) par rapport a Nut(i)
                  CALL Ajout (i, i, Aloc(1,2)*(DerEtaKL(1)*(p(Nut(i))-p(Nut(j)))+temp*Derivp(Nut(i))),Anut)
                  !! derivee de FF(i) par rapport a Nut(j)
                  CALL Ajout (i, j, Aloc(1,2)*(DerEtaKL(2)*(p(Nut(i))-p(Nut(j)))-temp*Derivp(Nut(j))),Anut)
               ENDIF  ! fin (i,j) 
               !
               IF (NtypS(k) /= 1 ) THEN
                  temp = EtaKL(Aloc(1,3), Nut(i), Nut(k))
                  FF(dim+i) = FF(dim+i) + Aloc(1,3)*temp*(p(Nut(i))-p(Nut(k)))
                  DerEtaKL = DerivEtaKL(Aloc(1,3), Nut(i), Nut(k))
                  !! derivee de FF(i) par rapport a Nut(i)
                  CALL Ajout (i, i, Aloc(1,3)*(DerEtaKL(1)*(p(Nut(i))-p(Nut(k)))+temp*Derivp(Nut(i))),Anut)
                  !! derivee de FF(i) par rapport a Nut(k)
                  CALL Ajout (i, k, Aloc(1,3)*(DerEtaKL(2)*(p(Nut(i))-p(Nut(k)))-temp*Derivp(Nut(k))),Anut)
               ENDIF ! fin (i,k) interieur
            END IF  ! fin i 
            
            ! contribution dans la ligne j 
            !-----------------------------
            
            IF ( NtypS(j) /= 1) THEN
               IF (NtypS(i) /= 1 ) THEN
                  temp = EtaKL(Aloc(2,1), Nut(j), Nut(i)) 
                  FF(dim+j) = FF(dim+j) + Aloc(2,1)*temp*(p(Nut(j)) - p(Nut(i)))
                  DerEtaKL = DerivEtaKL(Aloc(2,1), Nut(j), Nut(i))
                  !! derivee de FF(j) par rapport a Nut(j)
                  CALL Ajout (j, j, Aloc(2,1)*( DerEtaKL(1)*(p(Nut(j)) - p(Nut(i))) + temp*Derivp(Nut(j)) ), Anut )
                  !! derivee de FF(j) par rapport a Nut(i)
                  CALL Ajout (j, i, Aloc(2,1)*( DerEtaKL(2)*(p(Nut(j)) - p(Nut(i))) - temp*Derivp(Nut(i)) ), Anut )
               ENDIF  ! fin (j,i) 
               !
               IF (NtypS(k) /= 1 ) THEN
                  temp = EtaKL(Aloc(2,3), Nut(j), Nut(k))
                  FF(dim+j) = FF(dim+j) + Aloc(2,3)*temp*(p(Nut(j)) - p(Nut(k)))
                  DerEtaKL = DerivEtaKL(Aloc(2,3), Nut(j), Nut(k))
                  !! derivee de FF(j) par rapport a Nut(j)
                  CALL Ajout (j, j, Aloc(2,3)*( DerEtaKL(1)*(p(Nut(j)) - p(Nut(k))) + temp*Derivp(Nut(j)) ), Anut )
                  !! derivee de FF(j) par rapport a Nut(k)
                  CALL Ajout (j, k, Aloc(2,3)*( DerEtaKL(2)*(p(Nut(j)) - p(Nut(k))) - temp*Derivp(Nut(k)) ), Anut )
               ENDIF ! fin (j,k) 
            END IF  ! fin j
            
            ! contribution dans la ligne k
            !-----------------------------
            
            IF ( NtypS(k) /= 1) THEN
               IF (NtypS(i) /= 1 ) THEN
                  temp = EtaKL(Aloc(3,1), Nut(k), Nut(i))
                  FF(dim+k) = FF(dim+k) + Aloc(3,1)*temp*(p(Nut(k)) - p(Nut(i)))
                  DerEtaKL = DerivEtaKL(Aloc(3,1), Nut(k), Nut(i))
                  !! derivee de FF(k) par rapport a Nut(k)
                  CALL Ajout (k, k, Aloc(3,1)*( DerEtaKL(1)*(p(Nut(k)) - p(Nut(i))) + temp*Derivp(Nut(k)) ), Anut )
                  !! derivee de FF(k) par rapport a Nut(i)
                  CALL Ajout (k, i, Aloc(3,1)*( DerEtaKL(2)*(p(Nut(k)) - p(Nut(i))) - temp*Derivp(Nut(i)) ), Anut )
               ENDIF  ! fin (k,i) 
               !
               IF (NtypS(j) /= 1 ) THEN
                  temp = EtaKL(Aloc(3,2), Nut(k), Nut(j))
                  FF(dim+k) = FF(dim+k) + Aloc(3,2)*temp*(p(Nut(k)) - p(Nut(j)))
                  DerEtaKL = DerivEtaKL(Aloc(3,2), Nut(k), Nut(j))
                  !! derivee de FF(k) par rapport a Nut(k)
                  CALL Ajout (k, k, Aloc(3,2)*( DerEtaKL(1)*(p(Nut(k)) - p(Nut(j))) + temp*Derivp(Nut(k)) ), Anut )
                  !! derivee de FF(k) par rapport a Nut(j)
                  CALL Ajout (k, j, Aloc(3,2)*( DerEtaKL(2)*(p(Nut(k)) - p(Nut(j))) - temp*Derivp(Nut(j)) ), Anut )
               ENDIF ! fin (k,j) 
            END IF  !! fin k 
         END DO
         deallocate(Sxxk,SxyK,SyyK)

         !equation en ue:
         CoefTranspChi = chemo_endo; CoefDiffuAdeg = Diff_endo; 
        call tenseur(choixanisu,temps)

         DO jt =1,Nbt
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
    
            Aloc = - amatloc(jt) 
    
            ! contribution dans la ligne i
            !-----------------------------
           
            IF ( Ntyps(i) /= 1) THEN
               IF (Ntyps(j) /= 1 ) THEN 
                  temp = AKL(Aloc(1,2), Endo(i), Endo(j))
                  FF(2*dim+i) = FF(2*dim+i) + Aloc(1,2)*temp*(p(Endo(i)) - p(Endo(j)))
                  DerAKL = DerivAKL(Aloc(1,2), Endo(i), Endo(j))
                  !! derivee de FF(i) par rapport a Endo(i)
                  CALL Ajout (i, i, Aloc(1,2)*( DerAKL(1)*(p(Endo(i)) - p(Endo(j))) + temp*Derivp(Endo(i)) ), Aend )
                  !! derivee de FF(i) par rapport a Endo(j)
                  CALL Ajout (i, j, Aloc(1,2)*( DerAKL(2)*(p(Endo(i)) - p(Endo(j))) - temp*Derivp(Endo(j)) ), Aend )
               ENDIF  ! fin (i,j) 
               !
               IF (Ntyps(k) /= 1 ) THEN
                  temp = AKL(Aloc(1,3), Endo(i), Endo(k))
                  FF(2*dim+i) = FF(2*dim+i) + Aloc(1,3)*temp*(p(Endo(i)) - p(Endo(k)))
                  DerAKL = DerivAKL(Aloc(1,3), Endo(i), Endo(k))
                  !! derivee de FF(i) par rapport a Endo(i)
                  CALL Ajout (i, i, Aloc(1,3)*( DerAKL(1)*(p(Endo(i)) - p(Endo(k))) + temp*Derivp(Endo(i)) ), Aend )
                  !! derivee de FF(i) par rapport a Endo(k)
                  CALL Ajout (i, k, Aloc(1,3)*( DerAKL(2)*(p(Endo(i)) - p(Endo(k))) - temp*Derivp(Endo(k)) ), Aend )
               ENDIF ! fin (i,k) interieur
            END IF  ! fin i 
            !print*,'max FF apres i',maxval(FF)
            ! contribution dans la ligne j 
            !-----------------------------
            
            IF ( Ntyps(j) /= 1) THEN
               IF (Ntyps(i) /= 1 ) THEN
                  temp = AKL(Aloc(2,1), Endo(j), Endo(i)) 
                  FF(2*dim+j) = FF(2*dim+j) + Aloc(2,1)*temp*(p(Endo(j)) - p(Endo(i)))
                  DerAKL = DerivAKL(Aloc(2,1), Endo(j), Endo(i))
                  !! derivee de FF(j) par rapport a Endo(j)
                  CALL Ajout (j, j, Aloc(2,1)*( DerAKL(1)*(p(Endo(j)) - p(Endo(i))) + temp*Derivp(Endo(j)) ), Aend )
                  !! derivee de FF(j) par rapport a Endo(i)
                  CALL Ajout (j, i, Aloc(2,1)*( DerAKL(2)*(p(Endo(j)) - p(Endo(i))) - temp*Derivp(Endo(i)) ), Aend )
               ENDIF  ! fin (j,i) 
               !
               IF (Ntyps(k) /= 1 ) THEN
                  temp = AKL(Aloc(2,3), Endo(j), Endo(k))
                  FF(2*dim+j) = FF(2*dim+j) + Aloc(2,3)*temp*(p(Endo(j)) - p(Endo(k)))
                  DerAKL = DerivAKL(Aloc(2,3), Endo(j), Endo(k))
                  !! derivee de FF(j) par rapport a Endo(j)
                  CALL Ajout (j, j, Aloc(2,3)*( DerAKL(1)*(p(Endo(j)) - p(Endo(k))) + temp*Derivp(Endo(j)) ), Aend )
                  !! derivee de FF(j) par rapport a Endo(k)
                  CALL Ajout (j, k, Aloc(2,3)*( DerAKL(2)*(p(Endo(j)) - p(Endo(k))) - temp*Derivp(Endo(k)) ), Aend )
               ENDIF ! fin (j,k) 
            END IF  ! fin j
            !print*,'max FF apres j',maxval(FF)
            ! contribution dans la ligne k
            !-----------------------------
            
            IF ( Ntyps(k) /= 1) THEN
               IF (Ntyps(i) /= 1 ) THEN
                  temp = AKL(Aloc(3,1), Endo(k), Endo(i))
                  FF(2*dim+k) = FF(2*dim+k) + Aloc(3,1)*temp*(p(Endo(k)) - p(Endo(i)))
                  DerAKL = DerivAKL(Aloc(3,1), Endo(k), Endo(i))
                  !! derivee de FF(k) par rapport a Endo(k)
                  CALL Ajout (k, k, Aloc(3,1)*( DerAKL(1)*(p(Endo(k)) - p(Endo(i))) + temp*Derivp(Endo(k)) ), Aend )
                  !! derivee de FF(k) par rapport a Endo(i)
                  CALL Ajout (k, i, Aloc(3,1)*( DerAKL(2)*(p(Endo(k)) - p(Endo(i))) - temp*Derivp(Endo(i)) ), Aend )
               ENDIF  ! fin (k,i) 
               !
               IF (Ntyps(j) /= 1 ) THEN
                  temp = AKL(Aloc(3,2), Endo(k), Endo(j))
                  FF(2*dim+k) = FF(2*dim+k) + Aloc(3,2)*temp*(p(Endo(k)) - p(Endo(j)))
                  DerAKL = DerivAKL(Aloc(3,2), Endo(k), Endo(j))
                  !! derivee de FF(k) par rapport a Endo(k)
                  CALL Ajout (k, k, Aloc(3,2)*( DerAKL(1)*(p(Endo(k)) - p(Endo(j))) + temp*Derivp(Endo(k)) ), Aend )
                  !! derivee de FF(k) par rapport a Endo(j)
                  CALL Ajout (k, j, Aloc(3,2)*( DerAKL(2)*(p(Endo(k)) - p(Endo(j))) - temp*Derivp(Endo(j)) ), Aend )
               ENDIF ! fin (k,j) 
            END IF  !! fin k 
         END DO
         deallocate(Sxxk,SxyK,SyyK)

         !equation en V:
         CoefDiffV = VEGF_dif;
        call tenseur(choixanis,temps)

        DO jt =1,Nbt
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
    
            Aloc = - Amatloc(jt)*CoefDiffV !! We multiply by CoefDiffV
            
            !-----------------------------
            ! contribution dans la ligne i
            !-----------------------------
            
            IF ( NtypS(i) /= 1) THEN
               IF (NtypS(j) /= 1 ) THEN 
                  temp = EtaKL(Aloc(1,2), Vegf(i), Vegf(j))
                  FF(3*dim+i) = FF(3*dim+i) + Aloc(1,2)*temp*(p(Vegf(i)) - p(Vegf(j)))
                  DerEtaKL = DerivEtaKL(Aloc(1,2), Vegf(i), Vegf(j))
                  !! derivee de FF(i) par rapport a Vegf(i)
                  CALL Ajout (i, i, Aloc(1,2)*( DerEtaKL(1)*(p(Vegf(i)) - p(Vegf(j))) + temp*Derivp(Vegf(i)) ), Aveg )
                  !! derivee de FF(i) par rapport a Vegf(j)
                  CALL Ajout (i, j, Aloc(1,2)*( DerEtaKL(2)*(p(Vegf(i)) - p(Vegf(j))) - temp*Derivp(Vegf(j)) ), Aveg )
               ENDIF  ! fin (i,j) 
               !
               IF (NtypS(k) /= 1 ) THEN
                  temp = EtaKL(Aloc(1,3), Vegf(i), Vegf(k))
                  FF(3*dim+i) = FF(3*dim+i) + Aloc(1,3)*temp*(p(Vegf(i)) - p(Vegf(k)))
                  DerEtaKL = DerivEtaKL(Aloc(1,3), Vegf(i), Vegf(k))
                  !! derivee de FF(i) par rapport a Vegf(i)
                  CALL Ajout (i, i, Aloc(1,3)*( DerEtaKL(1)*(p(Vegf(i)) - p(Vegf(k))) + temp*Derivp(Vegf(i)) ), Aveg )
                  !! derivee de FF(i) par rapport a Vegf(k)
                  CALL Ajout (i, k, Aloc(1,3)*( DerEtaKL(2)*(p(Vegf(i)) - p(Vegf(k))) - temp*Derivp(Vegf(k)) ), Aveg )
               ENDIF ! fin (i,k) interieur
            END IF  ! fin i 
            
            ! contribution dans la ligne j 
            !-----------------------------
            
            IF ( NtypS(j) /= 1) THEN
               IF (NtypS(i) /= 1 ) THEN
                  temp = EtaKL(Aloc(2,1), Vegf(j), Vegf(i)) 
                  FF(3*dim+j) = FF(3*dim+j) + Aloc(2,1)*temp*(p(Vegf(j)) - p(Vegf(i)))
                  DerEtaKL = DerivEtaKL(Aloc(2,1), Vegf(j), Vegf(i))
                  !! derivee de FF(j) par rapport a Vegf(j)
                  CALL Ajout (j, j, Aloc(2,1)*( DerEtaKL(1)*(p(Vegf(j)) - p(Vegf(i))) + temp*Derivp(Vegf(j)) ), Aveg )
                  !! derivee de FF(j) par rapport a Vegf(i)
                  CALL Ajout (j, i, Aloc(2,1)*( DerEtaKL(2)*(p(Vegf(j)) - p(Vegf(i))) - temp*Derivp(Vegf(i)) ), Aveg )
               ENDIF  ! fin (j,i) 
               !
               IF (NtypS(k) /= 1 ) THEN
                  temp = EtaKL(Aloc(2,3), Vegf(j), Vegf(k))
                  FF(3*dim+j) = FF(3*dim+j) + Aloc(2,3)*temp*(p(Vegf(j)) - p(Vegf(k)))
                  DerEtaKL = DerivEtaKL(Aloc(2,3), Vegf(j), Vegf(k))
                  !! derivee de FF(j) par rapport a Vegf(j)
                  CALL Ajout (j, j, Aloc(2,3)*( DerEtaKL(1)*(p(Vegf(j)) - p(Vegf(k))) + temp*Derivp(Vegf(j)) ), Aveg )
                  !! derivee de FF(j) par rapport a Vegf(k)
                  CALL Ajout (j, k, Aloc(2,3)*( DerEtaKL(2)*(p(Vegf(j)) - p(Vegf(k))) - temp*Derivp(Vegf(k)) ), Aveg )
               ENDIF ! fin (j,k) 
            END IF  ! fin j
            
            ! contribution dans la ligne k
            !-----------------------------
            
            IF ( NtypS(k) /= 1) THEN
               IF (NtypS(i) /= 1 ) THEN
                  temp = EtaKL(Aloc(3,1), Vegf(k), Vegf(i))
                  FF(3*dim+k) = FF(3*dim+k) + Aloc(3,1)*temp*(p(Vegf(k)) - p(Vegf(i)))
                  DerEtaKL = DerivEtaKL(Aloc(3,1), Vegf(k), Vegf(i))
                  !! derivee de FF(k) par rapport a Vegf(k)
                  CALL Ajout (k, k, Aloc(3,1)*( DerEtaKL(1)*(p(Vegf(k)) - p(Vegf(i))) + temp*Derivp(Vegf(k)) ), Aveg )
                  !! derivee de FF(k) par rapport a Vegf(i)
                  CALL Ajout (k, i, Aloc(3,1)*( DerEtaKL(2)*(p(Vegf(k)) - p(Vegf(i))) - temp*Derivp(Vegf(i)) ), Aveg )
               ENDIF  ! fin (k,i) 
               !
               IF (NtypS(j) /= 1 ) THEN
                  temp = EtaKL(Aloc(3,2), Vegf(k), Vegf(j))
                  FF(3*dim+k) = FF(3*dim+k) + Aloc(3,2)*temp*(p(Vegf(k)) - p(Vegf(j)))
                  DerEtaKL = DerivEtaKL(Aloc(3,2), Vegf(k), Vegf(j))
                  !! derivee de FF(k) par rapport a Vegf(k)
                  CALL Ajout (k, k, Aloc(3,2)*( DerEtaKL(1)*(p(Vegf(k)) - p(Vegf(j))) + temp*Derivp(Vegf(k)) ), Aveg )
                  !! derivee de FF(k) par rapport a Vegf(j)
                  CALL Ajout (k, j, Aloc(3,2)*( DerEtaKL(2)*(p(Vegf(k)) - p(Vegf(j))) - temp*Derivp(Vegf(j)) ), Aveg )
               ENDIF ! fin (k,j) 
            END IF  !! fin k 
         END DO
         deallocate(Sxxk,SxyK,SyyK)

         !Termes de convection:
         !equation en u:
         CoefTranspChi = chi_u; CoefDiffuAdeg = Diff_u; 
        call tenseur(choixanisu,temps)

         DO jt = 1,Nbt
        
            Aloc =  - amatloc(jt)  
            
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
                    
            IF ( Ntyps(i) /= 1) THEN       
               IF (Ntyps(j) /= 1 ) THEN 
                  dVKL = Aloc(1,2)* (Nut(j)-Nut(i))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0) 
                 
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(1,2), Tum(i), Tum(j))
                    DerAKL = DerivAKL(Aloc(1,2), Tum(i), Tum(j))
                  end select
  
                  FF(i) = FF(i) + temp*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i))) )
                  !! derivee de FF(i) par rapport a Tum(i)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i))) )+ &
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(i)) + dVKLmoins*DerivMuKLDecroit(Tum(i)))
                  call ajout(i,i, CoefAjout, Atum)
                  !! derivee de FF(i) par rapport a Tum(j)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(j)) + dVKLmoins*DerivMuKLCroit(Tum(j)) )
                  call ajout(i,j, CoefAjout, Atum)
                  !! derivee de FF(i) par rapport a Nut(i) 
                  CoefAjout = -Aloc(1,2)*temp*( dVKLplus/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(i)) + & 
                  & MuKLDecroit(Tum(j)))+dVKLmoins/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i))))
                  call ajout(i,i,CoefAjout,Atumnut)
                  !! derivee de FF(i) par rapport a Nut(j) 
                  CoefAjout = Aloc(1,2)*temp*( dVKLplus/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(i)) + & 
                  & MuKLDecroit(Tum(j)))+dVKLmoins/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i))))
                  call ajout(i,j,CoefAjout,Atumnut)
               END IF
               !!
               IF (Ntyps(k) /= 1 ) THEN 
                  dVKL = Aloc(1,3)* (Nut(k)-Nut(i))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                  !
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(1,3), Tum(i), Tum(k))
                    DerAKL = DerivAKL(Aloc(1,3), Tum(i), Tum(k))
                  end select
    
                  FF(i) = FF(i) + temp*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i))) )
                  !! derivee de FF(i) par rapport a Tum(i)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(i)) + dVKLmoins*DerivMuKLDecroit(Tum(i)))
                  call ajout(i,i, CoefAjout, Atum)
                  !! derivee de FF(i) par rapport a Tum(k)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k))+MuKLDecroit(Tum(i))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(k)) + dVKLmoins*DerivMuKLCroit(Tum(k)) )
                  call ajout(i,k, CoefAjout, Atum)
                  !! derivee de FF(i) par rapport a Nut(i)
                  CoefAjout = -Aloc(1,3)*temp*( dVKLplus/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(i)) +&
                  & MuKLDecroit(Tum(k)) )+dVKLmoins/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i))))
                  call ajout(i,i,CoefAjout,Atumnut)
                  !! derivee de FF(i) par rapport a Nut(k)
                  CoefAjout = Aloc(1,3)*temp*( dVKLplus/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(i)) +&
                  & MuKLDecroit(Tum(k)) )+dVKLmoins/(abs(dVKL)+epsilon)*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i))))
                  call ajout(i,k,CoefAjout,Atumnut)
               END IF
            END IF
            ! fin i 
            !print*,'max FF apres i',maxval(FF)
            !------------------------------
            !  Contribution dans la ligne j
            !------------------------------
    
            IF ( Ntyps(j) /= 1) THEN 
               IF (Ntyps(i) /= 1 ) THEN 
                  dVKL = Aloc(2,1)* (Nut(i)-Nut(j))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                  
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(2,1), Tum(j), Tum(i))
                    DerAKL = DerivAKL(Aloc(2,1), Tum(j), Tum(i))
                  end select
    
                  FF(j) = FF(j) + temp*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j))) )
                  !! derivee de FF(j) par rapport a Tum(j)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j))) ) + &
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(j)) + dVKLmoins*DerivMuKLDecroit(Tum(j)))
                  call ajout(j,j, CoefAjout, Atum)
                  !! derivee de FF(j) par rapport a Tum(i)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(i)) + dVKLmoins*DerivMuKLCroit(Tum(i)) )
                  call ajout(j,i, CoefAjout, Atum)
                  !! derivee de FF(j) par rapport a Tum(j) 
                  CoefAjout = -Aloc(2,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(j)) &
                  &+ MuKLDecroit(Tum(i)) )+dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j))) ) 
                  call ajout(j,j,CoefAjout,Atumnut)
                  !! derivee de FF(j) par rapport a Tum(i) 
                  CoefAjout = Aloc(2,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(j)) &
                  &+ MuKLDecroit(Tum(i)) )+dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(j))) ) 
                  call ajout(j,i,CoefAjout,Atumnut)
               END IF
               IF (Ntyps(k) /= 1 ) THEN 
                  dVKL = Aloc(2,3)* (Nut(k)-Nut(j)) 
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                 
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(2,3), Tum(j), Tum(k))
                    DerAKL = DerivAKL(Aloc(2,3), Tum(j), Tum(k))
                  end select
    
                  FF(j) = FF(j) + temp*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j))) )
                  !! derivee de FF(j) par rapport a Tum(j)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k))+MuKLDecroit(Tum(j))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(j)) + dVKLmoins*DerivMuKLDecroit(Tum(j)))
                  call ajout(j,j, CoefAjout, Atum)
                  !! derivee de FF(j) par rapport a Tum(k)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(k)) + dVKLmoins*DerivMuKLCroit(Tum(k)) )
                  call ajout(j,k, CoefAjout, Atum)
                  !! derivee de FF(j) par rapport a Nut(j) 
                  CoefAjout = -Aloc(2,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(j)) &
                  &+ MuKLDecroit(Tum(k)) ) +dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j))) )
                  call ajout(j,j,CoefAjout,Atumnut)
                  !! derivee de FF(j) par rapport a Nut(k) 
                  CoefAjout = Aloc(2,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(j)) &
                  &+ MuKLDecroit(Tum(k)) ) +dVKLmoins*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j))) )
                  call ajout(j,k,CoefAjout,Atumnut)
               END IF
            END IF  ! fin j
           ! print*,'max FF apres j',maxval(FF)
    
            !-----------------------------
            ! Contribution dans la ligne k
            !-----------------------------
    
            IF ( Ntyps(k) /= 1) THEN 
               IF (Ntyps(i) /= 1 ) THEN 
                  dVKL = Aloc(3,1)* (Nut(i)-Nut(k)) 
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(3,1), Tum(k), Tum(i))
                    DerAKL = DerivAKL(Aloc(3,1), Tum(k), Tum(i))
                  end select
                  FF(k) = FF(k) + temp*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k))) )
                  !! derivee de FF(k) par rapport a Tum(k)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(k)) + dVKLmoins*DerivMuKLDecroit(Tum(k)))
                  call ajout(k,k, CoefAjout, Atum)
                  !! derivee de FF(k) par rapport a Tum(i)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(i)) + dVKLmoins*DerivMuKLCroit(Tum(i)) )
                  call ajout(k,i, CoefAjout, Atum)
                  !! derivee de FF(k) par rapport a Tum(k) 
                  CoefAjout = -Aloc(3,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(k)) &
                  &+ MuKLDecroit(Tum(i)) ) +dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k))) )
                  call ajout(k,k,CoefAjout,Atumnut)
                  !! derivee de FF(k) par rapport a Tum(i) 
                  CoefAjout = Aloc(3,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(k)) &
                  &+ MuKLDecroit(Tum(i)) ) +dVKLmoins*( MuKLCroit(Tum(i)) + MuKLDecroit(Tum(k))) )
                  call ajout(k,i,CoefAjout,Atumnut)
               END IF
               !print*,'max FF apres ki',maxval(FF)
               !
               IF (Ntyps(j) /= 1 ) THEN 
                  dVKL = Aloc(3,2)* (Nut(j)-Nut(k))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(3,2), Tum(k), Tum(j))
                    DerAKL = DerivAKL(Aloc(3,2), Tum(k), Tum(j))
                  end select
                  FF(k) = FF(k) + temp*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k))) )
                  !! derivee de FF(k) par rapport a Tum(k)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Tum(k)) + dVKLmoins*DerivMuKLDecroit(Tum(k)))
                  call ajout(k,k, CoefAjout, Atum)
                  !! derivee de FF(k) par rapport a Tum(j)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Tum(k)) + MuKLDecroit(Tum(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Tum(j))+MuKLDecroit(Tum(k))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Tum(j)) + dVKLmoins*DerivMuKLCroit(Tum(j)) )
                  call ajout(k,j, CoefAjout, Atum)
                  !! derivee de FF(k) par rapport a Nut(k)
                  CoefAjout = -Aloc(3,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(k)) &
                  &+ MuKLDecroit(Tum(j)) ) +dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k))) )
                  call ajout(k,k,CoefAjout,Atumnut)
                  !! derivee de FF(k) par rapport a Nut(j)
                  CoefAjout = Aloc(3,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Tum(k)) &
                  &+ MuKLDecroit(Tum(j)) ) +dVKLmoins*( MuKLCroit(Tum(j)) + MuKLDecroit(Tum(k))) )
                  call ajout(k,j,CoefAjout,Atumnut)
               END IF
            END IF  ! fin k
         END DO   !! fin boucle sur jt

         deallocate(Sxxk,SxyK,SyyK)

         !equation en ue:
         CoefTranspChi = chemo_endo; CoefDiffuAdeg = Diff_endo; 
         call tenseur(choixanisu,temps)

         DO jt = 1,Nbt
        
            Aloc =  - amatloc(jt)  
            
            i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
                    
            IF ( Ntyps(i) /= 1) THEN       
               IF (Ntyps(j) /= 1 ) THEN 
                  dVKL = Aloc(1,2)* (Vegf(j)-Vegf(i))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0) 
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(1,2), Endo(i), Endo(j))
                    DerAKL = DerivAKL(Aloc(1,2), Endo(i), Endo(j))
                  end select
  
                  FF(2*dim+i) = FF(2*dim+i) + temp*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i))) )
                  !! derivee de FF(i) par rapport a Endo(i)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i))) )+ &
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(i)) + dVKLmoins*DerivMuKLDecroit(Endo(i)))
                  call ajout(i,i, CoefAjout, Aend)
                  !! derivee de FF(i) par rapport a Endo(j)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(j)) + dVKLmoins*DerivMuKLCroit(Endo(j)) )
                  call ajout(i,j, CoefAjout, Aend)
                  !! derivee de FF(i) par rapport a Vegf(i) 
                  CoefAjout = -Aloc(1,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(i)) + &
                  & MuKLDecroit(Endo(j)) ) +dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i))) )
                  call ajout(i,i,CoefAjout,Aendveg)
                  !! derivee de FF(i) par rapport a Vegf(j)
                  CoefAjout = Aloc(1,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(i)) + &
                  & MuKLDecroit(Endo(j)) ) +dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i))) )
                  call ajout(i,j,CoefAjout,Aendveg)
               END IF
               !!
               IF (Ntyps(k) /= 1 ) THEN 
                  dVKL = Aloc(1,3)* (Vegf(k)-Vegf(i))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                  !
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(1,3), Endo(i), Endo(k))
                    DerAKL = DerivAKL(Aloc(1,3), Endo(i), Endo(k))
                  end select
    
                  FF(2*dim+i) = FF(2*dim+i) + temp*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i))) )
                  !! derivee de FF(i) par rapport a Endo(i)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(i)) + dVKLmoins*DerivMuKLDecroit(Endo(i)))
                  call ajout(i,i, CoefAjout, Aend)
                  !! derivee de FF(i) par rapport a Endo(k)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k))+MuKLDecroit(Endo(i))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(k)) + dVKLmoins*DerivMuKLCroit(Endo(k)) )
                  call ajout(i,k, CoefAjout, Aend)
                  !! derivee de FF(i) par rapport a Vegf(i)
                  CoefAjout = -Aloc(1,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(i)) +&
                  & MuKLDecroit(Endo(k)) ) +dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i))) )
                  call ajout(i,i,CoefAjout,Aendveg)
                  !! derivee de FF(i) par rapport a Vegf(k)
                  CoefAjout = Aloc(1,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(i)) +&
                  & MuKLDecroit(Endo(k)) ) +dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i))) )
                  call ajout(i,k,CoefAjout,Aendveg)
               END IF
            END IF
            ! fin i 
            !print*,'max FF apres i',maxval(FF)
            !------------------------------
            !  Contribution dans la ligne j
            !------------------------------
    
            IF ( Ntyps(j) /= 1) THEN 
               IF (Ntyps(i) /= 1 ) THEN 
                  dVKL = Aloc(2,1)* (Vegf(i)-Vegf(j))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                  
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(2,1), Endo(j), Endo(i))
                    DerAKL = DerivAKL(Aloc(2,1), Endo(j), Endo(i))
                  end select 
  
                  FF(2*dim+j) = FF(2*dim+j) + temp*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j))) )
                  !! derivee de FF(j) par rapport a Endo(j)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j))) ) + &
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(j)) + dVKLmoins*DerivMuKLDecroit(Endo(j)))
                  call ajout(j,j, CoefAjout, Aend)
                  !! derivee de FF(j) par rapport a Endo(i)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(i)) + dVKLmoins*DerivMuKLCroit(Endo(i)) )
                  call ajout(j,i, CoefAjout, Aend)
                  !! derivee de FF(j) par rapport a Vegf(j) 
                  CoefAjout = -Aloc(2,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(j)) +&
                  & MuKLDecroit(Endo(i)) ) +dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j))) )
                  call ajout(j,j,CoefAjout,Aendveg)
                  !! derivee de FF(j) par rapport a Vegf(i) 
                  CoefAjout = Aloc(2,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(j)) +&
                  & MuKLDecroit(Endo(i)) ) +dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(j))) )
                  call ajout(j,i,CoefAjout,Aendveg)
               END IF
               IF (Ntyps(k) /= 1 ) THEN 
                  dVKL = Aloc(2,3)* (Vegf(k)-Vegf(j)) 
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(2,3), Endo(j), Endo(k))
                    DerAKL = DerivAKL(Aloc(2,3), Endo(j), Endo(k))
                  end select
    
                  FF(2*dim+j) = FF(2*dim+j) + temp*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j))) )
                  !! derivee de FF(j) par rapport a Endo(j)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k))+MuKLDecroit(Endo(j))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(j)) + dVKLmoins*DerivMuKLDecroit(Endo(j)))
                  call ajout(j,j, CoefAjout, Aend)
                  !! derivee de FF(j) par rapport a Endo(k)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(k)) + dVKLmoins*DerivMuKLCroit(Endo(k)) )
                  call ajout(j,k, CoefAjout, Aend)
                  !! derivee de FF(j) par rapport a Vegf(j)
                  CoefAjout = -Aloc(2,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(j)) + &
                  &MuKLDecroit(Endo(k)) ) +dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j))) )
                  call ajout(j,j,CoefAjout,Aendveg)
                  !! derivee de FF(j) par rapport a Vegf(k)
                  CoefAjout = Aloc(2,3)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(j)) + &
                  &MuKLDecroit(Endo(k)) ) +dVKLmoins*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j))) )
                  call ajout(j,k,CoefAjout,Aendveg)
               END IF
            END IF  ! fin j
           ! print*,'max FF apres j',maxval(FF)
    
            !-----------------------------
            ! Contribution dans la ligne k
            !-----------------------------
    
            IF ( Ntyps(k) /= 1) THEN 
               IF (Ntyps(i) /= 1 ) THEN 
                  dVKL = Aloc(3,1)* (Vegf(i)-Vegf(k)) 
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(3,1), Endo(k), Endo(i))
                    DerAKL = DerivAKL(Aloc(3,1), Endo(k), Endo(i))
                  end select
    
                  FF(2*dim+k) = FF(2*dim+k) + temp*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k))) )
                  !! derivee de FF(k) par rapport a Endo(k)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(k)) + dVKLmoins*DerivMuKLDecroit(Endo(k)))
                  call ajout(k,k, CoefAjout, Aend)
                  !! derivee de FF(k) par rapport a Endo(i)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(i)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(i)) + dVKLmoins*DerivMuKLCroit(Endo(i)) )
                  call ajout(k,i, CoefAjout, Aend)
                  !! derivee de FF(k) par rapport a Vegf(k)
                  CoefAjout = -Aloc(3,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(k)) +&
                  & MuKLDecroit(Endo(i)) ) +dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k))) )
                  call ajout(k,k,CoefAjout,Aendveg)
                  !! derivee de FF(k) par rapport a Vegf(i)
                  CoefAjout = Aloc(3,1)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(k)) +&
                  & MuKLDecroit(Endo(i)) ) +dVKLmoins*( MuKLCroit(Endo(i)) + MuKLDecroit(Endo(k))) )
                  call ajout(k,i,CoefAjout,Aendveg)
               END IF
               !print*,'max FF apres ki',maxval(FF)
               !
               IF (Ntyps(j) /= 1 ) THEN 
                  dVKL = Aloc(3,2)* (Vegf(j)-Vegf(k))
                  dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
    
                  SELECT CASE(choixdegenere)
                  case(0)
                    temp = CoefDiffuAdeg 
                    DerAKL = 0.D0
                  case(1)
                    temp = AKL(Aloc(3,2), Endo(k), Endo(j))
                    DerAKL = DerivAKL(Aloc(3,2), Endo(k), Endo(j))
                  end select
    
                  FF(2*dim+k) = FF(2*dim+k) + temp*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k))) )
                  !! derivee de FF(k) par rapport a Endo(k)
                  CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k))) )+&
                       & temp*(dVKLplus*DerivMuKLCroit(Endo(k)) + dVKLmoins*DerivMuKLDecroit(Endo(k)))
                  call ajout(k,k, CoefAjout, Aend)
                  !! derivee de FF(k) par rapport a Endo(j)
                  CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Endo(k)) + MuKLDecroit(Endo(j)) ) + &
                       & dVKLmoins*( MuKLCroit(Endo(j))+MuKLDecroit(Endo(k))) )+&
                       & temp*( dVKLplus*DerivMuKLDecroit(Endo(j)) + dVKLmoins*DerivMuKLCroit(Endo(j)) )
                  call ajout(k,j, CoefAjout, Aend)
                  !! derivee de FF(k) par rapport a Vegf(k)
                  CoefAjout = -Aloc(3,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(k)) +&
                  & MuKLDecroit(Endo(j)) ) +dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k))) )
                  call ajout(k,k,CoefAjout,Aendveg)
                  !! derivee de FF(k) par rapport a Vegf(k)
                  CoefAjout = Aloc(3,2)/(abs(dVKL)+epsilon)*temp*( dVKLplus*( MuKLCroit(Endo(k)) +&
                  & MuKLDecroit(Endo(j)) ) +dVKLmoins*( MuKLCroit(Endo(j)) + MuKLDecroit(Endo(k))) )
                  call ajout(k,j,CoefAjout,Aendveg)
               END IF
              ! print*,'max FF apres kj',maxval(FF)
               !
            END IF  ! fin k
           ! print*,'max FF apres k',maxval(FF)
         END DO   

         deallocate(Sxxk,SxyK,SyyK)

         Tab_A(1,1)=Atum;Tab_A(1,2)=Atumnut;Tab_A(1,3)=Atumend;Tab_A(1,4)=Atumveg;
         Tab_A(2,1)=Anuttum;Tab_A(2,2)=Anut;Tab_A(2,3)=Anutend;Tab_A(2,4)=Anutveg;
         Tab_A(3,1)=Aendtum;Tab_A(3,2)=Aendnut;Tab_A(3,3)=Aend;Tab_A(3,4)=Aendveg;
         Tab_A(4,1)=Avegtum;Tab_A(4,2)=Avegnut;Tab_A(4,3)=Avegend;Tab_A(4,4)=Aveg;

         dX = BiCGSTAB(Tab_A,-FF,4,dim,TolerenceGradient)

         If (sqrt(dot_product(dX,dX)) <TolerenceNewton) exit
         print*,'kiter = ',kiter,'     tol = ',sqrt(dot_product(dX,dX))
         print*,MAXVAL(dX(1:dim))
         print*,MAXVAL(dX(dim+1:2*dim))
         print*,MAXVAL(dX(2*dim+1:3*dim))
         print*,MAXVAL(dX(3*dim+1:4*dim))

         Tum = Tum + dX(1:dim)
         Nut = Nut + dX(dim+1:2*dim)
         Endo = Endo + dX(2*dim+1:3*dim)
         Vegf = Vegf + dX(3*dim+1:4*dim)

    end do

    Tum = Tum + dX(1:dim)
    Nut = Nut + dX(dim+1:2*dim)
    Endo = Endo + dX(2*dim+1:3*dim)
    Vegf = Vegf + dX(3*dim+1:4*dim)

    end SUBROUTINE NewtonAngioP1sommets