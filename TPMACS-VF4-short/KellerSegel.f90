program KellerSegel

  USE longr
  USE parmmage
  USE imprime
  USE intmatvec
  USE algebrelineaire
  USE intbigradc  
  use plotvtkmod
  use plotvtkmodscheme
  Use fsourcemod
  use fsourcebreast
  use parameters
  IMPLICIT NONE
  !==========================
  ! Declaration des tableaux
  !==========================
  TYPE(MatCreux)       :: A,N,Mat_E,A_vegf, Big_A, Atum, Anut, Aend, Aveg
  type(Matcreux), dimension(:,:), allocatable :: Tab_A
  !==================================
  ! Declaration des variables locales
  !==================================
  INTEGER                             :: jt, i, j, is, jv,iseg,  kiter ,k,h,kplot, ls, ii, js, kt
  REAL(kind=long), DIMENSION(:), ALLOCATABLE  :: U,U0,dU, U_p , Uexacte, Nutriment, N0 , Endothelial, Cells_density, Vasegf&
  &, Nutact, Tum, NutrimOld, TumOld, EndothOld, VasegfOld, TotCells, Matieres, Areas, dX, Tum_mid,Nut_mid,Endo_mid,Vegf_mid,&
  & TotCells_mid
   real (kind=long), dimension(:,:), allocatable :: Tab_U, Tab_U_p
   character(len=6) , dimension(:), allocatable :: Tab_entity,Tab_equa
   integer, dimension(:), allocatable :: Tab_chemo
  REAL(kind = long)               :: tol, seuil, approx, resi, temps, resiold, maxresi, resinorm,&
  & resi1, resi2, resi3, resi4, resinorm1, resinorm2, resinorm3, resinorm4
  character(:), ALLOCATABLE :: str, strPb
  logical :: usemethodmul, dolikebreast, calc_implicit

  !===================
  ! Debut du programme
  !===================
  prefix = 'LAPLAC  '
  !
  CALL init                    !* Initialisation
  print*,  'init ok '

  print*, 'Keller Segel'

  ! Lecture du maillage � partir d'un fichier : MAILLAGEGEOx, x=1...6

  CALL readmesh
  !CALL readmatlab
   !call segments
  print*,  'readmesh ok '
  call NodeConnectivity
  print*, 'NodeConnectivity ok', NmaxT

  
  usemethodmul=.false.
   dolikebreast = .true.
   calc_implicit=.true.

   if (usemethodmul.eqv..true.) then
      ALLOCATE(Nutriment(Nbt), N0(Nbt),Nutact(1:Nbt))
      allocate(Endothelial(Nbt),Vasegf(Nbt))
      allocate(Cells_density(Nbt))
      ALLOCATE(VitesseSeg(1:2,1:Nseg))
      allocate(Tab_A(n_enty,n_enty))

      do k=1,n_enty
         do h=1,n_enty
            call matrixinitVF4 ( Tab_A(k,h) ) 
         end do
      end do
      call bigmatrix(Big_A,Tab_A(1,1))

      do i=1,Nbt 
         if (ntypt(i)==100) then
            Endothelial(i) = 8.D-1
         else 
            Endothelial(i) = 1.D-1
         end if
         !Endothelial(i) = fsource( coordK(1,i), coordK(2,i), choixpb ) 
      end do

         ! resolution du syst�me lin�aire  par la mathode du gradient conjugu�
      tol = 5.D-12
      N0(:)=1.D0

      allocate(U(1:Nbt),U0(1:Nbt),dU(1:Nbt),U_p(1:Nbt),Tum(1:Nbt))
      U0 = 0.D0 
      do i=1,Nbt
         if (ntypt(i)==100) then
            U0(i) = 0.D0
         else if (ntypt(i)==200) then
            U0(i) = 0.D0
         else if (ntypt(i)==300) then
            U0(i) = 1.D-1
         else if (ntypt(i)==400) then
            U0(i) = 4.D-1
         else if (ntypt(i)==500) then
            U0(i) = 2.5D-1
         end if
         !U0(i) = cond_ini(CoordK(1,i),CoordK(2,i))
         !Nutriment(i) = 0.5*fsource(CoordK(1,i),CoordK(2,i), ChoixPb)
      end do 
      U = U0

      call matrixinitVF4( N )
      call assembleNutri(N,Endothelial,U0)
      Nutriment = bigradient(N, N%Bg,N0,tol)
      call vide( N )
      call assembleNutri(N,U0,Endothelial)
      Vasegf = 0.D0!bigradient(N, N%Bg, N0,tol)

      allocate(Tab_U(n_enty,Nbt),Tab_U_p(n_enty,Nbt),Tab_entity(n_enty),Tab_equa(n_enty),Tab_chemo(n_enty))
      Tab_U(1,:) = U0;              Tab_U_p(1,:) = U0;            Tab_entity(1)='u_norm'; Tab_equa(1)='instat';
      Tab_chemo(1) = index_nut;
      Tab_U(2,:) = Nutriment;       Tab_U_p(2,:) = Nutriment;     Tab_entity(2)='Nutrim'; Tab_equa(2)='instat';
      Tab_chemo(2) = 0;
      Tab_U(3,:) = Endothelial;     Tab_U_p(3,:) = Endothelial;   Tab_entity(3)='Endoth'; Tab_equa(3)='instat';
      Tab_chemo(3) = index_vegf;
      Tab_U(4,:) = Vasegf;          Tab_U_p(4,:) = Vasegf;        Tab_entity(4)='VasEGF'; Tab_equa(4)='instat';
      Tab_chemo(4) = 0;

      Cells_density = Tab_U(index_norm,:) + Tab_U(index_endo,:) 

      strPb = num2string(WhichPb)

      call plot_vtk(U0,'TC_Pb'//strPb//'_0','Tumor_cells')
      call plot_vtk(Nutriment,'N_Pb'//strPb//'_0','Nutrients')
      call plot_vtk(Cells_density,'CD_Pb'//strPb//'_0','Cells_density')
      call plot_vtk(Tab_U(index_vegf,:),'VEGF_Pb'//strPb//'_0','VEGF')
      call plot_vtk(Endothelial,'E_Pb'//strPb//'_0','Endothelial_cells')
      U = U0
   
   
      print*, 'matrix init ok' 

      call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,Big_A,Tab_equa,Tab_chemo,'instat',tol,0.D0)
      Tab_U_p=Tab_U
      print*,'methnewton ok'
      i=1
      kplot = 1
      do while(i*dt<=Tf)
         call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,Big_A,Tab_equa,Tab_chemo,'instat',tol,i*delta)
         Cells_density = Tab_U(index_norm,:) + Tab_U(index_endo,:) 
         Tab_U_p = Tab_U
         i=i+1
         if (modulo(i,10)==1) then
            print*,i*delta
         end if
         if (modulo(i,100)==1) then 
            str = num2string(kplot)
            call plot_vtk(Tab_U(index_nut,:),'N_Pb'//strPb//'_'//str,'Nutrients')
            call plot_vtk(Tab_U(index_norm,:),'TC_Pb'//strPb//'_'//str,'Tumor_cells')
            call plot_vtk(Cells_density,'CD_Pb'//strPb//'_'//str,'Cells_density')
            call plot_vtk(Tab_U(index_vegf,:),'VEGF_Pb'//strPb//'_'//str,'VEGF')
            call plot_vtk(Tab_U(index_endo,:),'E_Pb'//strPb//'_'//str,'Endothelial_cells')
            print*,i*delta
            kplot = kplot + 1
         end if
      end do

   else if (dolikebreast.eqv..false.) then
      i=1
      kplot = 1
      do while(i*delta<=Tf) 
         !Actualisation des nutriments
         call assemInstaNutri(Nutact,Endothelial,U,Nutriment)
         !Actualisation des cellules tumorales
         call assemInstaTumor(Tum,Endothelial,U,Nutriment) 
         Nutriment = Nutact
         U = Tum
   
         !Les plots
         !if (modulo(i,10)==0) then
         !print*,i*delta
         !end if
         if (modulo(i,5000)==0) then 
            str = num2string(kplot)
            call plot_vtk(Nutriment,'N_Pb'//strPb//'_'//str,'Nutrients')
            call plot_vtk(U,'TC_Pb'//strPb//'_'//str,'Tumor_cells')
            !call plot_vtk(Cells_density,'CD_Pb'//strPb//'_'//str,'Cells_density')
            print*,i*delta
            kplot = kplot + 1
         end if
         i=i+1
      end do 
   else if (dolikebreast.eqv..true.) then
      call changevertices 

      strPb = num2string(WhichPb)
      NbInc = Nbs!NsInt
      NbTotal = Nbs

      print*,MINVAL(AireDSommet),MAXVAL(AireDSommet)

      allocate(Nutriment(NbInc),Endothelial(NbInc),Tum(NbInc),Vasegf(NbInc),TotCells(NbInc))
      allocate(NutrimOld(NbInc),TumOld(NbInc),EndothOld(NbInc),VasegfOld(NbInc),Matieres(Nbt))
      allocate(Gb(NbInc+1:NbTotal),isbomba(NbInc),Areas(Nbt),dX(4*NbInc))
      allocate(Tum_mid(NbInc),Nut_mid(NbInc),Endo_mid(NbInc),Vegf_mid(NbInc),TotCells_mid(NbInc))

       !Structure des matrices
      call matrixinitP1sommets(Atum)
      call matrixinitP1sommets(Anut)
      call matrixinitP1sommets(Aend)
      call matrixinitP1sommets(Aveg)
      print*, 'matrixinitP1Sommets done'


      !Conditions initiales :
      Tum = 0.D0
      Nutriment = 0.D0
      Endothelial = 0.D0
      Vasegf = 0.D0
      Matieres = 0.D0

      do i=1,Nbs
         if (ABS(CoordS(1,i)-0.516052)+ABS(CoordS(2,i)-0.31312533)<1.D-3) then
            print*,'le point 1 est le numero ',i 
            print*,CoordS(1,i),CoordS(2,i)
         else if (ABS(CoordS(1,i)-0.48372)+ABS(CoordS(2,i)-0.37072)<1.D-3) then
            print*,'le point 2 est le numero ',i 
            print*,CoordS(1,i),CoordS(2,i)
         else if (ABS(CoordS(1,i)-0.62517531)+ABS(CoordS(2,i)-0.27210263)<1.D-3) then
            print*,'le point 3 est le numero ',i 
            print*,CoordS(1,i),CoordS(2,i)
         end if
      end do
      do i=1,Nbt
         Areas(i) = Ntypt(i)
         select case(Ntypt(i))
         case(100)
            Matieres(i) = Diff_u
         case(200,300,400,500)
            Matieres(i) = 5*Diff_u
         end select
      end do 

      do i=1,Nbt
         do j=1,3
         if (Ntyps(NuSoK(j,i))/=1) then 
         select case(Ntypt(i))
         case(100)
            Tum(NuSoK(j,i)) = 0.D0!cond_ini(CoordS(1,NuSoK(j,i)),CoordS(2,NuSoK(j,i)))
            call random_number(Endothelial(NuSoK(j,i)))
            Endothelial(NuSoK(j,i)) = 0.3D0 + 0.1*Endothelial(NuSoK(j,i))
         case(200)
            Tum(NuSoK(j,i)) = 0.D0!cond_ini(CoordS(1,NuSoK(j,i)),CoordS(2,NuSoK(j,i)))
            call random_number(Endothelial(NuSoK(j,i)))
            Endothelial(NuSoK(j,i)) = 1.D-3!0.3D0 + 0.1*Endothelial(NuSoK(j,i))
         case(300)
            Tum(NuSoK(j,i)) = 5.D-4
            Endothelial(NuSoK(j,i)) = 1.D-3
         case(400)
            Tum(NuSoK(j,i)) = 0.23D0!cond_ini(CoordS(1,NuSoK(j,i)),CoordS(2,NuSoK(j,i)))
            Endothelial(NuSoK(j,i)) = 1.D-3!0.299D0
         case(500)
            Tum(NuSoK(j,i)) = cond_ini(CoordS(1,NuSoK(j,i)),CoordS(2,NuSoK(j,i)))
            call random_number(Endothelial(NuSoK(j,i)))
            Endothelial(NuSoK(j,i)) = 1.D-3!0.3D0 + 0.1*Endothelial(NuSoK(j,i))
         end select
         end if
         end do
      end do

      NutrimOld = 0.D0
      CoefDiffV = Coef_diffusion;
      call tenseur(choixanis,0.D0) !Milieu isotrope
      call Newtonconditioninitiale(Anut,NutrimOld, Nutriment,Tum,Endothelial, size(NutrimOld), choixpb, i*dt)
      deallocate(Sxxk,SxyK,SyyK)

      !open(unit=4456,file='init_nutrients.txt')
      !read(4456,*) Nutriment
      !close(4456)

      VasegfOld = 0.D0
      CoefDiffV = VEGF_dif;
      call tenseur(choixanis,0.D0) !Milieu isotrope
      call NewtonVegP1sommets(Aveg,VasegfOld,Vasegf,Tum,Endothelial,Nutriment,size(NutrimOld),choixpb,i*dt)
      deallocate(Sxxk,SxyK,SyyK)

      print*,'Conditions initiales faites'

      !Créer les plots initiaux
      call plot_vtk_scheme(utmax*Tum,'TC_Pb'//strPb//'_0','Tumor_cells',ChoixSchema)
      call plot_vtk_scheme(Nutriment,'N_Pb'//strPb//'_0','Nutrients',ChoixSchema)
      call plot_vtk_scheme(Vasegf,'VEGF_Pb'//strPb//'_0','VEGF',ChoixSchema)
      call plot_vtk_scheme(utmax*Endothelial,'E_Pb'//strPb//'_0','Endothelial_cells',ChoixSchema)
      call plot_vtk(Matieres,'Matiere_Pb'//strPb,'Matieres')
      call plot_vtk(Areas,'Areas_Pb'//strPb,'Areas')
      call plot_vtk_scheme(AireDSommet,'AireD_Pb'//strPb,'Diamond_cell_area',ChoixSchema)

      open(unit=443556,file='data_Pb'//strPb//'.txt')
      write(443556,*) 0,'     ',SUM(AireDSommet*Tum,NbInc),'     ',dt,'     ',Tum(182),'     '&
      &,Nutriment(182),'     ',Endothelial(182),'     ',Vasegf(182),'     ',Tum(208),'     '&
      &,Nutriment(208),'     ',Endothelial(208),'     ',Vasegf(208),'     ',Tum(903),'     '&
      &,Nutriment(903),'     ',Endothelial(903),'     ',Vasegf(903),'     ',SUM(AireDSommet*Nutriment,NbInc)&
      &,'     ',SUM(AireDSommet*Endothelial,NbInc),'     ',SUM(AireDSommet*Vasegf,NbInc)
      print*,start_chemo,end_chemo
      !Boucle en temps 
      temps = 0.D0
      i=1
      kplot = 1
      maxresi = 0.D0
      do while(temps<=Tf)
         TumOld=Tum; NutrimOld=Nutriment; EndothOld=Endothelial; VasegfOld=Vasegf
         TotCells = TumOld + EndothOld
         choixkscalaire=choixkscalaireu
         delta=deltau; deltax=deltaxu; deltaxy=deltaxyu; deltay=deltayu;
         if (time_surg_d>=temps.and.temps+dt>=time_surg_d) then
            print*,'surgery performed'
            do i=1,Nbt
               do j=1,3 
               select case(Ntypt(i))
               case(100,200)

               case(300,400)
                  Tum(NuSoK(j,i)) = 0.D0
                  Nutriment(NuSoK(j,i)) = 0.D0
                  Endothelial(NuSoK(j,i)) = 0.D0
                  Vasegf(NuSoK(j,i)) = 0.D0
               case(500)
                  if (ntyps(NuSoK(j,i))/=6) then
                     Tum(NuSoK(j,i)) = 0.D0
                     Nutriment(NuSoK(j,i)) = 0.D0
                     Endothelial(NuSoK(j,i)) = 0.D0
                     Vasegf(NuSoK(j,i)) = 0.D0
                  end if

               end select 
               end do
            end do

            !do k=1,NbInc
            !   if (Tum(k)>=seuil_surg) then 
            !   Tum(k) = 0.D0!Tum(k)/2000000
            !   Nutriment(k) = 0.D0!Nutriment(k)/2000000
            !   Endothelial(k) = 0.D0!Endothelial(k)/2000000
            !   Vasegf(k) = 0.D0!Vasegf(k)/2000000
            !   end if
            !end do
            print*,MINVAL(Tum),MAXVAL(Tum),MINVAL(Nutriment),MAXVAL(Nutriment)
            print*,MINVAL(Endothelial),MAXVAL(Endothelial),MINVAL(Vasegf),MAXVAL(Vasegf)
         else if (calc_implicit.eqv..false.) then

            CoefTranspChi = chi_u; CoefDiffuAdeg = Diff_u; 
            call tenseur(choixanisu,temps+dt)
            call NewtonTumP1sommets(Atum,NutrimOld,TumOld,Tum,NutrimOld,TotCells,size(NutrimOld),choixpb,temps+dt)
            deallocate(Sxxk,SxyK,SyyK)
            
            CoefDiffV = Coef_diffusion;
            call tenseur(choixanis,temps+dt) !Milieu isotrope
            call NewtonNutP1sommets(Anut,NutrimOld, Nutriment,TumOld,EndothOld, size(NutrimOld), choixpb, temps+dt)
            deallocate(Sxxk,SxyK,SyyK)

            CoefTranspChi = chemo_endo; CoefDiffuAdeg = Diff_endo; 
            call tenseur(choixanisu,temps+dt)
            call NewtonEndP1sommets(Aend,VasegfOld,EndothOld,Endothelial,0.D0,TotCells,size(NutrimOld),choixpb,temps+dt)
            deallocate(Sxxk,SxyK,SyyK)

            CoefDiffV = VEGF_dif;
            call tenseur(choixanis,temps+dt) !Milieu isotrope
            call NewtonVegP1sommets(Aveg,VasegfOld,Vasegf,TumOld,EndothOld,NutrimOld,size(NutrimOld),choixpb,temps+dt)
            deallocate(Sxxk,SxyK,SyyK)

         else if (calc_implicit.eqv..true.) then

            !call NewtonAngioP1sommets(Atum,Anut,Aend,Aveg,TumOld,NutrimOld,EndothOld,VasegfOld&
            !&,Tum,Nutriment,Endothelial,Vasegf,NbInc,choixpb,i*dt,dX)
            Tum_mid = TumOld;Nut_mid=NutrimOld;Endo_mid=EndothOld;Vegf_mid=VasegfOld;
            kiter = 1

            do while (kiter<=100)

               Tum_mid=Tum; Nut_mid=Nutriment; Endo_mid=Endothelial; Vegf_mid=Vasegf
               TotCells_mid = Tum_mid + Endo_mid
               choixkscalaire=choixkscalaireu
               delta=deltau; deltax=deltaxu; deltaxy=deltaxyu; deltay=deltayu;

               CoefTranspChi = chi_u; CoefDiffuAdeg = Diff_u; 
               call tenseur(choixanisu,temps+dt)
               call NewtonTumP1sommets(Atum,Nut_mid,TumOld,Tum,Nut_mid,TotCells_mid,size(Nut_mid),choixpb,temps+dt)
               deallocate(Sxxk,SxyK,SyyK)
            
               CoefDiffV = Coef_diffusion;
               call tenseur(choixanis,temps+dt) !Milieu isotrope
               call NewtonNutP1sommets(Anut,NutrimOld, Nutriment,Tum_mid,Endo_mid, size(NutrimOld), choixpb, temps+dt)
               deallocate(Sxxk,SxyK,SyyK)

               CoefTranspChi = chemo_endo; CoefDiffuAdeg = Diff_endo; 
               call tenseur(choixanisu,temps+dt)
               call NewtonEndP1sommets(Aend,Vegf_mid,EndothOld,Endothelial,0.D0,TotCells_mid,size(Nut_mid),choixpb,temps+dt)
               deallocate(Sxxk,SxyK,SyyK)

               CoefDiffV = VEGF_dif;
               call tenseur(choixanis,temps+dt) !Milieu isotrope
               call NewtonVegP1sommets(Aveg,VasegfOld,Vasegf,Tum_mid,Endo_mid,Nut_mid,size(Nut_mid),choixpb,temps+dt)
               deallocate(Sxxk,SxyK,SyyK)

               !resi = MAX(MAXVAL(ABS(Tum-Tum_mid)),MAXVAL(ABS(Nutriment-Nut_mid)),&
               !&MAXVAL(ABS(Endothelial-Endo_mid)),MAXVAL(ABS(Vasegf-Vegf_mid)))
               !resinorm = MAX(MAXVAL(ABS(Tum_mid)),MAXVAL(ABS(Nut_mid)),MAXVAL(ABS(Endo_mid)),MAXVAL(ABS(Vegf_mid)))
               resi1 = sqrt(dot_product(Tum-Tum_mid,Tum-Tum_mid))
               resi2 = sqrt(dot_product(Nutriment-Nut_mid,Nutriment-Nut_mid))
               resi3 = sqrt(dot_product(Endothelial-Endo_mid,Endothelial-Endo_mid))
               resi4 = sqrt(dot_product(Vasegf-Vegf_mid,Vasegf-Vegf_mid))
               resinorm1 = sqrt(dot_product(Tum_mid,Tum_mid))
               resinorm2 = sqrt(dot_product(Nut_mid,Nut_mid))
               resinorm3 = sqrt(dot_product(Endo_mid,Endo_mid))
               resinorm4 = sqrt(dot_product(Vegf_mid,Vegf_mid))
               resi = MAX(resi1/resinorm1,resi2/resinorm2,resi3/resinorm3,resi4/resinorm4)

               if (kiter==100.and.resi>TolerenceImplicit) then
                  maxresi = MAX(maxresi,resi)
               end if
               
               if (kiter==1) then
                  resiold = resi
                  kiter = kiter + 1
               else if (resi>resiold) then 
                  dt = dt/5.D0
                  Tum_mid = TumOld;Nut_mid=NutrimOld;Endo_mid=EndothOld;Vegf_mid=VasegfOld;
                  !print*,temps,kiter,dt,resi
                  kiter=1
               else 
                  resiold = resi
                  kiter = kiter + 1
               end if

               if (resi <TolerenceImplicit) then
                  dt = MIN(5*dt,5.D-2)
                  exit
               end if
               
               
            end do

            
         end if
         
         
         if (AINT(temps+dt)-AINT(temps)>0) then 
            write(443556,*) temps+dt,'     ',SUM(AireDSommet*Tum,NbInc),'     ',dt,'     ',Tum(182),'     '&
            &,Nutriment(182),'     ',Endothelial(182),'     ',Vasegf(182),'     ',Tum(208),'     '&
            &,Nutriment(208),'     ',Endothelial(208),'     ',Vasegf(208),'     ',Tum(903),'     '&
            &,Nutriment(903),'     ',Endothelial(903),'     ',Vasegf(903),'     ',SUM(AireDSommet*Nutriment,NbInc)&
            &,'     ',SUM(AireDSommet*Endothelial,NbInc),'     ',SUM(AireDSommet*Vasegf,NbInc)
            str = num2string(kplot)
            call plot_vtk_scheme(Tum,'TC_Pb'//strPb//'_'//str,'Tumor_cells',ChoixSchema)
            call plot_vtk_scheme(Nutriment,'N_Pb'//strPb//'_'//str,'Nutrients',ChoixSchema)
            call plot_vtk_scheme(Vasegf,'VEGF_Pb'//strPb//'_'//str,'VEGF',ChoixSchema)
            call plot_vtk_scheme(Endothelial,'E_Pb'//strPb//'_'//str,'Endothelial_cells',ChoixSchema)
            print*,temps+dt, dt, kiter, maxresi
            kplot = kplot+1
         end if
         if (AINT(temps+dt)==radio_time(1,1)-1.and.AINT(temps)==radio_time(1,1)-2) then
            isbomba = zone_bomba(Tum)
         end if
         

         i=i+1
         temps = temps+dt
      end do
      close(443556)
   end if


100 FORMAT(10(E10.3,2x))
200 FORMAT(6(E14.6,2x))
  CLOSE (uprint,status='delete')
  PRINT*,'fin du travail laplacien'
  print*,'CHOIX PB = ', choixpb


end program KellerSegel 