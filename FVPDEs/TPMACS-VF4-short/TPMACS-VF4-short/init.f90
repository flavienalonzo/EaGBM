SUBROUTINE  INIT
  !****************************************************************
  !     * Ce sous programme lit le fichier uread
  !     * le fichier d'entree
  !****************************************************************
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: i
  REAL(kind = long)    :: gbord
  CHARACTER(len=len_buffer) :: buffer
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'INIT'

  uread   = 9   ! fichier de donnee pourra etre ajoute par la suite 
  uprint = 12   ! Unite de listage pour verifier les donnees
  umesh= 39   ! unite  fichier maillage
  uele      =  15  ! unite de lecture du maillage
  uneigh    =  16  ! unite de lecture du maillage
  unode     =  17  ! unite de lecture du maillage
  uplotvtk  =  20  ! fichier de sauvgarde

  
  ! Le fichier UPRINT sert à l'ecriture de certainses données pour valider le code
  open(unit=uprint, file='UPRINT',status='unknown')

  !---------------------------------------------------------------------
  !     Donnees pour le probleme  de laplace sur un domaine qcq
  !     par  la methode des volumes finis sur un maillage admissible.
  !     equation :
  !       - c div( grad u) +theta * u + div(Vu)= F    dans omega
  !             u = Gbord        sur le bord
  !


  !-----------------------------------------------
  ! maillage triangulaire : fichier contenant le maillage  et les coonectivités

  ! il y a six maillages dispo 'MAILLAGEGEOx', x=1,...,6
  
   nom_mesh =  'MAILLAGEBRAIN'   
   !nom_mesh = 'Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt' 

  !--------------------------------------------------------------------

   iprint = 6  ! niveau d'impression


   n_enty = 4
   index_norm=1;index_nut=2;index_endo=3;index_vegf=4;
   !Pour le nutriment :
    Coef_diffusion = 3.8D-2!1.D-1
    !theta = 7.5D-1
    Coef_prod = 1.D2!5.D1
    seuil_hypo = 75.D0!6.D-2
    seuil_necro = 50.D0!3.D-2
    satur_nutri = 1.D2
    nut_degra = 3.75D-2!1.0D1
    Coef_cons = 6.D0!1.D3

   !Pour les cellules tumorales :
    Diff_u = 1.2D-6!1.D-3
   chi_u = 1D0!1.D-2
   rate = 0.27D0!1.D-1
   apop = 0.17D0!5.D-2
   rat_pop = rate - apop
    
    satur_norm = 0.D0
    choixkscalaireu = 1;
    deltau = 5.D-3; deltaxu = 3.5D-1; deltaxyu = 50.D0; deltayu = 3.5D-1;

  !Pour les cellules endotheliales
    Diff_endo = 2.16D-8!1.D-5
    chemo_endo = 1.D2!1.D-4
    satur_endo = 0.D-1
    rate_endo = 4.9D-3!5.D-4
    degr_endo = 3.1D-3!1.D-4

    !Pour le VEGF
  VEGF_prod = 3.4D2!Coef_prod
  VEGF_dif = 3.8D-3!Coef_diffusion
  VEGF_cons = 1.4!1.D1
  VEGF_degr = 15.6!nut_degra

  !Traitements
  !!Chirurgie 
  time_surg_d= 1000.0 ; time_surg_f=time_surg_d+dt; seuil_surg=5.D-4 ;
  !!Chimiotherapie
  dose_chemo = 0*1.96D-2; Upsilon=3.D0/3.D0;start_chemo=14.D0;end_chemo=56.0D0
  !!Radiotherapie
  Nradio = 6; radio_time(1,1:6) = (/14.0,21.0,28.0,35.0,42.0,49.0/); radio_time(2,1:6) = (/19.0,26.0,33.0,40.0,47.0,54.0/)
  dose_radio=2.D0 ; radio_beta=0.027D-1/30.D0 ; radio_alpha =0.027D0/30.D0 ;

  WhichPb = 66
  dt = 5.D-2
  Tf = 100
  TolerenceGradient = 1.D-13
  TolerenceNewton = 1.D-11
  TolerenceImplicit = 1.D-3
  utmax = 2.39D8
  ChoixPlot = 1 ; ChoixDegenere = 1 ; ChoixUpwind = 0 ; ChoixSchema = 2;
  choixanisu = 8 ; choixanis = 1; ChoixAdeg = 3 ; ChoixChi = 6;
  epsilon = 1.D-10
  
  !--------------------------------------------
  ! Choix du probleme à resoudre (ChoixPb)
  ! ChoixPb = 1, Uexacte =  1.
  ! ChoixPb = 2, Uexacte = x+y          
  ! ChoixPb = 3, Uexacte = x*x - y*y    
  ! ChoixPb = 4, Uexacte =  cos(5.*pi*(x+y)) 
  ! ChoixPb = 5, Uexacte =  x(1-x)y(1-y)
  ! Choixpb = 6, Uexacte = sin(pi*x)sin(pi*y)
  !----------------------------------------------
  !
  ! La fonction gbord contient la solution exacte 

  ChoixPb = 99        



  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres"

  STOP


  RETURN
END SUBROUTINE INIT








