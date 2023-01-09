!            **************************
!            **  SUBROUTINE MATRIXINIT**
!            **************************
!******************************************************************************

!     * Ce sous programme donne la structure generale de 
!     * de la matrice A%TMAT ainsi le remplissage de 
!     * A%IndPL  
!     * A%Indc
!     * la matrice elementaire locale Aloc du triangle jt
!     * Ajouter a la ligne I et la colone J la valeur Coef

!****************************************************************************
SUBROUTINE matrixinitP1Sommets(Mat)
  !--------
  ! Modules
  !--------
  USE longr
  USE parmmage
  USE imprime
  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux),         INTENT(OUT) :: Mat
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER                             :: iseg, is, js, i1, jv, kv,jmin,mloc
  INTEGER                             :: NcoefMat
  INTEGER, DIMENSION(NbInc)           :: IndPL

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MATINIT'
  !------
  NcoefMat = NbInc!NsInt 

  DO iseg = 1, Nseg
     is = NuSeg(1,iseg) ;   js = NuSeg(2,iseg)
     
     if(is >NbInc .or. js> NbInc) cycle
     !IF(is > NsInt .OR. js > NsInt) CYCLE
     NcoefMat = NcoefMat + 2 
     
  END DO
  !print*,NcoefMat
  !ALLOCATE(Mat%IndPL(1:NsInt+1), Mat%Indc(NcoefMat), Mat%TMat(NcoefMat))
  allocate(Mat%IndPL(1:NbInc+1), Mat%Indc(NcoefMat),Mat%TMat(NcoefMat))
  !ALLOCATE( Mat%F ( NsInt ), Mat%Bg ( NsInt ) )
  allocate(Mat%F(NbInc),Mat%Bg(NbInc))
  !write(*,*)'NcoefMat=',NcoefMat
  
  Mat%TMat = 0.0D0 ;  Mat%F = 0.D0 ;  Mat%Bg = 0.D0
  DO is = 1, NbInc+1!NsInt+1
     Mat%IndPL(is) = is
  END DO
  
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg) ; js = NuSeg(2, iseg)
     !write(*,*)is,js
     
     !IF(is > NsInt .OR. js > NsInt) CYCLE
     if(is >NbInc .or. js>NbInc) cycle

     !Mat%IndPL(is+1:NsInt+1) = Mat%IndPL(is+1:NsInt+1) + 1
     Mat%IndPL(is+1:NbInc+1) = Mat%IndPL(is+1:NbInc+1)+1
     !Mat%IndPL(js+1:NsInt+1) = Mat%IndPL(js+1:NsInt+1) + 1
     Mat%IndPL(js+1:NbInc+1) = Mat%IndPL(js+1:NbInc+1) + 1
     
  END DO
  
  !IndPL(1:NsInt) = Mat%IndPL(1:NsInt)
  IndPL(1:NbInc) = Mat%IndPL(1:NbInc)

  DO is = 1, NbInc!NsInt
     Mat%Indc( IndPL(is) ) = is
     IndPL(is)           = IndPL(is) + 1
  ENDDO
  DO iseg = 1, Nseg
   is = NuSeg(1,iseg) ;  js = NuSeg(2, iseg)

     !IF (is > NsInt .OR. js > NsInt) CYCLE
   if (is>NbInc.or.js>NbInc) cycle
     Mat%Indc(IndPL(is)) = js
     Mat%Indc(IndPL(js)) = is
     IndPL(is)           = IndPL(is) + 1
     IndPL(js)           = IndPL(js) + 1

  END DO

  DO is = 1, NbInc!NsInt 
     DO jv = Mat%IndPL(is), Mat%IndPL(is+1)-1 
        jmin = jv
        do   kv = jv +1,  Mat%IndPL(is+1)-1
           IF (Mat%Indc(kv) < Mat%Indc(jmin)) jmin = kv
        END DO
        mloc =  Mat%Indc(jmin )
        Mat%Indc(jmin ) = Mat%Indc(jv ) 
        Mat%Indc(jv ) = mloc 
     END DO
  END DO
  !write(*,*)'A%IndPL=',Mat%IndPL
  !write(*,*)'A%Indc=',Mat%Indc

  !------------
  ! Impressions
  !------------
  !IF (iprint >=5) THEN
   !  CALL prvari(uprint,'NcoefMat = ', NcoefMat )
   !  WRITE(uprint,*) ( Mat%Tmat(is) , is=1, NcoefMat )
   !  WRITE(uprint,*) ( Mat%Indc(is) , is=1, NcoefMat )
   !  WRITE(uprint,*) ( Mat%IndPL(is), is=1, NsInt+1 )
  !ENDIF
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE matrixinitP1Sommets



