!vdif_wqm.F
!************************************************************************
!**                                                                    **
!**                           FVCOM-ICM_4.0                            **
!**                                                                    **
!**               A Finite Volume Based Integrated Compartment         **
!**                         Water Quality Model                        **      
!**        The original unstructured-grid ICM code was developed by    ** 
!**    the FVCOM development team at the University of Massachusetts   ** 
!**         through a contract with U.S. Army Corps of Engineers       ** 
!**         [Dr. Changsheng Chen (PI), Dr. Jianhua Qi and              ** 
!**                      Dr. Geoffrey W. Cowles]                       **
!**                                                                    **
!**                Subsequent Development and Maintenance by           ** 
!**                   PNNL/UW Salish Sea Modeling Center               **
!**                                                                    **
!**                 Tarang Khangaonkar    :  PNNL (2008 - Present)     **
!**                 Lakshitha Premathilake:  PNNL (2019 - Present)     **
!**                 Adi Nugraha           :  PNNL/UW (2018 - Present)  **
!**                 Kurt Glaesmann        :  PNNL (2008 - Present)     **
!**                 Laura Bianucci        :  PNNL/DFO(2015 - Present)  **
!**                 Wen Long              :  PNNL (2012-2016)          **
!**                 Taeyum Kim            :  PNNL (2008-2011)          **
!**                 Rochelle G Labiosa    :  PNNL (2009-2010)          **
!**                                                                    **
!**                                                                    **
!**                     Adopted from CE-QUAL-ICM  Model                **
!**                           Developed by:                            **
!**                                                                    **
!**             Carl F. Cerco      : Water quality scheme              **
!**             Raymond S. Chapman : Numerical solution scheme         **
!**             Thomas M. Cole     : Computer algorithms & coding      **
!**             Hydroqual          : Sediment compartment              **
!**                                                                    **
!**                    Water Quality Modeling Group                    **
!**                    U.S. Army Corps of Engineers                    **
!**                    Waterways Experiment Station                    **
!**                    Vicksburg, Mississippi 39180                    **
!**                                                                    **
!************************************************************************
!
!Subroutine
! Subroutine VDIF_WQM()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!==============================================================================|
Subroutine VDIF_WQM (F)
  !
      Use MOD_PREC, Only: SP, DP !
      Use MOD_LIMS, Only: MTLOC, KBM1, MLOC, KB, KBM2, KBM1
  !
      Use MOD_TGE, Only: ISONB !

# if defined (WET_DRY)
	USE  MOD_WD, ONLY: ISWETN,ISWETNT
# endif
  !
      Use MOD_HYDROVARS, Only: DZ, DZ2D, DZZ, DZZ2D,D, KH
      Use MOD_SIZES, Only: NCP
      Use MOD_WQM, Only: NAC, AC, DLT, ZDFBCK, VDIFFLUX 
      Use MOD_CONTROL, Only: MSR
      Implicit None
      Real (SP), Dimension (0:MTLOC, KBM1, NCP) :: F
      Real (SP), Dimension (0:MTLOC, KBM1) :: FF
      Real (SP), Dimension (MLOC, KB) :: VHF, VHPF

      Real (SP), Dimension (MLOC, KB) :: AF, CF, RAD
  !    !!  changed to DP the following, to be as in fvcom 3jun2016
  !    Real (DP), Dimension (0:MTLOC, KBM1) :: FF
  !    Real (DP), Dimension (MLOC, KB) :: VHF, VHPF

  !    Real (DP), Dimension (MLOC, KB) :: AF, CF, RAD
!
      Real (SP), Dimension (MLOC) :: BENFLUX, WFSURF
      Real (SP), Dimension (MLOC) :: SOURCE1, SOURCE2, SOURCE3
      Real (SP), Dimension (MLOC) :: TBOT
      Real (SP) :: FKH, UMOLPR
      Real (SP) :: TEMPWUVBOT, TMP
      Real (SP) :: auxLBdiff ! 
      Integer :: I, K, J, KI, JCON, II
  !
      UMOLPR = ZDFBCK * 1.E0

      BENFLUX = 0
!
  !!VDIFFLUX = 0.0  !  initialize
  !
  !----------------------------------------------------------------
  !
  !  the following section solves the equation
  !  dti*(kh*f')' -f=-fb
  !
  !----------------------------------------------------------------
  !
      Do K = 2, KBM1
         Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
               FKH = KH (I, K)
               AF (I, K-1) = - DLT * (FKH+UMOLPR) / &
              & (DZ2D(I,K-1)*DZZ2D(I,K-1)*D(I)*D(I))
               CF (I, K) = - DLT * (FKH+UMOLPR) / &
              & (DZ2D(I,K)*DZZ2D(I,K-1)*D(I)*D(I))
            End If
         End Do
      End Do
  !
      WFSURF = 0.0
  !
  !------------------------------------------------
  !  Surface BCs; WFSURF
  !-----------------------------------------------
  !

  !  A BUNCH OF VARIABLES LOST A NCP DIMENSION
      Do JCON = 1, NAC
         II = AC (JCON)
         Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP)THEN
#  else
			IF(ISWETN(I) == 1)THEN
#  endif
               VHF (I, 1) = AF (I, 1) / (AF(I, 1)-1.)
               VHPF (I, 1) = - DLT * WFSURF (I) / (-DZ2D(I,1)*D(I)) - F (I, &
              & 1, II)
               VHPF (I, 1) = VHPF (I, 1) / (AF(I, 1)-1.)
            End If
         End Do
     !
         Do K = 2, KBM2
            Do I = 1, MLOC
#  if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
#  else
				IF(ISWETN(I) == 1)THEN
#  endif
                  VHPF (I, K) = 1. / (AF(I, K)+CF(I, K)*(1.-VHF(I, &
                 & K-1))-1.)
                  VHF (I, K) = AF (I, K) * VHPF (I, K)
                  VHPF (I, K) = (CF(I, K)*VHPF(I, K-1)-F(I, K, II)) * &
                 & VHPF (I, K)
               End If
            End Do
         End Do
     !
         Do K = 1, KBM1
            Do I = 1, MLOC
#  if !defined (WET_DRY)
				IF(D(I) > 0.0_SP)THEN
#  else
				IF(ISWETN(I) == 1)THEN
#  endif
                  FF (I, K) = F (I, K, II)
               End If
            End Do
         End Do
     !
         Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP .And. ISONB(I) /= 2)THEN
#  else
			IF(ISWETN(I) == 1 .And. ISONB(I) /= 2)THEN
#  endif
               FF (I, KBM1) = (CF(I, KBM1)*VHPF(I, KBM2)-FF(I, &
              & KBM1)-DLT*BENFLUX(I)/(D(I)*DZ2D(I,KBM1))) / (CF(I, &
              & KBM1)*(1.-VHF(I, KBM2))-1.)
               !
               auxLBdiff = (FF(I, KBM1)-F(I, KBM1, II)) / DLT * 86400.0 ! 
			   !VDIFFLUX(I, KBM1, II) = auxLBdiff ! : snapshot
               VDIFFLUX (I, KBM1, II) = VDIFFLUX (I, KBM1, II) + &
              & auxLBdiff ! : accumulated
            End If
         End Do

         Do K = 2, KBM1
            KI = KB - K
            Do I = 1, MLOC
#  if !defined (WET_DRY)
			IF(D(I) > 0.0_SP .And. ISONB(I) /= 2)THEN
#  else
			IF(ISWETN(I) == 1 .And. ISONB(I) /= 2)THEN
#  endif
                  FF (I, KI) = (VHF(I, KI)*FF(I, KI+1)+VHPF(I, KI))
              !
                  auxLBdiff = (FF(I, KI)-F(I, KI, II)) / DLT * 86400.0 ! 
				!VDIFFLUX(I, KI, II) = auxLBdiff ! : snapshot
                  VDIFFLUX (I, KI, II) = VDIFFLUX (I, KI, II) + &
                 & auxLBdiff ! : accumulated
               End If
            End Do
         End Do

         Do I = 1, MLOC
# if defined (WET_DRY)
		 IF(ISWETN(I)*ISWETNT(I)==1)THEN
#endif
			Do K = 1, KBM1
                  F (I, K, II) = FF (I, K)
             End Do
# if defined (WET_DRY)
		 ENDIF
#endif
         End Do
      End Do
  !
      Return
End Subroutine VDIF_WQM
!==============================================================================!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
