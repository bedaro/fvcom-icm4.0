!mod_ba.F
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
!************************************************************************
!**                                                                    **
!** Benthic Algage Module                                              **
!**                                                                    **
!**                                                                    **
!************************************************************************
!**                                                                    **
!** Inputs:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!** Outpus:                                                            **
!**                                                                    **
!**            Required inputs for benthic algae model                 **
!**                                                                    **
!************************************************************************
!
Module MOD_BA
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MTLOC, MLOC
      Use MOD_WQM, Only: AOCR, BAOFN, BB, BENDO, BENDOC, BENNH4, &
     & BENNO3, BENPO4, DLT, DOXG, KADPO4, NH4, NO3, PO4, SSI, T
!
      Use MOD_FILEINFO, Only: BAI, BAO
      Use MOD_HYDROVARS, Only: DZ,DZ2D, D
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR
  !
      Use MOD_OWQ, Only: IATBOT
  !
      Implicit None
      Save
  !
      Real (SP), Allocatable :: DTPOC1_BA (:), DTPOC2_BA (:), DTPOC3_BA &
     & (:), DTPON1_BA (:), DTPON2_BA (:), DTPON3_BA (:), DTPOP1_BA (:), &
     & DTPOP2_BA (:), DTPOP3_BA (:)
  !
      Real (SP), Allocatable :: DTAG3CFL_BA (:), DTAG3NFL_BA (:), &
     & DTAG3PFL_BA (:)!change of accumulated G3 Phosphorus flux
!
  !	DTAG3CFL_BA(:) !change of accumulated G3 Carbon flux during a sediment time step dlt_ba
  !   DTAG3NFL_BA(:) !change of accumulated G3 Nitrogen flux during a sediment time step dlt_ba
  !   DTAG3PFL_BA(:) !change of accumulated G3 Phosphorus flux during a sediment time step dlt_ba
!
      Real (SP), Allocatable :: HBA (:)!thickness of sediment layer that grows benthic algae (H1+H2 (m))
!
  !
  !variables that are used for as flags of model simulation
  !
  !
  !
      Logical :: BALGAE_CALC !flag for activating benthic algae
  !calculation
  !
  !
  !
  ! variables that are local to benthic algae module
  !
      Real (SP) :: DLT_BA !benthic algae biomass simulation time step (day)
  !==DLT/86400, where DLT is in seconds
  !
      Real (SP), Dimension (3) :: FRCPHB, FRNPHB, FRPPHB
  !FRCPHB    !Fraction of POC generated for G1,G2,G3 due to predation on benthic algae
  !FRNPHB    !Fraction of PON genereated for G1,G2,G3 due to predation on benthic algae
  !FRPPHB    !Fraction of POP generated for G1,G2,G3 due to predation on benthic algae
  !
  !Moved from subroutine SED_READ
      Character (Len=3) :: BALC
  !
      Real (SP) :: HBAALL !depth of sediments with benthic algae (m)
  !
      Integer :: INTBAC !steady state flag of bentic algae simulation
  !(1 for steady state simulation, 0 for time varying)
      Logical :: STEADY_STATE_BA !flag for steady state simulation of benthic algae
  !
      Logical :: BA_INITIALIZED !flag for if benthic algae model is initialized
  !
  !benthic algae growth kinetics
  !Moved from wqm_modules.F
      Real (SP) :: PMB, ANCB, APCB, KTGB1, KTGB2, TMB, KESED, KEBALG, &
     & KHNB, KHPB, KHRB, BMRB, BPRB, KTBB, TRB, BALGMIN, FNIB, FPIB
  !
  !
      Real (SP) :: ALPHB, CCHLB !C/CHLA ratio of benthic algae gC/gCHLA
  !
      Logical :: BA_OUTPUT !Flag for benthic algae output
  !
  !Checked the units of these settling velocities
  !	REAL(SP),ALLOCATABLE,DIMENSION(:) ::! WSSBNET,  & !netsettling rate of inorganic suspended sediments (m/d)
  ! WSLBNET,  & !net settling rate of LPOM (m/d)
  ! WSRBNET,  & !net settling rate of RPOM (m/d)
  ! WS1BNET,  & !net settling rate of alg 1 (m/d) (normal settling discounted by resuspension rate)
  ! WS2BNET,  & !net settling rate of alg 2 (m/d)
       ! WS3BNET,  & !net settling rate of alg 3 (m/d)
  ! WSUBNET     !net settling rate of particulate biogenic silica (m/d)
  !
  !Moved from subroutine SED_READ
  !	CHARACTER(LEN=8)  ::	SPVARB, &
  !							PRINTB
  !
  !
  !Moved from subroutine SED_READ
  !	CHARACTER(LEN=8)  :: SPVARLR,  &
  !						 PRINTLR
!
  !
  !other local variables that are used in the benthic algae module
  !
!
      Integer :: KWC
!
      Character (Len=20) :: SSNAME_BA (17)
  !
      Real (SP) :: NH4AVL, NO3AVL !,   &
!
      Real (SP) :: DF, PO4AVL
  !available PO4 for benthic algae growth (=DF*PO4) in bottom layer of water column (gP/L ?)
!
      Real (SP) :: IK
  !
      Real (SP) :: FTB, PRNB, FRDOB !
!
  !Benthic algae growth limitation and kinetic fluxes
  !Moved this to wqm_sed.F.
  !benthic algae nutrient limitation and kinetic fluxes
  !
      Real (SP), Allocatable, Dimension (:) :: FIB, NLB, PLB, NPPB, &
     & BBM, BLITE, BMB, PB, PRB, BANH4, BANO3, BAPO4, BADOC, BADO, &
     & BAPOC, BAPON, BAPOP !adjustment to benthic POP flux by benthic algae (gP/m^2/day)
  !
  !Moved here from wqm_modules.F
      Real (SP) :: TINTIM_BA !time (days) for steady state simulation of sediments
  !
  !
      Real (SP), Allocatable :: BBM_GL (:)
  !
  !
Contains
  !
  !subroutines:
  !subroutine BA_INIT()
  !subroutine BA_ALLOC()
  !subroutine BA_DEALLOC()
  !subroutine BA_READ()
  !subroutine BA_CALC()
  !subroutine BA_INT()
  !
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   BA_INIT                           **
  !********************************************************************************
      Subroutine BA_INIT
    !
         Integer :: I
    !Moved from wqm_modules.F
    !
         SSNAME_BA (15) = 'Benthic Algae       '
    !
    !allocate variables related to sediment diagenesis
         Call BA_ALLOC
    !
    !Set initial benthic algae sediment depth
         Do I = 1, MLOC
            HBA (I) = HBAALL * 0.01 !Convertion from cm to meter
         End Do
    !
         BA_INITIALIZED = .True.
    !
      End Subroutine BA_INIT
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   SED_ALLOC                         **
  !********************************************************************************
      Subroutine BA_ALLOC
    !
    !Moved this here from wqm_modules.F
         Allocate (FIB(MTLOC))
         FIB = 0.0
         Allocate (NLB(MTLOC))
         NLB = 0.0
         Allocate (PLB(MTLOC))
         PLB = 0.0
         Allocate (NPPB(MTLOC))
         NPPB = 0.0
         Allocate (BBM(MTLOC))
         BBM = 0.0
         Allocate (BLITE(MTLOC))
         BLITE = 0.0
         Allocate (BMB(MTLOC))
         BMB = 0.0
         Allocate (PB(MTLOC))
         PB = 0.0
         Allocate (PRB(MTLOC))
         PRB = 0.0
         Allocate (BANH4(MTLOC))
         BANH4 = 0.0
         Allocate (BANO3(MTLOC))
         BANO3 = 0.0
         Allocate (BAPO4(MTLOC))
         BAPO4 = 0.0
         Allocate (BADOC(MTLOC))
         BADOC = 0.0
         Allocate (BADO(MTLOC))
         BADO = 0.0
         Allocate (BAPOC(MTLOC))
         BAPOC = 0.0
         Allocate (BAPON(MTLOC))
         BAPON = 0.0
         Allocate (BAPOP(MTLOC))
         BAPOP = 0.0
    !
         Allocate (DTPOC1_BA(MTLOC))
         DTPOC1_BA = 0.0
         Allocate (DTPOC2_BA(MTLOC))
         DTPOC2_BA = 0.0
         Allocate (DTPOC3_BA(MTLOC))
         DTPOC3_BA = 0.0
         Allocate (DTPON1_BA(MTLOC))
         DTPON1_BA = 0.0
         Allocate (DTPON2_BA(MTLOC))
         DTPON2_BA = 0.0
         Allocate (DTPON3_BA(MTLOC))
         DTPON3_BA = 0.0
         Allocate (DTPOP1_BA(MTLOC))
         DTPOP1_BA = 0.0
         Allocate (DTPOP2_BA(MTLOC))
         DTPOP2_BA = 0.0
         Allocate (DTPOP3_BA(MTLOC))
         DTPOP3_BA = 0.0
    !
         Allocate (DTAG3CFL_BA(MTLOC))
         DTAG3CFL_BA = 0.0
         Allocate (DTAG3NFL_BA(MTLOC))
         DTAG3NFL_BA = 0.0
         Allocate (DTAG3PFL_BA(MTLOC))
         DTAG3PFL_BA = 0.0
    !
         Allocate (HBA(MTLOC))
         HBA = 0.0
    !
    !
      End Subroutine BA_ALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   BA_DEALLOC                        **
  !********************************************************************************
  !
      Subroutine BA_DEALLOC
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(FIB)) DEALLOCATE (FIB)
         If (ALLOCATED(NLB)) DEALLOCATE (NLB)
         If (ALLOCATED(PLB)) DEALLOCATE (PLB)
         If (ALLOCATED(NPPB)) DEALLOCATE (NPPB)
         If (ALLOCATED(BBM)) DEALLOCATE (BBM)
         If (ALLOCATED(BLITE)) DEALLOCATE (BLITE)
         If (ALLOCATED(BMB)) DEALLOCATE (BMB)
         If (ALLOCATED(PB)) DEALLOCATE (PB)
         If (ALLOCATED(PRB)) DEALLOCATE (PRB)
         If (ALLOCATED(BANH4)) DEALLOCATE (BANH4)
         If (ALLOCATED(BANO3)) DEALLOCATE (BANO3)
         If (ALLOCATED(BAPO4)) DEALLOCATE (BAPO4)
         If (ALLOCATED(BADOC)) DEALLOCATE (BADOC)
         If (ALLOCATED(BADO)) DEALLOCATE (BADO)
         If (ALLOCATED(BAPOC)) DEALLOCATE (BAPOC)
         If (ALLOCATED(BAPON)) DEALLOCATE (BAPON)
         If (ALLOCATED(BAPOP)) DEALLOCATE (BAPOP)
    !
         If (ALLOCATED(DTPOC1_BA)) DEALLOCATE (DTPOC1_BA)
         If (ALLOCATED(DTPOC2_BA)) DEALLOCATE (DTPOC2_BA)
         If (ALLOCATED(DTPOC3_BA)) DEALLOCATE (DTPOC3_BA)
         If (ALLOCATED(DTPON1_BA)) DEALLOCATE (DTPON1_BA)
         If (ALLOCATED(DTPON2_BA)) DEALLOCATE (DTPON2_BA)
         If (ALLOCATED(DTPON3_BA)) DEALLOCATE (DTPON3_BA)
         If (ALLOCATED(DTPOP1_BA)) DEALLOCATE (DTPOP1_BA)
         If (ALLOCATED(DTPOP2_BA)) DEALLOCATE (DTPOP2_BA)
         If (ALLOCATED(DTPOP3_BA)) DEALLOCATE (DTPOP3_BA)
!
         If (ALLOCATED(DTAG3CFL_BA)) DEALLOCATE (DTAG3CFL_BA)
         If (ALLOCATED(DTAG3NFL_BA)) DEALLOCATE (DTAG3NFL_BA)
         If (ALLOCATED(DTAG3PFL_BA)) DEALLOCATE (DTAG3PFL_BA)
!
         If (ALLOCATED(HBA)) DEALLOCATE (HBA)
!
         If (ALLOCATED(BBM_GL)) DEALLOCATE (BBM_GL)
!
      End Subroutine BA_DEALLOC
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   B A _ R E A D                     **
  !********************************************************************************
      Subroutine BA_READ
#  if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NGID, NHN, HN_LST
#  endif
         Implicit None
         Save
!
    !***** Variable declarations
    !
         Character (Len=24) :: FRNAME_BA (3)
         Integer :: I, JG ! , JT
    !
         Data FRNAME_BA / 'Benthic algal carbon    ', 'Benthic algal ni&
        &trogen  ', 'Benthic algal phosphorus' /!3   Benthic algae P
    !
    !
    !open file for diagnostics on input format check
    !
         If (BA_OUTPUT) Then
            If (MSR) Then
               Open (BAO, File=BAOFN)
            End If
         End If
    !********************************************************************************
    !**                                  Inputs                                    **
    !********************************************************************************
    !
    !
    ! benthic algae
         Read (BAI, 1000, Err=10100) HBAALL, INTBAC !depth of sediments with benthic algae,
    !steady state flag
         Read (BAI, 1015, Err=10100) PMB, ANCB, APCB, KTGB1, KTGB2, TMB
         Read (BAI, 1015, Err=10100) ALPHB, CCHLB, KESED, KEBALG, KHNB, &
        & KHPB, KHRB
         Read (BAI, 1015, Err=10100) BMRB, BPRB, KTBB, TRB, BALGMIN
         Read (BAI, 1015, Err=10100) FNIB, FPIB
         Read (BAI, 1010, Err=10100) FRPPHB !Fractions of P in benthic algae predation to POM for G1, G2, G3
         Read (BAI, 1010, Err=10100) FRNPHB !Fractions of N in benthic algae predation to POM for G1, G2, G3
         Read (BAI, 1010, Err=10100) FRCPHB !Fractions of C in benthic algae predation to POM for G1, G2, G3
    !
         Close (BAI)
    !
         STEADY_STATE_BA = INTBAC == 1
    !
         BA_INITIALIZED = .False.
    !
    !********************************************************************************
    !**                                 Outputs                                    **
    !********************************************************************************
    !
         If (BA_OUTPUT) Then
            If (MSR) Then
               Write (BAO, 2020) HBAALL
               If (STEADY_STATE_BA) Then
                  Write (BAO, 2022)
               Else
                  Write (BAO, 2025)
               End If
               Write (BAO, 2030)
               Write (BAO, 2110) FRNAME_BA (3), FRPPHB
               Write (BAO, 2110) FRNAME_BA (2), FRNPHB
               Write (BAO, 2110) FRNAME_BA (1), FRCPHB
          !
               Write (BAO, 2252) BALC
               Write (BAO, 2242) PMB, ANCB, APCB, KTGB1, KTGB2, TMB
               Write (BAO, 2244) ALPHB, CCHLB, KESED, KEBALG, KHNB, &
              & KHPB, KHRB
               Write (BAO, 2246) BMRB, BPRB, KTBB, TRB, BALGMIN
               Write (BAO, 2248) FNIB, FPIB
          !
               Close (BAO)
            End If
         End If
    !
    !
    !***** Input FORMAT'S
    !
1000     Format (:////8 X, F8.0, I8)
1010     Format (8 F10.0)
    !
1015     Format (/ / 8 X, 8 F8.1)
    !
    !***** Output FORMAT'S
2020     Format (/ ' ACTIVE LAYER DEPTH ', F8.3, ' CM')
2022     Format (/ ' STEADY-STATE VALUES OF BENTHIC ALGAE  COMPUTED' /)
2025     Format (/ ' NON STEADY-STATE VALUES OF BENTHIC ALGAE COMPUTED' &
        & /)
    !
2030     Format (/ // / 33 X, ' BENTHIC ALGAE    I N I T I A L   C O N &
        &D I T ', 'I O N S' /)
2110     Format (6 X, A24, 11 X, 3 F7.2)
    !
2242     Format (' PMB = ', F10.3/' ANCB = ', F10.3/' APCB = ', F10.3/'&
        & KTGB1 = ', F10.3/' KTGB2 = ', F10.3/' TMB = ', F10.3)
2244     Format (' ALPHB = ', F10.3/' CCHLB = ', F10.3/' KESED = ', &
        & F10.3/' KEBALG = ', F10.3, ' KHNB = ', F10.3/' KHPB = ', &
        & F10.3/' KHRB = ', F10.3)
2246     Format (' BMRB = ', F10.3/' BPRB = ', F10.3/' KTBB = ', &
        & F10.3/' TRB = ', F10.3/' BALGMIN = ', F10.3)
2248     Format (' FNIB = ', F10.3/' FPIB = ', F10.3)
2252     Format (/ ' BENTHIC ALGAE ARE ', A3)
         Return
    !
    !
10100    Continue
    !
         If (MSR) Then
            If (BA_OUTPUT) Then
               Write (BAO, 3010)
               Close (BAO)
            End If
         End If
3010     Format (/ ' Read error in benthic algae input deck')
         Stop 'BA_READ'
    !
         Return
      End Subroutine BA_READ
  !
  !
  !********************************************************************************
  !**                    S U B R O U T I N E   B A _ C A L C                     **
  !**                           Benthic Algae Calculations                       **
  !********************************************************************************
  !
      Subroutine BA_CALC
    !
         Implicit None
         Save
         Integer :: I

         If ( .Not. BA_INITIALIZED) Then
            Call BA_INIT ()
         End If
    !***** Pass MOD_WQM time-step (in days) to sediment subr
    !
         DLT_BA = DLT / 86400.
         If (STEADY_STATE_BA) TINTIM_BA = TINTIM_BA + DLT_BA

    !***** Update sediment concentrations
    !
         Do I = 1, MLOC

            If (BALGAE_CALC) Then
          !
          !******* Benthic algae algorithms start here
          !
          !******* Calculate mean light in algal mat
          !
               BLITE (I) = IATBOT (I, KWC) * Exp (-KESED) / &
              & (KEBALG+1.0E-8) / BBM (I) * &
              & (1.-Exp(-(KEBALG+1.0E-8)*BBM(I)))
          !
          !******* Temperature effects
          !
               If (T(I, KWC) < TMB) Then
                  FTB = Exp (-KTGB1*(T(I, KWC)-TMB)**2)
               Else
                  FTB = Exp (-KTGB2*(TMB-T(I, KWC))**2)
               End If
          !
          !******* Light effects
          !
               IK = PMB * FTB / ALPHB

          !******* Nutrient limitations
          !
          ! COMPUTE AVAILABLE AMMONIUM AND NITRATE
          !
               NH4AVL = BENNH4 (I) * DLT_BA + NH4 (I, KWC) * D (I) * DZ2D &
              & (I,KWC)
               NH4AVL = Max (0., NH4AVL)
               NO3AVL = BENNO3 (I) * DLT_BA + NO3 (I, KWC) * D (I) * DZ2D &
              & (I,KWC)
               NO3AVL = Max (0., NO3AVL)
          !
          ! COMPUTE NITROGEN LIMITATION
          !
               NLB (I) = (NH4AVL+NO3AVL) / (KHNB+NH4AVL+NO3AVL)
          !
          ! COMPUTE NITROGEN PREFERENCE
          !
               PRNB = NH4AVL * NO3AVL / ((KHNB+NH4AVL)*(KHNB+NO3AVL)) + &
              & NH4AVL * KHNB / ((1.E-30+NH4AVL+NO3AVL)*(KHNB+NO3AVL))
          !
          !******* Phosphorus available for algal growth
          !
               DF = 1. / (1.+KADPO4*SSI(I, KWC))
               PO4AVL = DF * PO4 (I, KWC) * D (I) * DZ2D (I,KWC)
               PO4AVL = PO4AVL + BENPO4 (I) * DLT_BA
               PO4AVL = Max (0., PO4AVL)
               PLB (I) = PO4AVL / (KHPB+PO4AVL)
          !
          !******* Base metabolism
          !
          ! IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET METABOLISM TO ZERO
          !
               If (BBM(I) > BALGMIN) Then
                  BMB (I) = BMRB * Exp (KTBB*(T(I, KWC)-TRB))
               Else
                  BMB (I) = 0.
               End If
          !
          !******* Production
          !
               PB (I) = PMB * FTB * AMIN1 (FIB(I), NLB(I), PLB(I)) / &
              & CCHLB
          !
          !******* Net primary production
          !
               NPPB (I) = (PB(I)-BMB(I)) * BBM (I)
          !
          !******* Predation
          !
          ! IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET PREDATION TO ZERO
          !
               If (BBM(I) > BALGMIN) Then
                  PRB (I) = BBM (I) * BPRB * Exp (KTBB*(T(I, KWC)-TRB))!
             !for PRB should have units 1/d
               Else
                  PRB (I) = 0.
               End If
          !
          ! ADJUST PREDATION SO BIOMASS DOESN'T GO NEGATIVE
          !
               PRB (I) = Min (PRB(I), PB(I)-BMB(I)+0.99/DLT_BA)
          !
          !******* Compute effects of algal activity on benthic flux
          !
               BANH4 (I) = ANCB * (BMB(I)*FNIB-PRNB*PB(I)+PRB(I)*FNIB) &
              & * BBM (I)
               BANO3 (I) = - (1.-PRNB) * PB (I) * ANCB * BBM (I)
               BAPO4 (I) = APCB * (BMB(I)*FPIB-PB(I)+PRB(I)*FPIB) * BBM &
              & (I)
               FRDOB = 1. - KHRB / (DOXG(I, KWC)+KHRB)
               BADO (I) = ((1.3-0.3*PRNB)*PB(I)-FRDOB*BMB(I)) * AOCR * &
              & BBM (I)
               BADOC (I) = (1.-FRDOB) * BMB (I) * BBM (I)
          !
          !stuff to pass to water column due to change to NH4, NO3, PO4, DOC, DOXG in water column
               BENNH4 (I) = BENNH4 (I) + BANH4 (I)!gN/m^2/day
               BENNO3 (I) = BENNO3 (I) + BANO3 (I)!gN/m^2/day
               BENPO4 (I) = BENPO4 (I) + BAPO4 (I)!gP/m^2/day
               BENDOC (I) = BENDOC (I) + BADOC (I)!gC/m^2/day !LDOC
               BENDO (I) = BENDO (I) + BADO (I)!gO2/m^2/day
          !
          ! COMPUTE EFFECTS OF ALGAL ACTIVITY ON ORGANIC PARTICULATES (MG/M**3)
          !stuff to pass to sediment module
               BAPOC (I) = PRB (I) * BBM (I)
               BAPON (I) = ANCB * (1.-FNIB) * (BMB(I)+PRB(I)) * BBM (I)
               BAPOP (I) = APCB * (1.-FPIB) * (BMB(I)+PRB(I)) * BBM (I)
          !
               DTPOC1_BA (I) = 1000. * BAPOC (I) * FRCPHB (1) * DLT_BA &
              & / HBA (I)!1000*(gC/m^2/day) * (day) /(m) ==> mgC/m^3
               DTPOC2_BA (I) = 1000. * BAPOC (I) * FRCPHB (2) * DLT_BA &
              & / HBA (I)
               DTPOC3_BA (I) = 1000. * BAPOC (I) * FRCPHB (3) * DLT_BA &
              & / HBA (I)
               DTPON1_BA (I) = 1000. * BAPON (I) * FRNPHB (1) * DLT_BA &
              & / HBA (I)
               DTPON2_BA (I) = 1000. * BAPON (I) * FRNPHB (2) * DLT_BA &
              & / HBA (I)
               DTPON3_BA (I) = 1000. * BAPON (I) * FRNPHB (3) * DLT_BA &
              & / HBA (I)
               DTPOP1_BA (I) = 1000. * BAPOP (I) * FRPPHB (1) * DLT_BA &
              & / HBA (I)
               DTPOP2_BA (I) = 1000. * BAPOP (I) * FRPPHB (2) * DLT_BA &
              & / HBA (I)
               DTPOP3_BA (I) = 1000. * BAPOP (I) * FRPPHB (3) * DLT_BA &
              & / HBA (I)
          !
          !******* Accumulate fluxes for steady-state computation
          !
               If (STEADY_STATE_BA) Then
             !

                  DTAG3CFL_BA (I) = 1000. * PRB (I) * FRCPHB (3) * BBM &
                 & (I) * DLT_BA
                  DTAG3NFL_BA (I) = 1000. * PRB (I) * FRNPHB (3) * ANCB &
                 & * BBM (I) * DLT_BA
                  DTAG3PFL_BA (I) = 1000. * PRB (I) * FRPPHB (3) * APCB &
                 & * BBM (I) * DLT_BA
               End If
          !
          !******* Change in benthic algal biomass
          !
               BBM (I) = BBM (I) * (1.+DLT_BA*(PB(I)-BMB(I)-PRB(I)))

            End If

         End Do

         Return
      End Subroutine BA_CALC
  !
      Subroutine BA_INT ()
    !
         Integer :: I
    !
         Return
      End Subroutine BA_INT

End Module MOD_BA
