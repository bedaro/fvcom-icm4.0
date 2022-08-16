!mod_sed_sav_exchange_vars.F
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
!
Module MOD_SED_SAV_EXCHANGE_VARS
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MTLOC
      Implicit None
      Save
  !
      Real (SP), Allocatable :: NH4T2TM1S_SHARE (:)!NH4T2TM1S to be shared between SAV and sediment module
      Real (SP), Allocatable :: PO4T2TM1S_SHARE (:)!PO4T2TM1S to be shared between SAV and sediment mnodule
      Real (SP), Allocatable :: HST2TM1S_SHARE (:)!
      Real (SP) :: M2_SHARE !kg/L
      Real (SP) :: PIE2HS_SHARE !partitioning coef. of H2S in sediment layer 2 (L/kg)
  !
Contains
  !
      Subroutine SED_SAV_EXCHANGE_ALLOC
    !
         Allocate (NH4T2TM1S_SHARE(MTLOC))
         NH4T2TM1S_SHARE = 0.0
         Allocate (PO4T2TM1S_SHARE(MTLOC))
         PO4T2TM1S_SHARE = 0.0
         Allocate (HST2TM1S_SHARE(MTLOC))
         HST2TM1S_SHARE = 0.0
    !
      End Subroutine SED_SAV_EXCHANGE_ALLOC
  !
      Subroutine SED_SAV_EXCHANGE_DEALLOC
    !
         If (ALLOCATED(NH4T2TM1S_SHARE)) DEALLOCATE (NH4T2TM1S_SHARE)
         If (ALLOCATED(PO4T2TM1S_SHARE)) DEALLOCATE (PO4T2TM1S_SHARE)
         If (ALLOCATED(HST2TM1S_SHARE)) DEALLOCATE (HST2TM1S_SHARE)
    !
      End Subroutine SED_SAV_EXCHANGE_DEALLOC
End Module MOD_SED_SAV_EXCHANGE_VARS
!
