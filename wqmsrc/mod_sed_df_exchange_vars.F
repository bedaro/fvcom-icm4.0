!mod_sed_df_exchange_vars.F
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
Module MOD_SED_DF_EXCHANGE_VARS
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MTLOC
      Implicit None
      Save
  !
      Real (SP) :: M1_SED_DF
      Real (SP) :: M2_SED_DF
      Real (SP), Allocatable :: POC1TM1S_SED_DF (:), POC2TM1S_SED_DF &
     & (:)!POC1TM1S and POC2TM1S shared between DF and sediment module
  !
  !
  !
Contains
  !
      Subroutine SED_DF_EXCHANGE_ALLOC
    !
         Allocate (POC1TM1S_SED_DF(MTLOC))
         POC1TM1S_SED_DF = 0.0
         Allocate (POC2TM1S_SED_DF(MTLOC))
         POC2TM1S_SED_DF = 0.0
    !
      End Subroutine SED_DF_EXCHANGE_ALLOC
  !
      Subroutine SED_DF_EXCHANGE_DEALLOC
    !
         If (ALLOCATED(POC1TM1S_SED_DF)) DEALLOCATE (POC1TM1S_SED_DF)
         If (ALLOCATED(POC2TM1S_SED_DF)) DEALLOCATE (POC2TM1S_SED_DF)
    !
      End Subroutine SED_DF_EXCHANGE_DEALLOC
End Module MOD_SED_DF_EXCHANGE_VARS
!
