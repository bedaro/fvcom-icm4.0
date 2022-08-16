!mod_sed_sf_exchange_vars.F
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
Module MOD_SED_SF_EXCHANGE_VARS
  !
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: MTLOC
      Implicit None
      Save
  !
      Real (SP), Allocatable :: JLPOC_SED_SF (:), JLPON_SED_SF (:), &
     & JLPOP_SED_SF (:), JRPOC_SED_SF (:), JRPON_SED_SF (:), &
     & JRPOP_SED_SF (:), JNH4_SED_SF (:), JPO4_SED_SF (:), SOD_SED_SF &
     & (:), JSA_SED_SF (:), JSU_SED_SF (:), SSI_SED_SF (:), SU_SED_SF &
     & (:), SA_SED_SF (:), PIP_SED_SF (:)
  !
Contains
  !
      Subroutine SED_SF_EXCHANGE_ALLOC
    !
         Allocate (JLPOC_SED_SF(MTLOC))
         JLPOC_SED_SF = 0.0
         Allocate (JLPON_SED_SF(MTLOC))
         JLPON_SED_SF = 0.0
         Allocate (JLPOP_SED_SF(MTLOC))
         JLPOP_SED_SF = 0.0
         Allocate (JRPOC_SED_SF(MTLOC))
         JRPOC_SED_SF = 0.0
         Allocate (JRPON_SED_SF(MTLOC))
         JRPON_SED_SF = 0.0
         Allocate (JRPOP_SED_SF(MTLOC))
         JRPOP_SED_SF = 0.0
    !
         Allocate (JNH4_SED_SF(MTLOC))
         JNH4_SED_SF = 0.0
         Allocate (JPO4_SED_SF(MTLOC))
         JPO4_SED_SF = 0.0
         Allocate (SOD_SED_SF(MTLOC))
         SOD_SED_SF = 0.0
         Allocate (JSA_SED_SF(MTLOC))
         JSA_SED_SF = 0.0
         Allocate (JSU_SED_SF(MTLOC))
         JSU_SED_SF = 0.0
         Allocate (SSI_SED_SF(MTLOC))
         SSI_SED_SF = 0.0
    !
         Allocate (SU_SED_SF(MTLOC))
         SU_SED_SF = 0.0
         Allocate (SA_SED_SF(MTLOC))
         SA_SED_SF = 0.0
         Allocate (PIP_SED_SF(MTLOC))
         PIP_SED_SF = 0.0
    !
    !
    !
      End Subroutine SED_SF_EXCHANGE_ALLOC
  !
      Subroutine SED_SF_EXCHANGE_DEALLOC
    !
         If (ALLOCATED(JLPOC_SED_SF)) DEALLOCATE (JLPOC_SED_SF)
         If (ALLOCATED(JLPON_SED_SF)) DEALLOCATE (JLPON_SED_SF)
         If (ALLOCATED(JLPOP_SED_SF)) DEALLOCATE (JLPOP_SED_SF)
         If (ALLOCATED(JRPOC_SED_SF)) DEALLOCATE (JRPOC_SED_SF)
         If (ALLOCATED(JRPON_SED_SF)) DEALLOCATE (JRPON_SED_SF)
         If (ALLOCATED(JRPOP_SED_SF)) DEALLOCATE (JRPOP_SED_SF)
    !
         If (ALLOCATED(JNH4_SED_SF)) DEALLOCATE (JNH4_SED_SF)
         If (ALLOCATED(JPO4_SED_SF)) DEALLOCATE (JPO4_SED_SF)
         If (ALLOCATED(SOD_SED_SF)) DEALLOCATE (SOD_SED_SF)
         If (ALLOCATED(JSA_SED_SF)) DEALLOCATE (JSA_SED_SF)
         If (ALLOCATED(JSU_SED_SF)) DEALLOCATE (JSU_SED_SF)
         If (ALLOCATED(SSI_SED_SF)) DEALLOCATE (SSI_SED_SF)
         If (ALLOCATED(SU_SED_SF)) DEALLOCATE (SU_SED_SF)
         If (ALLOCATED(SA_SED_SF)) DEALLOCATE (SA_SED_SF)
         If (ALLOCATED(PIP_SED_SF)) DEALLOCATE (PIP_SED_SF)
    !
      End Subroutine SED_SF_EXCHANGE_DEALLOC
End Module MOD_SED_SF_EXCHANGE_VARS
!
