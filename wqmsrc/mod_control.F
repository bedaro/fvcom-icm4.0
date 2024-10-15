!mod_control.F
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
!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|
!
Module MOD_CONTROL
  !
      Use MOD_PREC, Only: SP !
      Implicit None
      Save
  !
      Logical SERIAL !!TRUE IF SINGLE PROCESSOR
      Logical MSR !!TRUE IF MASTER PROCESSOR (MYID==1)
      Logical PAR !!TRUE IF MULTIPROCESSOR RUN
      Logical :: verbose = .False. !!Tells setup subroutines to print messages
  !
      Character (Len=80) CASENAME !!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
      Character (Len=120) CASETITLE !!CASE TITLE
  !
      Real (SP) :: HMAX !!GLOBAL MAXIMUM DEPTH
      Real (SP) :: HMIN !!GLOBAL MINIMUM DEPTH
      Real (SP) :: UMOL !!VERTICAL DIFFUSION COEFFICIENT
      Real (SP) :: HORCON !!HORIZONTAL DIFFUSION COEFFICIENT
      Real (SP) :: DTI !!internal time step
  !
      Character (Len=80) HORZMIX !!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
      Character (Len=100) :: FILENUMBER, PREFF, INPDIR, GEOAREA, &
     & RIV_FILENUMBER
      Character (Len=80) INFLOW_TYPE !!SPECIFIED RIVER INFLOW TYPE (edge/node)
      Character (Len=80) POINT_ST_TYPE !!(calculated/specified)
      Character (Len=80) PNT_SOURCE !!point_source
      Integer :: DAY, in_jday
  !
End Module MOD_CONTROL
