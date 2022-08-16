!mod_types.F
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
!  Converted all POINTER to ALLOCATABLE, since its faster
Module MOD_TYPES
      Use MOD_PREC, Only: SP
      Implicit None
  !
      Type GMAP
         Integer NSIZE
         Integer, Allocatable, Dimension (:) :: LOC_2_GL
      End Type GMAP
  !
      Type COMM
     !----------------------------------------------------------
     ! SND: TRUE IF YOU ARE TO SEND TO PROCESSOR               |
     ! RCV: TRUE IF YOU ARE TO RECEIVE FROM PROCESSOR          |
     ! NSND: NUMBER OF DATA TO SEND TO PROCESSOR               |
     ! NRCV: NUMBER OF DATA TO RECEIVE FROM PROCESSOR          |
     ! SNDP: ARRAY POINTING TO LOCATIONS TO SEND TO PROCESSOR  |
     ! RCVP: ARRAY POINTING TO LOCATIONS RECEIVED FROM PROCESS |
     ! RCPT: PONTER TO LOCATION IN RECEIVE BUFFER             |
     !----------------------------------------------------------
     !
     !  LOGICAL :: SND,RCV
         Integer NSND, NRCV, RCPT
         Integer, Allocatable, Dimension (:) :: SNDP, RCVP
         Real (SP), Allocatable, Dimension (:) :: MLTP
      End Type COMM
  !
      Type BC
         Integer NTIMES
         Real (SP), Allocatable, Dimension (:) :: TIMES
         Character (Len=80) :: LABEL
      End Type BC
  !
End Module MOD_TYPES
!
