!mod_lims.f
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
!==============================================================================|
!   GLOBAL LIMITS AND ARRAY SIZING PARAMETERS                                  !
!==============================================================================|
!
Module MOD_LIMS
  !
      Implicit None
      Save
  !
      Integer NLOC !!NUMBER OF ELEMENTS
      Integer MLOC !!NUMBER OF NODES
      Integer NISBCE_1 !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
      Integer NISBCE_2 !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
      Integer NISBCE_3 !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3
  !
      Integer KB !!NUMBER OF SIGMA LEVELS
      Integer KBM1 !!NUMBER OF SIGMA LEVELS-1
      Integer KBM2 !!NUMBER OF SIGMA LEVELS-2
      Integer MYID !!UNIQUE PROCESSOR ID (1 => NPROCS)
      Integer NPROCS !!NUMBER OF PROCESSORS
      Integer NE !!NUMBER OF UNIQUE EDGES
      Integer NCV !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
  ! 
      Integer IINT !!
      Integer NCV_I !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
      Integer NTLOC !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
      Integer MTLOC !!TOTAL OF LOCAL INTERNAL + HALO NODES
      Integer NCT !!(NTLOC) *3
      Integer MX_NBR_ELEM !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE
  !
      Integer NUMQBC_GL, NUMPNT_GL
      Integer NUMQBC, NUMPNT
      Integer NstationMax
      Parameter (NstationMax=200)
  !
  !Maximum number of stations
  !note this is predifined here because this is going into
  !defintion of NstationNum_GL for NAMELIST /wqm_stations/
  !fortran90 does not support dynamic arrays in namelist yet
  !
End Module MOD_LIMS
!
