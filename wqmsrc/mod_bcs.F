!mod_bcs.F
!!************************************************************************
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
Module MOD_BCS
  !
      Use MOD_TYPES, Only: BC
      Use MOD_PREC, Only: SP
  !
      Use MOD_CONTROL, Only: PAR
  !
      Implicit None
      Save
  !----------------boundary conditions: fresh water discharge-------------------------!
      Integer, Allocatable :: INODEQ (:)!!LOCAL FRESH WATER INFLOW NODES
      Integer, Allocatable :: ICELLQ (:)!!LOCAL FRESH WATER INFLOW ELEMENTS
      Integer, Allocatable :: RIV_GL2LOC (:)!!GLOBAL TO LOCAL MAP OF FW NODES
      Integer, Allocatable :: INOPNT (:)!!LOCAL NON-POINT SOURCE NODES
  !
      Type (BC) :: QBC_TM !!TIME MAP FOR RIVER DISCHARGE DATA
      Type (BC) :: PNT_TM !!TIME MAP FOR NON POINT SOURCE DATA
  !
  !Moved this to mod_bcmap.F
  !TYPE(BC)              :: NUT_TM           !!TIME MAPPING FOR NUTRIENT OBC  !should be
  !
  !These are never used for we are not solving hydrodynamics
  !REAL(SP), ALLOCATABLE     ::  VQDIST(:,:)     !!DISCHARGE VERTICAL DISTRIBUTION
  !REAL(SP), ALLOCATABLE     ::   DQDIS(:,:)     !!WATER FLUX DISCHARGE DATA
  !
      Real (SP), Allocatable :: DWDIS (:, :, :)!!WATER QUALITY DISCHARGE DATA
      Real (SP), Allocatable :: WDIS (:, :)!!CURRENT TIME FRESH WATER QUALITY
      Real (SP), Allocatable :: WVQDIST (:, :)!!DISCHARGE VERTICAL DISTRIBUTION for point source
  !
  !
      Real (SP), Allocatable :: PQDIS (:)!!Current DISCHARGE at non-point source node
      Real (SP), Allocatable :: PDQDIS (:, :)!!DISCHARGE at non-point source node
      Real (SP), Allocatable :: PWQDIS (:, :)!!Current WATER QUALITY at non-point source node
  !
      Real (SP), Allocatable :: PDWQDIS (:, :, :)!!WATER QUALITY DATA at non-point source node
  !
  !Moved this to mod_bcmap.F
  !REAL(SP), ALLOCATABLE     :: WQOBC(:,:,:,:)   !!WATER QUALITY DATA AT BOUNDARY	!
  !===================================================================================!
  !
      Integer, Allocatable :: N_ICELLQ (:, :)!!node number of the edge for element edge source
  !
Contains
  !
  !Subroutine BCS_ALLOC()
  !Subroutine BCS_DEALLOC()
  !
      Subroutine BCS_ALLOC
    !===================================================================================!
    !
      End Subroutine BCS_ALLOC
  !
      Subroutine BCS_DEALLOC
    !
         If (ALLOCATED(INODEQ)) DEALLOCATE (INODEQ)
         If (ALLOCATED(ICELLQ)) DEALLOCATE (ICELLQ)
#if defined (MULTIPROCESSOR)
         If (PAR) Then
            If (ALLOCATED(RIV_GL2LOC)) DEALLOCATE (RIV_GL2LOC)
         End If
#endif
         If (ALLOCATED(INOPNT)) DEALLOCATE (INOPNT)
         If (ALLOCATED(DWDIS)) DEALLOCATE (DWDIS)
         If (ALLOCATED(WVQDIST)) DEALLOCATE (WVQDIST)
         If (ALLOCATED(PDQDIS)) DEALLOCATE (PDQDIS)
         If (ALLOCATED(PDWQDIS)) DEALLOCATE (PDWQDIS)

         If (ALLOCATED(QBC_TM%TIMES)) DEALLOCATE (QBC_TM%TIMES)
         If (ALLOCATED(PNT_TM%TIMES)) DEALLOCATE (PNT_TM%TIMES)

         If (ALLOCATED(N_ICELLQ)) DEALLOCATE (N_ICELLQ)
    !
      End Subroutine BCS_DEALLOC
  !
End Module MOD_BCS
