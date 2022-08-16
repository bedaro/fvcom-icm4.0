!domdec.F
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
!subroutine DOMDEC()
!
!==============================================================================|
!  DECOMPOSE THE DOMAIN BY ELEMENTS USING METIS GRAPH PARTITIONING TOOL        |
!    RETURNS[EL_PID(1:NGL)]                                                    |
!    EL_PID(I) = Processor ID Of Owner of Element I                            |
!==============================================================================|
!
Subroutine DOMDEC (NGL, NPROCS, EL_PID, MSR)
# if defined (MULTIPROCESSOR)
  !
      Use MOD_HYDROVARS, Only: NVG
  !
  !==============================================================================|
      Implicit None
      Include "mpif.h"
      Integer, Intent (In) :: NGL, NPROCS
      Integer, Intent (Out) :: EL_PID (NGL)
      Logical, Intent (In) :: MSR
      Integer, Allocatable :: NVT (:)
      Integer :: I, IERR, ii !NTEMP,!
  !==============================================================================|
  !
  !
  !----------------READ IN NODE LIST FROM ***_grd.dat FILE-----------------------!
  !
  !
      If (MSR) Then
         Allocate (NVT(3*NGL))
     !
     !fix NVT creation
     !  II = 0
     !  DO I=1,NGL
     !    II = II+1
     !    NVT(II)   = NVG(I,1)    !TEMP(I,1)
     !    NVT(II+1) = NVG(I,1)    !TEMP(I,3)
     !    NVT(II+2) = NVG(I,1)    !TEMP(I,2)
     !  ENDDO
         ii = - 2
         Do I = 1, NGL
            ii = ii + 3
            NVT (ii) = NVG (I, 1)!TEMP(I,1)
            NVT (ii+1) = NVG (I, 2)!TEMP(I,3)
            NVT (ii+2) = NVG (I, 3)!TEMP(I,2)
         End Do
     !
     !
     !-------------DECOMPOSE ELEMENTS USING METIS GRAPH PARTITIONING ---------------!
     !
         Call PARTITION (NPROCS, NGL, maxval(NVT), loc(NVT), &
        & loc(EL_PID))
     !
         EL_PID = EL_PID + 1
         Deallocate (NVT)
     !
      End If
  !
      Call MPI_BCAST (EL_PID, NGL, MPI_INTEGER, 0, MPI_COMM_WORLD, &
     & IERR)
  !
  !
# endif
End Subroutine DOMDEC
!==============================================================================|
