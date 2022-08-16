!mod_prec.F
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
!**                     Adopted from CE-QUAL-ICM  Model                 **
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
!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
Module MOD_PREC
      Implicit None
#  if defined (MULTIPROCESSOR)
      Include "mpif.h"
#  endif
  !
  !--Single Precision Coding------------------------------------------------------!
#  if !defined (DOUBLE_PRECISION)
      Integer, Parameter :: SP = Selected_Real_Kind (6, 30)
#  if defined (MULTIPROCESSOR)
      Integer, Parameter :: MPI_F = MPI_REAL
#  endif
  !
#  else
  !--Double Precision Coding------------------------------------------------------!
      Integer, Parameter :: SP = Selected_Real_Kind (12, 300)
  !   INTEGER, PARAMETER :: SP = KIND(1.d0)   !Guarantee precision up to that of the machine-compiler-specific double precision 
!
#  if defined (MULTIPROCESSOR)
      Integer, Parameter :: MPI_F = MPI_DOUBLE_PRECISION
#  endif
#  endif
  !
  !
  !
      Integer, Parameter :: DP = Selected_Real_Kind (12, 300)
  !   INTEGER, PARAMETER :: DP     = KIND(1.d0)   !Guarantee precision up to that of the machine-compiler-specific double precis
!
#  if defined (MULTIPROCESSOR)
      Integer, Parameter :: MPI_DP = MPI_DOUBLE_PRECISION
#  endif
  !
  ! USE CORRECT PRECISION WHEN READ NETCDF
      Integer, Parameter :: CDF_PREC = 4
  !
  ! USE CORRECT PRECISION WHEN USING MPI ON CDF DATA
#  if defined (MULTIPROCESSOR)
      Integer, Parameter :: MPI_CDF = MPI_REAL
#  endif
  !
End Module MOD_PREC
!
