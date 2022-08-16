!mod_filenames.F
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
Module MOD_FILENAMES
  !
  !INTEGER ifindext
  !CHARACTER *(*) fnameprefix
  !CHARACTER *(*) fnameext
Contains
  !
  !functions:
  !	function ifindext()
  !	function fnameprefix()
  !	function fnameext()
  !
  !: added a few functions that deals with filenames and extensions etc
  !
  !
      Integer Function ifindext (FILENAME)

    !
         Implicit None
         Integer (4) :: index_dot
         Character (*) :: FILENAME
         index_dot = index (TRIM(FILENAME), '.', .True.)
         ifindext = index_dot
         Return
      End Function ifindext
  !
      Character (1024) Function fnameprefix (FILENAME)
    !
    !Fortran function to find the filename after removing the extension (.dot etc)
    !
         Implicit None
         Integer (4) :: index_dot
         Character (*) :: FILENAME
         Character (Len=Len(FILENAME)) :: FILENAME_PREFIX, FILENAME_NEW
    !
         FILENAME_NEW = TRIM (FILENAME)
         index_dot = index (TRIM(FILENAME), '.', .True.)
         If (index_dot > 0) Then
            FILENAME_PREFIX = TRIM (FILENAME_NEW(1:index_dot-1))
         Else
            FILENAME_PREFIX = TRIM (FILENAME_NEW)
         End If
         fnameprefix = FILENAME_PREFIX
         Return
      End Function fnameprefix
  !
      Character (1024) Function fnameext (FILENAME)
    !
    !fortran function to find extension name if any
    !
         Implicit None
         Integer (4) :: index_dot, length
         Character (*) :: FILENAME
         Character (Len=Len(FILENAME)) :: FILENAME_EXT, FILENAME_NEW
    !
         FILENAME_NEW = TRIM (FILENAME)
         index_dot = index (TRIM(FILENAME), '.', .True.)
         length = len (TRIM(FILENAME_NEW))
    !
         If (index_dot > 0) Then
            FILENAME_EXT = TRIM (FILENAME_NEW(index_dot:length))
         Else
            FILENAME_EXT = ''
         End If
         fnameext = FILENAME_EXT
         Return
      End Function fnameext
  !
  !
End Module MOD_FILENAMES
!
