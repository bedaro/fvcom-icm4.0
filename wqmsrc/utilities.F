!utilities.F
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
!
!Subroutines :
!	Subroutine PSTOP
!
!
!==============================================================================|
!     UTILITIES FILE                                                           |
!      PSTOP:   HALTS PROGRAM CORRECTLY                                        |
!      WRITE_BANNER: WRITE FVCOM BANNER TO OUTPUT                              |
!==============================================================================|
!
!
!==============================================================================|
!Subroutine PSTOP
!  !==============================================================================|
!#if defined (MULTIPROCESSOR)
!  !
!      Include "mpif.h"
!      Integer IERR
!  !  CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
!      Call MPI_FINALIZE (IERR)
!#endif
!  !
!      Stop
!End Subroutine PSTOP

 SUBROUTINE PSTOP
!==============================================================================|
#if defined (MULTIPROCESSOR)
  Include "mpif.h"
 !USE MPI
  INTEGER IERR, ecode

  ecode = -1
  CALL MPI_ABORT(MPI_COMM_WORLD,ecode,IERR)

  CALL MPI_FINALIZE(IERR)
#endif

  STOP
  END SUBROUTINE PSTOP

SUBROUTINE WRITE_BANNER(IUNIT)

  INTEGER, INTENT(IN) :: IUNIT

   WRITE(IUNIT,*)'!==============================================================================!'
   WRITE(IUNIT,*)'!######  #       #    ####   ####    #       #       #####  #### #      #      !'
   WRITE(IUNIT,*)'!#        #     #    #      #    #   # #   # #         #   #     # #  # #      !'
   WRITE(IUNIT,*)'!####      #   #     #      #    #   #   #   #  #####  #   #     #   #  #      !'
   WRITE(IUNIT,*)'!#          # #      #      #    #   #       #         #   #     #      #      !'
   WRITE(IUNIT,*)'!#           #        ####   ####    #       #       ###### #### #      #      !'
   WRITE(IUNIT,*)'!==============================================================================!'
   WRITE(IUNIT,*)'!                                                                              !'
   WRITE(IUNIT,*)'!=================DOMAIN DECOMPOSITION USING: METIS 4.0.1 =====================!'
   WRITE(IUNIT,*)'!============Copyright 1998, Regents of University of Minnesota================!'
   WRITE(IUNIT,*)'!                                                                              !'

   RETURN
END SUBROUTINE WRITE_BANNER
!==============================================================================|
!
