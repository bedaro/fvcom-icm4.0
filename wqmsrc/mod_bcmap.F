!mod_bcmap.F
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
Module MOD_BCMAP
  !
      Use MOD_PREC, Only: SP
      Use MOD_TYPES, Only: BC
  !
  !Moved the following here from mod_obcs.F
      Integer :: IOBCN_GL !!GLOBAL NUMBER OF OPEN BOUNDARY NODES
      Integer :: IOBCN !!LOCAL NUMBER OF OPEN BOUNDARY NODES
      Integer, Allocatable :: I_OBC_GL (:)!!GLOBAL ID OF OPEN BOUNDARY NODES
      Integer, Allocatable :: I_OBC_N (:)!!OPEN BOUNDARY NODE LIST
  !
  !Moved here from mod_bcs.F
      Real (SP), Allocatable :: WQOBC (:, :, :, :)!!WATER QUALITY DATA AT BOUNDARY
      Type (BC) :: NUT_TM !!TIME MAPPING FOR NUTRIENT OBC  !should be
  !
Contains
  !
  !Subroutine BCMAP()
  !Subroutine BCMAP_DEALLOC
  !
  !--------------------------------------------
      Subroutine BCMAP
    !--------------------------------------------
         Use MOD_SIZES, Only: NCP, MGL
    !
         Use MOD_LIMS, Only: KBM1, NPROCS
!
         Use MOD_CONTROL, Only: SERIAL, MSR, PAR, CASENAME
    !
         Use MOD_PREC, Only: SP
    	 USE MOD_FILEINFO, ONLY: INOBC
         Use MOD_UTILS, Only: FOPEN
# if defined (MULTIPROCESSOR)
         Use MOD_PAR, Only: NLID
# endif
         Implicit None
    !
    !
#if defined (MULTIPROCESSOR)
         Include "mpif.h"
#endif
         Integer :: I, I1, NCNT, IERR, J, JJ, JT, K 
         Integer, Allocatable :: TEMP1 (:)

         Real (SP), Allocatable :: TEMP9 (:, :, :, :)
         Real (SP) :: FTEMP1
         Real (SP), Allocatable :: WQOBC_GL (:, :, :, :)
    !
         Integer ::  NTMP
         Character (Len=100) :: ISTR
    !

    !------------------------------------------------------------------------------!
    !
    !==============================================================================|
    !   OPEN BOUNDARY CONDITION NODES                                              |
    !==============================================================================|
    !
    !----------------------------REPORT--------------------------------------------!
    !
         If (MSR) WRITE (*,*) '!'
         If (MSR) WRITE (*,*) '!           SETTING UP OPEN BOUNDARY NOD&
        &ES  '
         If (MSR) WRITE (*,*) '!'
    !
         IOBCN = 0
    !   IBCN  = 0
    !
         ISTR = "./inputs/" // trim (CASENAME)
         Call FOPEN (INOBC, trim(ISTR)//'_obc_wq.dat', "cfr")
         Rewind (INOBC)
    !
         Read (INOBC,*) IOBCN_GL
    !
         If (IOBCN_GL > 0) Then
       !
       !------------Read in Open Boundary Nodes and Temperature/Salinity Conditions---!
       !
       !
       !----Input Number of Data Times for T/S at Open Boundary Every Sigma Layers---!
       !
            Allocate (I_OBC_GL(IOBCN_GL))
       !
            Do I = 1, IOBCN_GL
               Read (INOBC,*) I1, I_OBC_GL (I)
            End Do
            NCNT = 0
            Do WHILE (.True.)
               Read (INOBC,*, END=10) FTEMP1
               Do JJ = 1, NCP
                  Do J = 1, IOBCN_GL
                     Read (INOBC,*)
                  End Do
               End Do
               NCNT = NCNT + 1
            End Do
10          Continue
            Rewind (INOBC)
       !
            If (NCNT == 0) WRITE (*,*) 'NO WQ DATA at BND'
       !
            NUT_TM%NTIMES = NCNT
       !
       !    READ nutrient open boundary values
       !
            Allocate (NUT_TM%TIMES(NUT_TM%NTIMES))
            NUT_TM%TIMES = 0.
            Allocate (WQOBC_GL(IOBCN_GL, KBM1, NUT_TM%NTIMES, NCP))
            WQOBC_GL = 0.
       !
            Read (INOBC,*)
            Do I = 1, IOBCN_GL
               Read (INOBC,*)
            End Do
       !
            Do JT = 1, NUT_TM%NTIMES
               Read (INOBC,*) NUT_TM%TIMES(JT)!time(hours)
               Do JJ = 1, NCP
                  Do I = 1, IOBCN_GL
                     Read (INOBC,*) NTMP, (WQOBC_GL(I, K, JT, JJ), K=1, &
                    & KBM1)!temp in OB sigma layer
                  End Do
               End Do
          !          NUT_TM%TIMES = 24.0*(NUT_TM%TIMES) !shift to model hours
            End Do
            Close (INOBC)
       !
       !---------Note in log file if TALK is a function of Salinity at open boundary---------!
       !            (done in bcond_wqm.F, not here) -                      !
#ifdef TA_OBC_INI_FROM_S
            Write (*,*) '!  TALK at OBC and INI      :    function of s&
           &alinity'
#endif
#ifdef TA_INI_FROM_S
            Write (*,*) '!  TALK at INI      :    function of salinity'
#endif
       !--------- Jan 7 2016---------------------------------------------------------------!
       !
       !----------------------Make Sure It Is In Global Domain------------------------!
       !
            Do I = 1, IOBCN_GL
               If ((I_OBC_GL(I) > MGL)) Then
                  Write (*,*) '==================ERROR=================&
                 &================='
                  Write (*,*) 'OPEN BOUNDARY NODE NUMBER', I, 'IS NOT I&
                 &N THE'
                  Write (*,*) 'GLOBAL DOMAIN'
                  Write (*,*) 'CHECK INPUT FILE AND ENSURE OPEN BOUNDAR&
                 &Y NODES <= ', MGL
                  Write (*,*) '========================================&
                 &================='
                  Call PSTOP
               End If
            End Do
       !
       !----------Shift Open Boundary Node List,Type,Salt,and Temp to Local-----------!
       !
            If (SERIAL) Then
          !
               IOBCN = IOBCN_GL
               Allocate (I_OBC_N(IOBCN))
               I_OBC_N = I_OBC_GL
          !
               Allocate (WQOBC(IOBCN, KBM1, NUT_TM%NTIMES, NCP))
               WQOBC = 0.0
               WQOBC (:, :, :, :) = WQOBC_GL (:, :, :, :)
          !
            End If
       !
#    if defined (MULTIPROCESSOR)
            If (PAR) Then
               Allocate (TEMP1(IOBCN_GL))
               Allocate (TEMP9(IOBCN_GL, KBM1, NUT_TM%NTIMES, NCP))
               TEMP9 = 0 
               NCNT = 0
          !!SET UP LOCAL OPEN BOUNDARY NODES
               Do I = 1, IOBCN_GL
                  I1 = NLID (I_OBC_GL(I))
                  If (I1 /= 0) Then
                     NCNT = NCNT + 1
                     TEMP1 (NCNT) = I1
                     TEMP9 (NCNT, :, :, :) = WQOBC_GL (I, :, :, :)
                  End If
               End Do
               IOBCN = NCNT
          !
               If (NCNT > 0) Then
                  Allocate (I_OBC_N(NCNT))
                  I_OBC_N = TEMP1 (1:NCNT)
                  Allocate (WQOBC(NCNT, KBM1, NUT_TM%NTIMES, NCP))
                  WQOBC = 0.0 
                  WQOBC (1:NCNT, :, :, :) = TEMP9 (1:NCNT, :, :, :)
               End If
          !
               Deallocate (TEMP9)
               Deallocate (TEMP1)
            End If
#    endif
       !
            Deallocate (WQOBC_GL)
            Write (*,*) '!  TEMP/SALT OBC READ      :    COMPLETE'
         End If !!IOBCN_GL > 0
    !
    !==============================================================================|
    !   REPORT AND CHECK RESULTS                                                   |
    !==============================================================================|
         Allocate (TEMP1(NPROCS))
         TEMP1 (1) = IOBCN
    !
# if defined (MULTIPROCESSOR)
         Call MPI_GATHER (IOBCN, 1, MPI_INTEGER, TEMP1, 1, MPI_INTEGER, &
        & 0, MPI_COMM_WORLD, IERR)
# endif
    !
    !
         If (MSR) WRITE (*, 100) '!  IOBCN                 :', &
        & IOBCN_GL, (TEMP1(I), I=1, NPROCS)
         Deallocate (TEMP1)
    !
         Return
100      Format (1 X, A26, I6, " =>", 2 X, 4(I5, 1 H,))
      End Subroutine BCMAP
  !==============================================================================|
  !
      Subroutine BCMAP_DEALLOC
    !
    !Moved here from wqm_main.F
         If (ALLOCATED(I_OBC_N)) DEALLOCATE (I_OBC_N)
         If (ALLOCATED(I_OBC_GL)) DEALLOCATE (I_OBC_GL)
    !
    !Added deallocation WQOBC
    !
         If (ALLOCATED(WQOBC)) DEALLOCATE (WQOBC)
    !
    !These deallocations were missing
         If (ALLOCATED(NUT_TM%TIMES)) DEALLOCATE (NUT_TM%TIMES)
    !
    !
    !
      End Subroutine BCMAP_DEALLOC
  !
End Module MOD_BCMAP
