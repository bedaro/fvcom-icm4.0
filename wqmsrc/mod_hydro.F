!mod_hydro.F
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
Module MOD_HYDRO
  !
Contains
  !subroutine HYDRO ()
  !
  !: To be combined with mod_hydrovars.F in the future
  !
  !************************************************************************
  !**                  S U B R O U T I N E   H Y D R O                   **
  !************************************************************************
  !
      Subroutine HYDRO ()
         Use MOD_PREC, Only: SP
         Use MOD_HYDROVARS, Only: UNC2, VNC2, WTSNC2, UARD_OBCNNC2, &
        & XFLUX_OBCNC2, DTFANC2, KHNC2, TNC2, SNC2, ELNC2, &
        & num_hyd_ints, NCFILE_PREFIX, NCFILE_SUFFIX, NCFILE_NUMBER, &
        & FORMAT_STR, hydro_dir, hydro_filenumwidth, &
        & hydro_filenumstart, hydro_Nrec, IFNC, NTRECNC, NTHYDRO

         Use MOD_CONTROL, Only: MSR, PAR
    !
         Use MOD_LIMS, Only: MTLOC, KBM1
         Use MOD_SIZES, Only: NCP
!
         Implicit None
#if defined (MULTIPROCESSOR)
         Include "mpif.h"
#endif
         Save
         Integer F, SB, L, JCON, I, J, NC_ID2
         Logical END_OF_FILE
         Real (SP) NXDAY, TDUM
         Real (SP) MASS (0:MTLOC, KBM1, NCP)!
         Character (Len=1024) :: NCFILE
    !
    !************************************************************************
    !**                              Inputs                                **
    !************************************************************************

         num_hyd_ints = hydro_Nrec

         If (NTRECNC .Eq. num_hyd_ints) Then
            NTRECNC = 0
            IFNC = IFNC + 1
         End If
         NTRECNC = NTRECNC + 1 !NTRECNC is time record # within one hydro file
         NTHYDRO = NTHYDRO + 1 !NTHYDRO is overall record # of all hydro files that have been used
    !
         If (NTRECNC .Eq. 1) Then !open new file (IFNC) when trying to read the first record
       !
            NCFILE_NUMBER = ''
            Write (NCFILE_NUMBER(1:hydro_filenumwidth), &
           & TRIM(FORMAT_STR)) (IFNC+hydro_filenumstart-1)
            NCFILE = TRIM (hydro_dir) // TRIM (NCFILE_PREFIX) // TRIM &
           & (NCFILE_NUMBER) // TRIM (NCFILE_SUFFIX)
            If (MSR) WRITE (*,*) 'NCFILE=', TRIM (NCFILE)

            Call NCD_READ_OPEN (NCFILE, UNC2, VNC2, WTSNC2, &
           & UARD_OBCNNC2, XFLUX_OBCNC2, DTFANC2, KHNC2, ELNC2, TNC2, &
           & SNC2, NTRECNC)

         Else
       !
            Call NCD_READ (NCFILE, UNC2, VNC2, WTSNC2, UARD_OBCNNC2, &
           & XFLUX_OBCNC2, DTFANC2, KHNC2, ELNC2, TNC2, SNC2, NTRECNC)
       !
         End If
    !
         Return
      End Subroutine HYDRO
  !
  !
#if defined (MULTIPROCESSOR)
      Subroutine BROADCAST_HYDRO (MPI_SRC, UL, VL, WTSL, UARD_OBCNL, &
     & XFLUX_OBCL, DTFAL, KHL, ELL, TL, SL)
    !
    !
    !Subroutine to broadcast hydrodynamic fields to child processes from master process
    !
         Use MOD_LIMS, Only: NLOC, MLOC, KB, KBM1, NTLOC, MTLOC
    !
         Use MOD_WQM, Only: UL_GL, VL_GL, WTSL_GL, UARD_OBCN_GL, &
        & XFLUX_OBC_GL, DTFAL_GL, KHL_GL, ELL_GL, SL_GL, TL_GL

         Use MOD_SIZES, Only: MGL, NGL, NOBTY

         Use MOD_CONTROL, Only: PAR

         Use MOD_PAR, Only: NLID, EGID, NGID, NHE, HE_LST, NHN, HN_LST

         Use MOD_BCMAP, Only: IOBCN, I_OBC_GL

         Use MOD_PREC, Only: MPI_F, SP, MPI_CDF
         Implicit None
         Include "mpif.h"
    !
    !passed variables
         Integer, Intent (In) :: MPI_SRC ! ID of the source process for broadcasting
         Real (SP), Dimension (0:NTLOC, KB), Intent (Out) :: UL, VL
         Real (SP), Dimension (0:MTLOC, KB), Intent (Out) :: KHL, WTSL
         Real (SP), Dimension (0:MTLOC, KBM1), Intent (Out) :: TL, SL 
         Real (SP), Dimension (0:MTLOC), Intent (Out) :: ELL, DTFAL
         Real (SP), Dimension (0:NOBTY+1), Intent (Out) :: UARD_OBCNL
         Real (SP), Dimension (0:NOBTY, KBM1), Intent (Out) :: &
        & XFLUX_OBCL
    !
    !local variables
    !
         Integer :: MPI_COUNT_TMP = 0 !count to send
         Integer :: I, J, K
         Integer :: IERR
    !Broadcast global variables
         If (PAR) Then

            MPI_COUNT_TMP = NGL * KBM1
            Call MPI_BCAST (UL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = NGL * KBM1
            Call MPI_BCAST (VL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL * KB
            Call MPI_BCAST (WTSL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = NOBTY + 1
            Call MPI_BCAST (UARD_OBCN_GL, MPI_COUNT_TMP, MPI_CDF, &
           & MPI_SRC, MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = NOBTY * KBM1
            Call MPI_BCAST (XFLUX_OBC_GL, MPI_COUNT_TMP, MPI_CDF, &
           & MPI_SRC, MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL
            Call MPI_BCAST (DTFAL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL * KB
            Call MPI_BCAST (KHL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL
            Call MPI_BCAST (ELL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL * KBM1
            Call MPI_BCAST (SL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL * KBM1
            Call MPI_BCAST (TL_GL, MPI_COUNT_TMP, MPI_CDF, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
         End If
    !
    !---then push the global arrays into local arays including halo nodes and halo elements
         If (PAR) Then
            Do I = 1, MLOC
               ELL (I) = ELL_GL (NGID(I))
            End Do
       !
            Do I = 1, NHN
               ELL (I+MLOC) = ELL_GL (HN_LST(I))
            End Do
       !
            Do K = 1, KBM1
               Do I = 1, MLOC
                  SL (I, K) = SL_GL (NGID(I), K)
                  IF (SL(I,K) < 0.0000001) THEN 
                    SL(I,K) = 0.0
                  END IF
               End Do
          !
               Do I = 1, NHN
                  SL (I+MLOC, K) = SL_GL (HN_LST(I), K)
                  IF (SL(I,K) < 0.0000001) THEN 
                    SL(I,K) = 0.0
                  END IF
               End Do
            End Do
       !
            Do K = 1, KBM1
               Do I = 1, MLOC
                  TL (I, K) = TL_GL (NGID(I), K)
               End Do
          !
               Do I = 1, NHN
                  TL (I+MLOC, K) = TL_GL (HN_LST(I), K)
               End Do
            End Do
       !
            Do K = 1, KBM1
               Do I = 1, NLOC
                  UL (I, K) = UL_GL (EGID(I), K)
               End Do
          !
               Do I = 1, NHE
                  UL (I+NLOC, K) = UL_GL (HE_LST(I), K)
               End Do
            End Do
       !
            Do K = 1, KBM1
               Do I = 1, NLOC
                  VL (I, K) = VL_GL (EGID(I), K)
               End Do
          !
               Do I = 1, NHE
                  VL (I+NLOC, K) = VL_GL (HE_LST(I), K)
               End Do
            End Do
       !
            Do K = 1, KB
               Do I = 1, MLOC
                  WTSL (I, K) = WTSL_GL (NGID(I), K)
               End Do
          !
               Do I = 1, NHN
                  WTSL (I+MLOC, K) = WTSL_GL (HN_LST(I), K)
               End Do
            End Do
       !
            UARD_OBCNL = - 99999
            K = 0
            Do I = 1, NOBTY
               J = NLID (I_OBC_GL(I))
               If (J .Ne. 0) Then
                  K = K + 1
                  UARD_OBCNL (K) = UARD_OBCN_GL (I)
               End If
            End Do
            If (K .Ne. IOBCN) Then
               Write (*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
               Call pstop
            End If
       !
            XFLUX_OBCL = - 99999
            K = 0
            Do I = 1, NOBTY
               J = NLID (I_OBC_GL(I))
               If (J .Ne. 0) Then
                  K = K + 1
                  XFLUX_OBCL (K, :) = XFLUX_OBC_GL (I, :)
               End If
            End Do
       !
            If (K .Ne. IOBCN) Then
               Write (*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
               Call pstop
            End If
       !
            Do I = 1, MLOC
               DTFAL (I) = DTFAL_GL (NGID(I))
            End Do
       !
            Do I = 1, NHN
               DTFAL (I+MLOC) = DTFAL_GL (HN_LST(I))
            End Do
       !
            Do K = 1, KB
               Do I = 1, MLOC
                  KHL (I, K) = KHL_GL (NGID(I), K)
               End Do
          !
               Do I = 1, NHN
                  KHL (I+MLOC, K) = KHL_GL (HN_LST(I), K)
               End Do
            End Do
       !
         End If
    !
         Return
      End Subroutine BROADCAST_HYDRO
#endif
  !
#if defined (MULTIPROCESSOR)
      Subroutine BROADCAST_HYDRO_REDUCED (MPI_SRC, UARD_OBCNL, &
     & XFLUX_OBCL, DTFAL, ELL)
    !
    !Subroutine to broadcast hydrodynamic fields to child processes from master process
    !
    !
         Use MOD_LIMS, Only: NTLOC, MLOC, MTLOC, KBM1
    !
         Use MOD_PREC, Only: MPI_F, SP
    !
         Use MOD_SIZES, Only: MGL, NOBTY
    !
         Use MOD_WQM, Only: UARD_OBCN_GL, XFLUX_OBC_GL, DTFAL_GL, &
        & ELL_GL
    !Wen Long took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
         Use MOD_CONTROL, Only: PAR
    !
         Use MOD_BCMAP, Only: IOBCN, I_OBC_GL
    !
         Use MOD_PAR, Only: NLID, NGID, NHN, HN_LST
         Implicit None
         Include "mpif.h"
    !passed variables
         Integer, Intent (In) :: MPI_SRC ! ID of the source process for broadcasting
         Real (SP), Dimension (0:MTLOC), Intent (Out) :: ELL, DTFAL
         Real (SP), Dimension (0:NOBTY+1), Intent (Out) :: UARD_OBCNL
         Real (SP), Dimension (0:NOBTY, KBM1), Intent (Out) :: &
        & XFLUX_OBCL
    !local variables
         Integer :: MPI_COUNT_TMP = 0 !count to send
         Integer :: I, J, K
         Integer :: IERR
    !Broadcast global variables
         If (PAR) Then
            MPI_COUNT_TMP = NOBTY + 1
            Call MPI_BCAST (UARD_OBCN_GL(:), MPI_COUNT_TMP, MPI_F, &
           & MPI_SRC, MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = NOBTY * KBM1
            Call MPI_BCAST (XFLUX_OBC_GL(:, :), MPI_COUNT_TMP, MPI_F, &
           & MPI_SRC, MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL
            Call MPI_BCAST (DTFAL_GL(:), MPI_COUNT_TMP, MPI_F, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
            MPI_COUNT_TMP = MGL
            Call MPI_BCAST (ELL_GL(:), MPI_COUNT_TMP, MPI_F, MPI_SRC, &
           & MPI_COMM_WORLD, IERR)
         End If
    !---then push the global arrays into local arays including halo nodes and halo elements
         If (PAR) Then
            Do I = 1, MLOC
          !WRITE(*,*)'ELL_GL Debug here here 03 '
               ELL (I) = ELL_GL (NGID(I))
          !WRITE(*,*)'ELL_GL Debug here here 3 '
            End Do
       !
            Do I = 1, NHN
          !WRITE(*,*)'ELL_GL Debug here here 04 '
               ELL (I+MLOC) = ELL_GL (HN_LST(I))
          !WRITE(*,*)'ELL_GL Debug here here 4 '
            End Do
       !
            UARD_OBCNL = - 99999
            K = 0
            Do I = 1, NOBTY
               J = NLID (I_OBC_GL(I))
               If (J .Ne. 0) Then
                  K = K + 1
                  UARD_OBCNL (K) = UARD_OBCN_GL (I)
               End If
            End Do
            If (K .Ne. IOBCN) Then
               Write (*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
               Call pstop
            End If
       !
            XFLUX_OBCL = - 99999
            K = 0
            Do I = 1, NOBTY
               J = NLID (I_OBC_GL(I))
               If (J .Ne. 0) Then
                  K = K + 1
                  XFLUX_OBCL (K, :) = XFLUX_OBC_GL (I, :)
               End If
            End Do
       !
            If (K .Ne. IOBCN) Then
               Write (*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
               Call pstop
            End If
       !
            Do I = 1, MLOC
               DTFAL (I) = DTFAL_GL (NGID(I))
            End Do
       !
            Do I = 1, NHN
               DTFAL (I+MLOC) = DTFAL_GL (HN_LST(I))
            End Do
       !
         End If
    !
         Return
    !
      End Subroutine BROADCAST_HYDRO_REDUCED
#endif
  !
#if defined (MULTIPROCESSOR)
      Subroutine SCATTER (N1, NL, NT, NTG, KT, MYID, NPROCS, GM, CM, &
     & NID, NID_X, AG, AL, BG, BL, CG, CL)
    !======================================================================!
    ! PASS global array AG from master proc to child procs local array AL  !
    ! Also can do the same for BG -->BL, CG->CL if present                 !
    !-----------------------------------------------------------------------
    !
         Use MOD_TYPES, Only: GMAP, COMM
         Use MOD_PAR, Only: EXCHANGE
         Use MOD_CONTROL, Only: PAR !,SERIAL,MSR
         Use MOD_PREC, Only: MPI_F, CDF_PREC, SP, DP, MPI_DP
    !
         Use MOD_BUFFERS, Only: SNDBUF, RCVBUF
    !
         Implicit None
         Include "mpif.h"
    !
         Integer, Intent (In) :: N1, NL, NT, NTG, KT, MYID, NPROCS !total number of procs
    !
         Type (GMAP), Intent (In) :: GM (NPROCS)!mapping (element or node)
         Real (CDF_PREC), Intent (In) :: AG (1:NTG, KT)!global arry AG always starts at 1
         Real (SP), Intent (Out) :: AL (N1:NT, KT)!local array AL, size includs halos, NT-NL = number of halo elemetns/nodes
         Real (CDF_PREC), Intent (In), Optional :: BG (1:NTG, KT)
         Real (SP), Intent (Out), Optional :: BL (N1:NT, KT)
         Real (CDF_PREC), Intent (In), Optional :: CG (1:NTG, KT)
         Real (SP), Intent (Out), Optional :: CL (N1:NT, KT)
         Type (COMM), Intent (In) :: CM (NPROCS)!communicator (for elements or nodes)
         Integer, Intent (In) :: NID (0:NTG)!local index of a global element/node,
    !      0 if out of current proc
         Integer, Intent (In) :: NID_X (0:NTG)!local index of a global element/node including halos,
    !      0 if out of current proca
    !local variables
         Integer :: IERR
         Integer :: MPI_COUNT_TMP = 0 !count to send
         Integer :: I, J, K
         Integer :: IP
         Integer :: TSND
         Integer :: TRCV
         Integer :: ISNDTAG
         Integer :: IRCVTAG
         Integer :: ISND
         Integer :: IRCV

         Integer :: NMSG
         Integer :: IREQR (NPROCS), IREQS (NPROCS)
         Integer :: STAT (MPI_STATUS_SIZE), ISTATR (MPI_STATUS_SIZE, &
        & NPROCS)
         Integer :: INDX
         Integer :: NVAR_PASS = 1 !number of variables to pass
    !
         Logical :: BYES = .False., CYES = .False. !Yes if CG and CL are passed in
    !
#if(0)  /* skip debugging after testing */
    !debugging local variables
         Real (SP) :: SLEEPTIME
         Character (Len=1024) :: COMMAND_STR, SLEEPTIME_STR
#endif
    !
         If (PAR) Then
       !
       !master proc knows every child procs' global ID (element or nodes) through mapping GM
       !GM contains the following attributes:
       !        NSIZE             --- size of the footprint on a proc
       !        LOC_2_GL          --- global index of a local arraya
       !Note that GM is only avaiable on master proc
       !
       !on the other hand, all procs know the correspoonding global ID through NID and NID_X
       !NID(I)   gives local value of global element or nodes  (0:NTG) (RETURNS 0 IF I IS NOT AN INTERIOR element or node)
       !NID_X(I) gives local value of global element or nodes including halos  (0:NTG) (RETURNS 0 IF I IS NOT AN INTERIOR OR HALO el
       !
       !
       !sanity check on optional argments
            If (PRESENT(BG) .neqv. PRESENT(BL)) Then
               Write (*,*) 'ERROR calling SCATTER(), BG and BL must be &
              &passed at the same time'
               Call pstop
            End If
            If (PRESENT(CG) .neqv. PRESENT(CL)) Then
               Write (*,*) 'ERROR calling SCATTER(), CG and CL must be &
              &passed at the same time'
               Call pstop
            End If
       !
            If (PRESENT(BG) .And. PRESENT(BL)) BYES = .True.
            If (PRESENT(CG) .And. PRESENT(CL)) CYES = .True.
       !
            If (CYES .Eqv. .True. .And. BYES .Eqv. .False.) Then
               Write (*,*) 'ERROR calling SCATTER(), BG and BL must be &
              &avaiable when CG and CL are used'
               Call pstop
            End If
       !
            If (BYES) NVAR_PASS = NVAR_PASS + 1
            If (CYES) NVAR_PASS = NVAR_PASS + 1
       !
            If ( .Not. ALLOCATED(RCVBUF)) ALLOCATE &
           & (RCVBUF(1:NTG*KT*NVAR_PASS))!actually can find maximum length across all procs so do not
       !really need size 1:NTG, however child does not have GM
       !so this buffer size is hard for child to allocate
       !hence we just use NTG as an overkill but simple
            If ( .Not. ALLOCATED(SNDBUF)) ALLOCATE &
           & (SNDBUF(1:NTG*KT*NVAR_PASS))
       !
       !Child receives  data
            If (MYID /= 1) Then
               TRCV = 0
               MPI_COUNT_TMP = NL * NVAR_PASS * KT ! number of data points to receive (locally)
               RCVBUF (:) = 0.0
               TRCV = TRCV + 1
               IRCVTAG = 1000 + MYID
          !non-blocking receiving method with receiving status in IREQR to querry later
               Call MPI_IRECV (RCVBUF(1), MPI_COUNT_TMP, MPI_F, 0, &
              & IRCVTAG, MPI_COMM_WORLD, IREQR(TRCV), IERR)
            End If
       !
       !
       !Master sends out data to child processes
            If (MYID == 1) Then
               TSND = 0 !count number of sends by master proc
               Do IP = 2, NPROCS !send data to child procs (2:NPROCS)
             !Master finds chunk of data for each child
             !find size of buffer to send
                  SNDBUF (:) = 0.0
             !fill out buffer
                  ISND = 0
                  Do I = 1, GM(IP)%NSIZE !loop through all values for proc IP
                     Do K = 1, KT
                        ISND = ISND + 1
                        SNDBUF (ISND) = AG (GM(IP)%LOC_2_GL(I), K)!master fetches global data for this child
                   !and puts it in the sending buffer
                     End Do
                  End Do
             !
                  If (BYES) Then !pass array BG
                     Do I = 1, GM(IP)%NSIZE !loop through all values for proc IP
                        Do K = 1, KT
                           ISND = ISND + 1
                           SNDBUF (ISND) = BG (GM(IP)%LOC_2_GL(I), K)!master fetches global data for this child
                      !and puts it in the sending buffer
                        End Do
                     End Do
                  End If
             !
                  If (CYES) Then !pass array CG
                     Do I = 1, GM(IP)%NSIZE !loop through all values for proc IP
                        Do K = 1, 1
                           ISND = ISND + 1
                           SNDBUF (ISND) = CG (GM(IP)%LOC_2_GL(I), K)!master fetches global data for this child
                      !and puts it in the sending buffer
                        End Do
                     End Do
                  End If
             !
                  MPI_COUNT_TMP = ISND !number of data points to send to proc IP.
             !Master sends chunk of data to each child
                  ISNDTAG = IP + 1000 !Note possible bug here with 1000, we need to make sure ISNDTAG is unique for each process IP
                  TSND = TSND + 1
             !
             !Send with non-blocking method, with request in IREQS to querry later
                  Call MPI_ISEND (SNDBUF(1), MPI_COUNT_TMP, MPI_F, &
                 & IP-1, ISNDTAG, MPI_COMM_WORLD, IREQS(TSND), IERR)
             !
             !or use the blocking send
             !CALL MPI_SEND(SNDBUF(1),MPI_COUNT_TMP, MPI_F, IP-1, ISNDTAG,MPI_COMM_WORLD,IERR)
               End Do
            End If
       !
       !           !Child receives data
       !           IF(MYID/=1)THEN
       !                 TRCV=0
       !                 MPI_COUNT_TMP=NL*NVAR_PASS*KT    ! number of data points to receive (locally)
       !                 RCVBUF(:)=0.0
       !                 TRCV=TRCV+1
       !                 IRCVTAG= 1000+MYID
       !                !non-blocking receiving method with receiving status in IREQR to querry later
       !                 CALL MPI_IRECV(RCVBUF(1),MPI_COUNT_TMP,MPI_F,0,IRCVTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
       !
       !                !or blocking receve
       !                ! CALL MPI_RECV(RCVBUF(1),MPI_COUNT_TMP,MPI_F,0,IRCVTAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,IERR)
       !           ENDIF
       !
       !Loop over procs until a message is received and unpack the data
            If (MYID /= 1) Then
               Do NMSG = 1, TRCV
                  Call MPI_WAITANY (TRCV, IREQR, INDX, STAT, IERR, &
                 & RCVBUF, SNDBUF)!use this when use the non-blocking recieve method
             !Unpack the message
                  IRCV = 0 !index in the buffer
                  Do I = 1, NL !NL must be equal to GM(MYID)%NSIZE although GM is only available in master
                     Do K = 1, KT
                        IRCV = IRCV + 1
                        AL (I, K) = RCVBUF (IRCV)!Child retrieves the data from the receiving buffer
                     End Do
                  End Do
             !
                  If (BYES) Then !retrive array BL
                     Do I = 1, NL !NL must be equal to GM(MYID)%NSIZE although GM is onlya vailable in master
                        Do K = 1, KT
                           IRCV = IRCV + 1
                           BL (I, K) = RCVBUF (IRCV)!Child retrieves the data from the receiving buffer
                        End Do
                     End Do
                  End If
             !
                  If (CYES) Then !retrieve array CL
                     Do I = 1, NL !NL must be equal to GM(MYID)%NSIZE although GM is onlya vailable in master
                        Do K = 1, KT
                           IRCV = IRCV + 1
                           CL (I, K) = RCVBUF (IRCV)!Child retrieves the data from the receiving buffer
                        End Do
                     End Do
                  End If
             !
               End Do
            End If
       !
       !WAIT FOR COMPLETION OF NON-BLOCKING SENDS (turn on this when using non-blocking sends)
            If (MYID == 1) Then
               Call MPI_WAITALL (TSND, IREQS, ISTATR, IERR, RCVBUF, &
              & SNDBUF)! Add RBUF,SBUF to fool compiler into NOT optimizing into a bug
            End If
       !
            If (ALLOCATED(RCVBUF)) DEALLOCATE (RCVBUF)
            If (ALLOCATED(SNDBUF)) DEALLOCATE (SNDBUF)
       !
       !Master gets its own copy of the global data into local array
            If (MYID == 1) Then
               Do K = 1, KT
                  Do I = 1, GM(MYID)%NSIZE
                     AL (I, K) = AG (GM(MYID)%LOC_2_GL(I), K)
                     If (BYES) BL (I, K) = BG (GM(MYID)%LOC_2_GL(I), K)
                     If (CYES) CL (I, K) = CG (GM(MYID)%LOC_2_GL(I), K)
                  End Do
               End Do
            End If
       !
            Call MPI_BARRIER (MPI_COMM_WORLD, IERR)
       !
       !Note here CM is type COMM for elements and nodes, i.e. element communicator and node communicator
       !They have the followinng properties:
       !
       !SND: TRUE if you are to send to processor
       !RCV: TRUE IF you are to receive from processor
       !NSND: Number of data to send to processor
       !NRCV: Number of data to receive from processor
       !SNDP: Array pointing to locations to send to processor
       !RCVP: Array pointing to locations received from process
       !RCPT: Pointer to location in receive buffer
       !
       !The EXCHANGE() subroutine will prepare these properties and then send or receive
       !and then set the halo elements/nodes with values in the received buffer
       !
       !
       !Call the EXCHANGE() subroutine after scattering to make sure
       !each process gets its halo node/element information from neighbors
       !
       !
            Call EXCHANGE (CM, NT, KT, MYID, NPROCS, AL)
       !
            If (BYES) Call EXCHANGE (CM, NT, KT, MYID, NPROCS, BL)
            If (CYES) Call EXCHANGE (CM, NT, KT, MYID, NPROCS, CL)
       !
       !
#if(0) /*skip debugging and testing */
       !
       !Test : master gathers the same data that was scattered to child processes
       !       and check if the gathered data is the same as what was sent out
       !
       !Master prints the global data
            Do IP = 1, NPROCS
               If (IP .Eq. 1 .And. MYID .Eq. 1) Then
                  Do K = 1, KT
                     Do I = 1, NTG
                        Write (*, '(A9,1x,I10,1x,A2,1x,I10,1x,A8,F20.5)&
                       &') 'M-A,GIND=', I, 'K=', K, 'AG(I,K)=', AG (I, &
                       & K)
                        If (BYES) WRITE (*, '(A9,1x,I10,1x,A2,1x,I10,1x&
                       &,A8,F20.5)') 'M-B,GIND=', I, 'K=', K, 'BG(I,K)=&
                       &', BG (I, K)
                        If (CYES) WRITE (*, '(A9,1x,I10,1x,A2,1x,I10,1x&
                       &,A8,F20.5)') 'M-C,GIND=', I, 'K=', K, 'CG(I,K)=&
                       &', CG (I, K)
                     End Do
                  End Do
               End If
               SLEEPTIME = Mod (MYID, 2)!(MYID+1)*1  ! sleep time in microseconds
               SLEEPTIME_STR = ''
               Write (SLEEPTIME_STR, '(I0)') Int (SLEEPTIME)! put number of microseconds in integer in a string
               COMMAND_STR = 'usleep' // ' ' // SLEEPTIME_STR ! command line string for system command
               Call SYSTEM (TRIM(COMMAND_STR))
               Call MPI_BARRIER (MPI_COMM_WORLD, IERR)
            End Do
       !
       !every one prints local data
            Do IP = 1, NPROCS !printing is serialized by nesting MPI_BARRIER in this IP loop
               If (IP .Eq. MYID) Then
                  Do K = 1, KT
                     Do I = 1, NTG !loop globally
                        If (NID(I) .Ne. 0) Then
                           Write (*, '(A9,1x,I10,1x,A2,1x,I10,1x,A8,1x,&
                          &F20.5)') 'C-A,GIND=', I, 'K=', K, 'AL(I,K)='&
                          & , AL (NID(I), K), MYID
                           If (BYES) WRITE (*, '(A9,1x,I10,1x,A2,1x,I10&
                          &,1x,A8,1x,F20.5)') 'C-B,GIND=', I, 'K=', K, &
                          & 'BL(I,K)=', BL (NID(I), K), MYID
                           If (CYES) WRITE (*, '(A9,1x,I10,1x,A2,1x,I10&
                          &,1x,A8,1x,F20.5)') 'C-C,GIND=', I, 'K=', K, &
                          & 'CL(I,K)=', CL (NID(I), K), MYID
                        End If
                        If (NID_X(I) .Ne. 0 .And. NID(I) .Eq. 0) Then !print the HALO data only, this is to make sure HALO data are correct
                           Write (*, '(A10,1x,I10,1x,A2,1x,I10,1x,A8,1x&
                          &,F20.5)') 'CH-A,GIND=', I, 'K=', K, 'AL(I,K)&
                          &=', AL (NID_X(I), K), MYID
                           If (BYES) WRITE (*, '(A10,1x,I10,1x,A2,1x,I1&
                          &0,1x,A8,1x,F20.5)') 'CH-B,GIND=', I, 'K=', &
                          & K, 'BL(I,K)=', BL (NID_X(I), K), MYID
                           If (CYES) WRITE (*, '(A10,1x,I10,1x,A2,1x,I1&
                          &0,1x,A8,1x,F20.5)') 'CH-C,GIND=', I, 'K=', &
                          & K, 'CL(I,K)=', CL (NID_X(I), K), MYID
                        End If
                     End Do
                  End Do
               End If
          !
               SLEEPTIME = Mod (MYID, 2)!Sleep shortly so that different procs spend different amount time in IP loop
               SLEEPTIME_STR = '' !to make sure this printing is serialized
               Write (SLEEPTIME_STR, '(I0)') Int (SLEEPTIME)! put number of microseconds in integer in a string
               COMMAND_STR = 'usleep' // ' ' // SLEEPTIME_STR ! command line string for system command
               Call SYSTEM (TRIM(COMMAND_STR))
               Call MPI_BARRIER (MPI_COMM_WORLD, IERR)
            End Do
            Call pstop
#endif
         End If
         Return
      End Subroutine SCATTER
#endif
  !
End Module MOD_HYDRO
