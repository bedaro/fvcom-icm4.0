!genmap.F
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
!subroutine GENMAP()
!
!==============================================================================!
!  SET UP GLOBAL <--> LOCAL NODE AND ELEMENT NUMBERING MAPS                    !
!  SET UP INTERPROCESSOR COMMUNICATION                                         !
!==============================================================================!
!
Subroutine GENMAP
# if defined (MULTIPROCESSOR)
  !==============================================================================!
      Use MOD_PREC, Only: SP
      Use MOD_LIMS, Only: NLOC, MLOC, MYID, NPROCS, KB, NTLOC, MTLOC
      Use MOD_HYDROVARS, Only: NVG
  !
  !Commented out the following CONTROL module for it is already included in MOD_HYDROVARS
      Use MOD_CONTROL, Only: SERIAL, MSR
      Use MOD_SIZES, Only: MGL, NGL
  !
      Use MOD_PAR, Only: EL_PID, ELID, NLID, ELID_X, NLID_X, EGID, &
     & NGID, NHE, HE_LST, HE_OWN, NBN, MX_MLT, BN_LST, BN_LOC, BN_MLT, &
     & BN_NEY, NDE_ID, NHN, HN_LST, HN_OWN, EC, NC, BNC, EMAP, NMAP, &
     & PNE, PNN, PNHE, PNBN, PMBM, PNHN, NODE_MATCH, EXCHANGE, GATHER, &
     & SORT
  !
      Use MOD_SIZES, Only: MGL, NGL
      Implicit None
# if defined (MULTIPROCESSOR)
      Include "mpif.h"
# endif
      Integer I, J, K, N1, IERR, NGLOB, NN, I1, I2, I3, NSZE, ISND, &
     & IRCV, OWNER ! ,N2,N3,CHECK!LBcleanup
      Integer, Allocatable :: NTEMP (:), NTEMP2 (:), NTEMP3 (:, :), &
     & NTEMP4 (:, :)!NP(:,:),!LBcleanup
      Integer, Allocatable :: NODEMARK (:)
      Integer :: NNODES_MINE
      Integer IGLID, TERCV !, IC1, IC2 !LBcleanup
      Integer STAT (MPI_STATUS_SIZE)
      Logical :: INELEM, ISMYN !NEW,!LBcleanup
      Real (SP), Allocatable :: TESTA (:, :), TESTB (:, :), TESTC (:, &
     & :)
      Real (SP) DIFF, VAL, NEAVE, NNAVE !,TIME_PAR!LBcleanup
      Real (SP) :: FCHECK
      Intrinsic SYSTEM_CLOCK
  !
  !
      If (MSR) WRITE (*,*) '!'
      If (MSR) WRITE (*,*) '!            SETTING UP MAPPING/MESSAGE PAS&
     &SING '
      If (MSR) WRITE (*,*) '!'
  !
  !==============================================================================|
  !   COUNT NUMBER OF LOCAL DOMAIN ELEMENTS                        NLOC             |
  !   CREATE LOCAL--> GLOBAL ELEMENT MAP                           EGID(1:NLOC)     |
  !   EGID(I) = GLOBAL ELEMENT NUMBER OF LOCAL ELEMENT I                         |
  !==============================================================================|
      Allocate (NTEMP(NGL))
      NTEMP = 0
      NLOC = 0
      Do I = 1, NGL
         If (EL_PID(I) == MYID) Then
            NLOC = NLOC + 1
            NTEMP (NLOC) = I
         End If
      End Do
      Allocate (EGID(0:NLOC))
      EGID = 0
      EGID (1:NLOC) = NTEMP (1:NLOC)
      Deallocate (NTEMP)
      If (MSR) WRITE (*,*) '!  ELEMENT MAP           :    COMPLETE'
  !==============================================================================|
  !   COUNT NUMBER OF LOCAL NODES                                  MLOC             |
  !   CREATE LOCAL--> GLOBAL NODE MAP                              NGID(1:MLOC)     |
  !   NGID(I) = GLOBAL NODE NUMBER OF LOCAL NODE I                               |
  !==============================================================================|
      Allocate (NODEMARK(MGL))
      NODEMARK = 0
      Allocate (NTEMP(MGL))
      Do I = 1, NGL
         If (EL_PID(I) == MYID) Then
            NODEMARK (NVG(I, 1:3)) = 1
         End If
      End Do
  !
      MLOC = 0
      Do I = 1, MGL
         If (NODEMARK(I) == 1) Then
            MLOC = MLOC + 1
            NTEMP (MLOC) = I
         End If
      End Do
  !   DEALLOCATE(NODEMARK)
  !
  !
      Allocate (NGID(MLOC))
      NGID (1:MLOC) = NTEMP (1:MLOC)
      Deallocate (NTEMP)
  !
      If (MSR) WRITE (*,*) '!  NODE MAP              :    COMPLETE'
  !
  !==============================================================================|
  !   DETERMINE NUMBER OF HALO ELEMENTS                            NHE           |
  !   DETERMINE GLOBAL MAPPING OF HALO ELEMENTS                    HE_LST(NHE)   |
  !   DETERMINE LOCAL MAPPING OF HALO ELEMENTS                     HE_LOC(NHE)   |
  !==============================================================================|
  !
      Allocate (NTEMP(NGL))
      Allocate (NTEMP2(NGL))
      NHE = 0
      Do I = 1, NGL
         NNODES_MINE = SUM (NODEMARK(NVG(I, 1:3)))
     !     IF(NNODES_MINE == 1 .OR. NNODES_MINE ==2)THEN
         If (NNODES_MINE > 0 .And. EL_PID(I) /= MYID) Then
            NHE = NHE + 1
            NTEMP (NHE) = I
            NTEMP2 (NHE) = EL_PID (I)
         End If
      End Do
      Deallocate (NODEMARK)
  !
      Allocate (HE_LST(NHE), HE_OWN(NHE))
      HE_LST (1:NHE) = NTEMP (1:NHE)
      HE_OWN (1:NHE) = NTEMP2 (1:NHE)
      Deallocate (NTEMP, NTEMP2)
      If (MSR) WRITE (*,*) '!  HALO ELEMENTS         :    COMPLETE'
  !
  !==============================================================================|
  !   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL ONLY)     ELID          |
  !   ELID(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                             |
  !             RETURNS 0 IF I IS NOT AN INTERIOR ELEMENT                        |
  !                                                                              |
  !   DETERMINE LOCAL MAPPING OF GLOBAL ELMNTS (INTERNAL + HALO)   ELID_X        |
  !   ELID_X(I) = LOCAL ELEMENT ID OF GLOBAL ELEMENT I                           |
  !             RETURNS 0 IF I IS NOT AN INTERIOR OR HALO ELEMENT                |
  !==============================================================================|
  !
      Allocate (ELID(0:NGL))
      ELID = 0
      Allocate (ELID_X(0:NGL))
      ELID = 0
  !
      Do I = 1, NLOC
         ELID (EGID(I)) = I
      End Do
      ELID_X = ELID
  !
      Do I = 1, NHE
         ELID_X (HE_LST(I)) = I + NLOC
      End Do
  !
  !
  !==============================================================================|
  !   DETERMINE NUMBER OF INTERNAL BOUNDARY NODES             NBN                |
  !   DETERMINE GLOBAL MAPPING OF INTERNAL BOUNDARY NODES     BN_LST(NBN)        |
  !   DETERMINE LOCAL MAPPING OF INTERNAL BOUNDARY NODES      BN_LOC(NBN)        |
  !   DETERMINE MULTIPLICITY (HOW MANY OWNERS) OF IBNs        BN_MLT(NBN)        |
  !   MARK OTHER OWNERS OF IBNs                               BN_NEY(NBN,NPROCS) |
  !==============================================================================|
  !
      Allocate (NTEMP(MLOC))
      NTEMP = 0
      Allocate (NTEMP3(MLOC, NPROCS))
      NTEMP3 = 0
      Allocate (NDE_ID(MLOC))
      NDE_ID = 0
  !
      Do I = 1, MLOC
         Do J = 1, NGL
            NGLOB = NGID (I)
            FCHECK = DBLE (NVG(J, 1)-NGLOB) * DBLE (NVG(J, 2)-NGLOB) * &
           & DBLE (NVG(J, 3)-NGLOB)
            If (Abs(FCHECK) < 1.0 .And. MYID /= EL_PID(J)) Then
               NTEMP (I) = 1
               NTEMP3 (I, EL_PID(J)) = 1
               NDE_ID (I) = 1
            End If
         End Do
      End Do
  !
      NBN = SUM (NTEMP)
      Allocate (BN_LST(NBN))
      Allocate (BN_MLT(NBN))
      Allocate (BN_LOC(NBN))
      Allocate (BN_NEY(NBN, NPROCS))
      BN_NEY = 0
      BN_MLT = 0
      BN_LST = 0
      NN = 1
      Do I = 1, MLOC
         If (NTEMP(I) == 1) Then
            BN_LST (NN) = NGID (I)
            BN_LOC (NN) = I
            BN_NEY (NN, MYID) = 1
            Do J = 1, NPROCS
               If (NTEMP3(I, J) == 1) BN_NEY (NN, J) = 1
            End Do
            BN_MLT (NN) = SUM (NTEMP3(I, :))
            NN = NN + 1
         End If
      End Do
      BN_MLT = BN_MLT + 1
      MX_MLT = 1
      If (NBN > 0) MX_MLT = MAXVAL (BN_MLT)
      Deallocate (NTEMP, NTEMP3)
      If (MSR) WRITE (*,*) '!  BOUNDARY NODES        :    COMPLETE'
  !
  !==============================================================================|
  !   DETERMINE NUMBER OF HALO NODES                               NHN           |
  !   DETERMINE GLOBAL NUMBERING OF HALO NODES                     HN_LST(NHN)   |
  !   DETERMINE PRIMARY OWNER OF HALO NODE                         HN_OWN(NHN)   !
  !==============================================================================|
  !
      Allocate (NTEMP(MGL))
      NTEMP = 0
      NHN = 0
      Do I = 1, MGL
         Do J = 1, NHE
            I1 = NVG (HE_LST(J), 1)
            I2 = NVG (HE_LST(J), 2)
            I3 = NVG (HE_LST(J), 3)
            If (I1 == I .Or. I2 == I .Or. I3 == I) Then
               INELEM = .True.
               ISMYN = .False.
               Do K = 1, MLOC
                  NGLOB = NGID (K)
                  If (I == NGLOB) ISMYN = .True.
               End Do
               If (INELEM .And. .Not. ISMYN) Then
                  NHN = NHN + 1
                  NTEMP (NHN) = I
                  Exit
               End If
            End If
         End Do
      End Do
      Allocate (HN_LST(NHN))
      HN_LST = NTEMP (1:NHN)
      Deallocate (NTEMP)
  !
      Allocate (HN_OWN(NHN))
      HN_OWN = 0
      Do I = 1, NHN
         IGLID = HN_LST (I)
         Do J = 1, NHE
            I1 = NVG (HE_LST(J), 1)
            I2 = NVG (HE_LST(J), 2)
            I3 = NVG (HE_LST(J), 3)
            If (I1 == IGLID .Or. I2 == IGLID .Or. I3 == IGLID) Then
               HN_OWN (I) = HE_OWN (J)
               Exit
            End If
         End Do
         If (HN_OWN(I) == 0) Then
            Write (*,*) 'NO OWNER FOUND FOR HALO NODE', I, 'OF DOMAIN', &
           & MYID
         End If
      End Do
  !
      If (MSR) WRITE (*,*) '!  HALO NODES            :    COMPLETE'
  !
  !==============================================================================|
  !   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL ONLY)     NLID          |
  !   DETERMINE LOCAL MAPPING OF GLOBAL NODES  (INTERNAL + HALO)   NLID_X        |
  !==============================================================================|
  !
      Allocate (NLID(0:MGL))
      NLID = 0
      Allocate (NLID_X(0:MGL))
      NLID_X = 0
  !
      Do I = 1, MLOC
         NLID (NGID(I)) = I
      End Do
      NLID_X = NLID
  !
      Do I = 1, NHN
         NLID_X (HN_LST(I)) = I + MLOC
      End Do
  !
  !
  !==============================================================================|
  !   SEND INFORMATION TO PROCESSOR ZERO FOR STATISTICS OUTPUT                   |
  !   PNE  :: NUMBER OF ELEMENTS IN EACH DOMAIN                                  |
  !   PNN  :: NUMBER OF NODES IN EACH DOMAIN                                     |
  !   PNHE :: NUMBER OF HALO ELEMENTS IN EACH DOMAIN                             |
  !   PNBN :: NUMBER OF BOUNDARY NODES IN EACH DOMAIN                            |
  !   PMBM :: MAXIMUM MULTIPLICITY OF BOUNDARY NODES IN EACH DOMAIN              |
  !   PNHN :: NUMBER OF HALO NODES IN EACH DOMAIN                                |
  !==============================================================================|
  !
      Allocate (PNE(NPROCS))
      Allocate (PNN(NPROCS))
      Allocate (PNHE(NPROCS))
      Allocate (PNBN(NPROCS))
      Allocate (PMBM(NPROCS))
      Allocate (PNHN(NPROCS))
  !
      Call MPI_GATHER (NLOC, 1, MPI_INTEGER, PNE, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      Call MPI_GATHER (MLOC, 1, MPI_INTEGER, PNN, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      If (MSR) Then
         Write (*,*) '         PROC     #ELEMENTS   IMBALANCE      #NOD&
        &ES  IMBALANCE'
         NEAVE = FLOAT (NGL) / FLOAT (NPROCS)
         NNAVE = FLOAT (SUM(PNN(1:NPROCS))) / FLOAT (NPROCS)!average number of local nodes per proc
         Do I = 1, NPROCS
            Write (*, '(I12,I14,F12.6,I14,F12.6)') I, PNE (I), FLOAT &
           & (PNE(I)) / NEAVE, PNN (I), FLOAT (PNN(I)) / NNAVE
         End Do
      End If
      Call MPI_GATHER (NHE, 1, MPI_INTEGER, PNHE, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      Call MPI_GATHER (NBN, 1, MPI_INTEGER, PNBN, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      Call MPI_GATHER (MX_MLT, 1, MPI_INTEGER, PMBM, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      Call MPI_GATHER (NHN, 1, MPI_INTEGER, PNHN, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
      If (MSR) WRITE (*,*) '         PROC     #HALO ELMNTS   #BNDRY NOD&
     &ES  MAX MLTPLCTY  #HALO NODES'
      Do I = 1, NPROCS
         If (MSR) WRITE (*, '(I12,4I14)') I, PNHE (I), PNBN (I), PMBM &
        & (I), PNHN (I)
      End Do
  !
  !
  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)       |
  !   ELEMENT MAP :: EMAP(NPROCS)                                                |
  !     EMAP(1-->NPROCS)%NSIZE  :: NUMBER OF ELEMENTS IN EACH DOM                |
  !     EMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
  !==============================================================================|
  !
      Allocate (EMAP(NPROCS))
      EMAP(:)%NSIZE = 0
      If (MSR) Then
     !
     !--Determine Number of Elements for Each Processor
         Do I = 1, NGL
            EMAP (EL_PID(I))%NSIZE = EMAP (EL_PID(I))%NSIZE + 1
         End Do
     !
     !--Allocate Mapping Array for Each Processor
         Do I = 1, NPROCS
            Allocate (EMAP(I)%LOC_2_GL(0:EMAP(I)%NSIZE))
            EMAP(I)%LOC_2_GL(0) = 0
         End Do
     !
     !--Construct Mapping Array for Each Processor
         Allocate (NTEMP(NPROCS))
         NTEMP = 0
         Do I = 1, NGL
            I1 = EL_PID (I)
            NTEMP (I1) = NTEMP (I1) + 1
            EMAP(I1)%LOC_2_GL(NTEMP(I1)) = I
         End Do
         Deallocate (NTEMP)
      End If
  !
  !==============================================================================|
  !   TEST ELEMENT MAPPING USING GATHER                                          |
  !==============================================================================|
  !
      If (MSR) Then
         Allocate (TESTA(NLOC, KB))
         Allocate (TESTB(NGL, KB))
      Else
         Allocate (TESTA(NLOC, KB))
         Allocate (TESTB(1, 1))
      End If
  !
      TESTA = FLOAT (MYID)
      Call GATHER (LBOUND(TESTA, 1), UBOUND(TESTA, 1), NLOC, NGL, KB, &
     & MYID, NPROCS, EMAP, TESTA, TESTB)
      If (MSR) Then
         Do J = 1, KB
            Do I = 1, NGL
               If (TESTB(I, J) /= EL_PID(I)) Then
                  Write (*,*) 'GLOBAL ARRAY NOT CORRECT', I, J, TESTB &
                 & (I, J), EL_PID (I)
               End If
            End Do
         End Do
      End If
      Deallocate (TESTA, TESTB)
  !
  !
      If (MSR) WRITE (*,*) '!  ELEMENT GATHER TEST   :    PASSED'
  !
  !
  !==============================================================================|
  !   SET UP NODE MAPPING FOR GLOBAL REASSEMBLY OF ARRAYS (MASTER ONLY)          |
  !   NODAL   MAP :: NMAP(NPROCS)                                                |
  !     NMAP(1-->NPROCS)%NSIZE  :: NUMBER OF NODES IN EACH DOM                   |
  !     NMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM  |
  !==============================================================================|
  !
      Allocate (NMAP(NPROCS))
      NMAP(:)%NSIZE = 0
      If (MSR) Then
         NMAP(:)%NSIZE = PNN (:)
     !
         Do I = 1, NPROCS
            Allocate (NMAP(I)%LOC_2_GL(0:NMAP(I)%NSIZE))
            NMAP(I)%LOC_2_GL(0) = 0
         End Do
      End If
  !
      If ( .Not. MSR) Then
         ISND = MYID + 1000
         Call MPI_SEND (NGID, MLOC, MPI_INTEGER, 0, ISND, &
        & MPI_COMM_WORLD, IERR)
      End If
  !
      If (MSR) Then
         NMAP(1)%LOC_2_GL(1:MLOC) = NGID (1:MLOC)
         Do I = 2, NPROCS
            NSZE = NMAP(I)%NSIZE
            Allocate (NTEMP(NSZE))
            IRCV = I + 1000
            Call MPI_RECV (NTEMP, NSZE, MPI_INTEGER, I-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
            NMAP(I)%LOC_2_GL(1:NSZE) = NTEMP
            Deallocate (NTEMP)
         End Do
      End If
  !WL   debug3602 = -1 ! not found yet
  !WL   do j=1,MLOC
  !WL      if(NGID(j) .eq. 3602) debug3602=j
  !WL   enddo
  !
  !==============================================================================|
  !   TEST NODE MAPPING USING GATHER                                             |
  !==============================================================================|
  !
      Allocate (TESTA(MLOC, KB))
      Allocate (TESTB(MGL, KB))
  !
  !
      TESTA = 100.0
      Call GATHER (LBOUND(TESTA, 1), UBOUND(TESTA, 1), MLOC, MGL, KB, &
     & MYID, NPROCS, NMAP, TESTA, TESTB)
      If (MSR) Then
         Do J = 1, KB
            Do I = 1, MGL
               If (TESTB(I, J) /= 100.0) Then
                  Write (*,*) 'GLOBAL ARRAY NOT CORRECT', I, J, TESTB &
                 & (I, J)
               End If
            End Do
         End Do
      End If
      Deallocate (TESTA, TESTB)
  !
      If (MSR) WRITE (*,*) '!  NODE GATHER TEST      :    PASSED'
      Call MPI_BARRIER (MPI_COMM_WORLD, IERR)
  !
  !
  !
  !==============================================================================|
  !   SET UP COMMUNICATION: HALO ELEMENTS                                        |
  !==============================================================================|
  !
      Allocate (EC(NPROCS))
  !
  !
      EC(:)%NRCV = 0
      Allocate (NTEMP3(NPROCS, NGL))
  !
  !--Set up Recieve Information First
      Do I = 1, NHE
         OWNER = HE_OWN (I)
         EC(OWNER)%NRCV = EC(OWNER)%NRCV + 1
         NTEMP3 (OWNER, EC(OWNER)%NRCV) = I + NLOC
      End Do
  !
      NN = 0
      Do I = 1, NPROCS
         EC(I)%RCPT = NN
         NSZE = EC(I)%NRCV
         NN = NN + NSZE
     !     ALLOCATE(EC(I)%RCVP(NSZE))
         If (NSZE > 0) ALLOCATE (EC(I)%RCVP(NSZE))
         If (NSZE == 0) ALLOCATE (EC(I)%RCVP(1))
         If (NSZE > 0) EC(I)%RCVP = NTEMP3 (I, 1:NSZE)
      End Do
      Deallocate (NTEMP3)
      TERCV = SUM (EC(:)%NRCV)
      If (TERCV /= NHE) Then
         Write (*,*) 'TOTAL NUMBER OF ELEMENTS SET UP TO RECEIVE DATA'
         Write (*,*) 'NOT EQUAL TO TOTAL NUMBER OF HALO ELEMENTS'
         Write (*,*) TERCV, NHE
         Call MPI_FINALIZE (IERR)
         Call PSTOP
      End If
  !
  !--Set up Send
      EC(:)%NSND = 0
  !
      Do I = 1, NPROCS
         If (MYID /= I) Then
            ISND = MYID + 1000 * I
            Call MPI_SEND (EC(I)%NRCV, 1, MPI_INTEGER, I-1, ISND, &
           & MPI_COMM_WORLD, IERR)
         End If
      End Do
      Do I = 1, NPROCS
         If (MYID /= I) Then
            IRCV = 1000 * MYID + I
            Call MPI_RECV (EC(I)%NSND, 1, MPI_INTEGER, I-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
         End If
     !   ALLOCATE(EC(I)%SNDP(EC(I)%NSND))
         If (EC(I)%NSND > 0) ALLOCATE (EC(I)%SNDP(EC(I)%NSND))
         If (EC(I)%NSND == 0) ALLOCATE (EC(I)%SNDP(1))
      End Do
      Do I = 1, NPROCS
         If (MYID /= I .And. EC(I)%NRCV > 0) Then
            Allocate (NTEMP(EC(I)%NRCV))
            Do J = 1, EC(I)%NRCV
               NTEMP (J) = HE_LST (EC(I)%RCVP(J)-NLOC)
            End Do
            ISND = MYID + 2000 * I
            Call MPI_SEND (NTEMP, EC(I)%NRCV, MPI_INTEGER, I-1, ISND, &
           & MPI_COMM_WORLD, IERR)
            Deallocate (NTEMP)
         End If
      End Do
      Do I = 1, NPROCS
         If (MYID /= I .And. EC(I)%NSND > 0) Then
            Allocate (NTEMP(EC(I)%NSND))
            IRCV = 2000 * MYID + I
            Call MPI_RECV (NTEMP, EC(I)%NSND, MPI_INTEGER, I-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
            N1 = 0
            Do J = 1, EC(I)%NSND
               Do K = 1, NLOC
                  If (EGID(K) == NTEMP(J)) Then
                     N1 = N1 + 1
                     EC(I)%SNDP(N1) = K
                  End If
               End Do
            End Do
            Deallocate (NTEMP)
         End If
      End Do
  !
  !
  !==============================================================================|
  !   CHECK ELEMENT COMMUNICATION                                                |
  !==============================================================================|
  !--Set Up Testing Array
      Allocate (TESTA(0:NLOC+NHE, KB))
      TESTA = 0.0
      Allocate (TESTB(0:NLOC+NHE, KB))
      TESTB = 0.0
      Allocate (TESTC(0:NLOC+NHE, KB))
      TESTC = 0.0
  !
  !--Initialize with Function of Global Element Index in Interior, Zero in Halo
      Do I = 1, KB
         TESTA (1:NLOC, I) = (FLOAT(EGID(1:NLOC))+10000*I)
         TESTB (1:NLOC, I) = - (FLOAT(EGID(1:NLOC))+10000*I)
         TESTC (1:NLOC, I) = (FLOAT(EGID(1:NLOC))+10000*I) + .5
      End Do
  !
  !
  !--Perform Communication
  !
      Call EXCHANGE (EC, NLOC+NHE, KB, MYID, NPROCS, TESTA, TESTB, &
     & TESTC)
  !
  !--Check Results For Interior on Array A
      Do I = 1, NLOC
         Do J = 1, KB
            VAL = FLOAT (EGID(I)+10000*J)
            DIFF = Abs (TESTA(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL ELEMENT VALUE CHANGED on A', I, &
              & J, TESTA (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array A
      Do I = 1, NHE
         Do J = 1, KB
            VAL = FLOAT (HE_LST(I)+10000*J)
            DIFF = Abs (TESTA(I+NLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO ELEMENT CHANGED on A', I, J, TESTA &
              & (I+NLOC, J), VAL, I + NLOC
               Call PSTOP
            End If
         End Do
      End Do
  !
  !
  !--Check Results For Interior on Array B
      Do I = 1, NLOC
         Do J = 1, KB
            VAL = - FLOAT (EGID(I)+10000*J)
            DIFF = Abs (TESTB(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL ELEMENT VALUE CHANGED on B', I, &
              & J, TESTB (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array B
      Do I = 1, NHE
         Do J = 1, KB
            VAL = - FLOAT (HE_LST(I)+10000*J)
            DIFF = Abs (TESTB(I+NLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO ELEMENT CHANGED on B', I, J, TESTB &
              & (I+NLOC, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Interior on Array C
      Do I = 1, NLOC
         Do J = 1, KB
            VAL = FLOAT (EGID(I)+10000*J) + .5
            DIFF = Abs (TESTC(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL ELEMENT VALUE CHANGED on C', I, &
              & J, TESTC (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array C
      Do I = 1, NHE
         Do J = 1, KB
            VAL = FLOAT (HE_LST(I)+10000*J) + .5
            DIFF = Abs (TESTC(I+NLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO ELEMENT CHANGED on C', I, J, TESTC &
              & (I+NLOC, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
      Deallocate (TESTA, TESTB, TESTC)
  !
      If (MSR) WRITE (*,*) '!  ELEMENT COMM TEST     :    PASSED'
      Call MPI_BARRIER (MPI_COMM_WORLD, IERR)
  !
  !==============================================================================|
  !   SET UP COMMUNICATION: BOUNDARY NODES                                       |
  !==============================================================================|
  !
      Allocate (NTEMP3(NBN, NPROCS), NTEMP4(NBN, NPROCS))
      Allocate (BNC(NPROCS))
      BNC%NRCV = 0
      BNC%NSND = 0
  !
  !--Count Receives from Each Processor
      Do I = 1, NBN
         Do J = 1, NPROCS
            If (BN_NEY(I, J) == 1 .And. MYID /= J) Then
               BNC(J)%NRCV = BNC(J)%NRCV + 1
               NTEMP3 (BNC(J)%NRCV, J) = BN_LOC (I)
               NTEMP4 (BNC(J)%NRCV, J) = BN_MLT (I)
            End If
         End Do
      End Do
  !
  !--Sent up Indices for Receive Buffer
      NN = 0
      Do I = 1, NPROCS
         BNC(I)%RCPT = NN
         NSZE = BNC(I)%NRCV
         NN = NN + NSZE
      End Do
  !
  !--Sends = Receives for Boundary Nodes
      BNC%NSND = BNC%NRCV
  !
  !--Array of Receives/Sends/Multiplicities for Each Processor
      Do I = 1, NPROCS
         NSZE = BNC(I)%NRCV
         Allocate (BNC(I)%RCVP(NSZE))
         Allocate (BNC(I)%SNDP(NSZE))
         Allocate (BNC(I)%MLTP(NSZE))
         BNC(I)%RCVP = NTEMP3 (1:NSZE, I)
         BNC(I)%SNDP = NTEMP3 (1:NSZE, I)
         BNC(I)%MLTP = NTEMP4 (1:NSZE, I)
      End Do
      Deallocate (NTEMP3)
  !
  !
  !--Sort Send/Recive Arrays According to Global Identity to Match Order
  !
      Do I = 1, NPROCS
         If (BNC(I)%NRCV > 0) Then
            NSZE = BNC(I)%NRCV
            Allocate (NTEMP(NSZE), NTEMP2(NSZE))
            NTEMP (:) = NGID (BNC(I)%RCVP(:))
            Call SORT (NTEMP, NTEMP2, NSZE)
            NTEMP = BNC(I)%RCVP
            Do J = 1, NSZE
               BNC(I)%RCVP(J) = NTEMP (NTEMP2(J))
            End Do
            NTEMP = BNC(I)%MLTP
            Do J = 1, NSZE
               BNC(I)%MLTP(J) = NTEMP (NTEMP2(J))
            End Do
            BNC(I)%SNDP = BNC(I)%RCVP
            Deallocate (NTEMP, NTEMP2)
         End If
      End Do
  !
  !
  !==============================================================================|
  !   TEST COMMUNICATION: HALO NODES                                             |
  !==============================================================================|
  !
  !--Set Up Testing Array
      Allocate (TESTA(0:MLOC+NHN, KB))
      Allocate (TESTB(0:MLOC+NHN, KB))
      Allocate (TESTC(0:MLOC+NHN, KB))
      TESTA = FLOAT (MYID)
      TESTB = FLOAT (MYID) * 100.0
      TESTC = FLOAT (MYID) * 10000.0
      Do J = 1, KB
         Do I = 1, NBN
            TESTA (BN_LOC(I), J) = J
            TESTB (BN_LOC(I), J) = 100 * J
            TESTC (BN_LOC(I), J) = 10000 * J
         End Do
      End Do
  !
  !--Perform Communication
  !
      Call NODE_MATCH (1, NBN, BN_MLT, BN_LOC, BNC, MLOC+NHN, KB, MYID, &
     & NPROCS, TESTA, TESTB, TESTC)
  !
  !
  !--Check Results For Interior and Boundary on Array A
      Do J = 1, KB
         Do I = 1, MLOC
            If (NDE_ID(I) == 0) Then
               If (Abs(TESTA(I, J)-FLOAT(MYID)) > .01) Then
                  Write (*,*) 'INTERNAL NODE CHANGED', I, J, TESTA (I, &
                 & J), FLOAT (MYID)
                  Call PSTOP
               End If
            Else
               If (Abs(TESTA(I, J)-J) > .01) Then
                  Write (*,*) 'BOUNDARY PROBLEM', MYID, I, J, TESTA (I, &
                 & J)
                  Write (*,*) 'GLOBAL NODE: ', NGID (I)
                  Call PSTOP
               End If
            End If
         End Do
      End Do
  !
  !--Check Results For Interior and Boundary on Array B
      Do J = 1, KB
         Do I = 1, MLOC
            If (NDE_ID(I) == 0) Then
               If (Abs(TESTB(I, J)-FLOAT(MYID)*100.0) > .01) Then
                  Write (*,*) 'INTERNAL NODE CHANGED', I, J, TESTB (I, &
                 & J), FLOAT (MYID)
                  Call PSTOP
               End If
            Else
               If (Abs(TESTB(I, J)-100*J) > .01) Then
                  Write (*,*) 'BOUNDARY PROBLEM', MYID, I, J, TESTB (I, &
                 & J)
                  Call PSTOP
               End If
            End If
         End Do
      End Do
  !
  !--Check Results For Interior and Boundary on Array C
      Do J = 1, KB
         Do I = 1, MLOC
            If (NDE_ID(I) == 0) Then
               If (Abs(TESTC(I, J)-FLOAT(MYID)*10000.0) > .01) Then
                  Write (*,*) 'INTERNAL NODE CHANGED', I, J, TESTC (I, &
                 & J), FLOAT (MYID)
                  Call PSTOP
               End If
            Else
               If (Abs(TESTC(I, J)-10000.0*J) > .01) Then
                  Write (*,*) 'BOUNDARY PROBLEM', MYID, I, J, TESTC (I, &
                 & J)
                  Call PSTOP
               End If
            End If
         End Do
      End Do
  !
      Deallocate (TESTA, TESTB, TESTC)
  !
  !
      If (MSR) WRITE (*,*) '!  BNDRY NODE COMM TEST  :    PASSED'
  !
  !==============================================================================|
  !   SET UP COMMUNICATION: HALO NODES                                           |
  !==============================================================================|
  ! NHN,HN_LST --HN_OWN (owner of halo node i, like he_own)
  !
      Allocate (NC(NPROCS))
      NC(:)%NRCV = 0
  !
      Allocate (NTEMP3(NPROCS, MGL))
  !
  !--Set up Recieve Information First
      Do I = 1, NHN
         OWNER = HN_OWN (I)
         NC(OWNER)%NRCV = NC(OWNER)%NRCV + 1
         NTEMP3 (OWNER, NC(OWNER)%NRCV) = I + MLOC
      End Do
  !
      NN = 0
      Do I = 1, NPROCS
         NC(I)%RCPT = NN
         NSZE = NC(I)%NRCV
         NN = NN + NSZE
         Allocate (NC(I)%RCVP(NSZE))
         NC(I)%RCVP = NTEMP3 (I, 1:NSZE)
      End Do
      Deallocate (NTEMP3)
      TERCV = SUM (NC(:)%NRCV)
      If (TERCV /= NHN) Then
         Write (*,*) 'TOTAL NUMBER OF NODES SET UP TO RECEIVE DATA'
         Write (*,*) 'NOT EQUAL TO TOTAL NUMBER OF HALO NODES'
         Write (*,*) TERCV, NHN
         Call MPI_FINALIZE (IERR)
         Call PSTOP
      End If
  !
  !--Set up Send
      NC(:)%NSND = 0
  !
      Do I = 1, NPROCS
         If (MYID /= I) Then
            ISND = MYID + 1000 * I
            Call MPI_SEND (NC(I)%NRCV, 1, MPI_INTEGER, I-1, ISND, &
           & MPI_COMM_WORLD, IERR)
         End If
      End Do
      Do I = 1, NPROCS
         If (MYID /= I) Then
            IRCV = 1000 * MYID + I
            Call MPI_RECV (NC(I)%NSND, 1, MPI_INTEGER, I-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
         End If
         Allocate (NC(I)%SNDP(NC(I)%NSND))
      End Do
  !
      Do I = 1, NPROCS
         If (MYID /= I .And. NC(I)%NRCV > 0) Then
            Allocate (NTEMP(NC(I)%NRCV))
            Do J = 1, NC(I)%NRCV
               NTEMP (J) = HN_LST (NC(I)%RCVP(J)-MLOC)
            End Do
            ISND = MYID + 1000 * I
            Call MPI_SEND (NTEMP, NC(I)%NRCV, MPI_INTEGER, I-1, ISND, &
           & MPI_COMM_WORLD, IERR)
            Deallocate (NTEMP)
         End If
      End Do
  !
      Do I = 1, NPROCS
         If (MYID /= I .And. NC(I)%NSND > 0) Then
            Allocate (NTEMP(NC(I)%NSND))
            IRCV = 1000 * MYID + I
            Call MPI_RECV (NTEMP, NC(I)%NSND, MPI_INTEGER, I-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
            N1 = 0
            Do J = 1, NC(I)%NSND
               Do K = 1, MLOC
                  If (NGID(K) == NTEMP(J)) Then
                     N1 = N1 + 1
                     NC(I)%SNDP(N1) = K
                  End If
               End Do
            End Do
            Deallocate (NTEMP)
         End If
      End Do
  !
  !==============================================================================|
  !   CHECK HALO NODE COMMUNICATION                                              |
  !==============================================================================|
  !
  !--Set Up Testing Array
      Allocate (TESTA(0:MLOC+NHN, KB))
      TESTA = 0.0
      Allocate (TESTB(0:MLOC+NHN, KB))
      TESTB = 0.0
      Allocate (TESTC(0:MLOC+NHN, KB))
      TESTC = 0.0
  !
  !--Initialize with Global Nodal ID in Interior, Zero in Halo
      Do I = 1, KB
         TESTA (1:MLOC, I) = FLOAT (NGID) + 10000 * I
         TESTB (1:MLOC, I) = - (FLOAT(NGID)+10000*I)
         TESTC (1:MLOC, I) = (FLOAT(NGID)+10000*I) + .50
      End Do
  !
  !--Perform Communication
      Call EXCHANGE (NC, MLOC+NHN, KB, MYID, NPROCS, TESTA, TESTB, &
     & TESTC)
  !--Check Results For Interior on Array A
      Do I = 1, MLOC
         Do J = 1, KB
            VAL = FLOAT (NGID(I)+10000*J)
            DIFF = Abs (TESTA(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL NODE VALUE CHANGED on A', I, J, &
              & TESTA (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array A
      Do I = 1, NHN
         Do J = 1, KB
            VAL = FLOAT (HN_LST(I)+10000*J)
            DIFF = Abs (TESTA(I+MLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO NODE CHANGED on A', I, J, TESTA &
              & (I+MLOC, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !
  !--Check Results For Interior on Array B
      Do I = 1, MLOC
         Do J = 1, KB
            VAL = - FLOAT (NGID(I)+10000*J)
            DIFF = Abs (TESTB(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL NODE VALUE CHANGED on B', I, J, &
              & TESTB (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array B
      Do I = 1, NHN
         Do J = 1, KB
            VAL = - FLOAT (HN_LST(I)+10000*J)
            DIFF = Abs (TESTB(I+MLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO NODE CHANGED on B', I, J, TESTB &
              & (I+MLOC, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Interior on Array C
      Do I = 1, MLOC
         Do J = 1, KB
            VAL = FLOAT (NGID(I)+10000*J) + .5
            DIFF = Abs (TESTC(I, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'INTERNAL NODE VALUE CHANGED on C', I, J, &
              & TESTC (I, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
  !--Check Results For Halo on Array C
      Do I = 1, NHN
         Do J = 1, KB
            VAL = FLOAT (HN_LST(I)+10000*J) + .5
            DIFF = Abs (TESTC(I+MLOC, J)-VAL)
            If (DIFF > .01) Then
               Print *, MYID, 'HALO ELEMENT CHANGED on C', I, J, TESTC &
              & (I+MLOC, J), VAL
               Call PSTOP
            End If
         End Do
      End Do
  !
      Deallocate (TESTA, TESTB, TESTC)
  !
      If (MSR) WRITE (*,*) '!  NODE COMM TEST        :    PASSED'
  !
  !==============================================================================|
  !   SET TOTAL FOR EACH PROCESSOR                                               |
  !   NTLOC:   TOTAL NUMBER OF ELEMENTS = INTERIOR (NLOC) + HALO (NHE)                 |
  !   MTLOC:   TOTAL NUMBER OF NODES    = INTERIOR (MLOC) + HALO (NHN)                 |
  !==============================================================================|
      NTLOC = NLOC + NHE
      MTLOC = MLOC + NHN
  !
      Return
# endif
End Subroutine GENMAP
!==============================================================================!
