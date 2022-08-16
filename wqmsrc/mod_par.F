!mod_par.F
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
Module MOD_PAR
# if defined (MULTIPROCESSOR)
      Use MOD_PREC, Only: MPI_F, SP
      Use MOD_TYPES, Only: GMAP, COMM
	  USE MOD_FILEINFO, ONLY: IREPORT
      Implicit None
      Include "mpif.h"
      Save
  !
  !--Global Information
  !
      Integer, Allocatable :: EL_PID (:)!!PROCESSOR OWNER OF GLOBAL ELEMENT
      Integer, Allocatable :: ELID (:)!!LOCAL VALUE OF GLOBAL ELEMENT
      Integer, Allocatable :: NLID (:)!!LOCAL VALUE OF GLOBAL NODE
      Integer, Allocatable :: ELID_X (:)!!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
      Integer, Allocatable :: NLID_X (:)!!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS
  !
  !--Global to Local Boundary Condition Transforms
  !
      Integer, Allocatable :: OBN_GL2LOC (:)!!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
  !
  !--Internal Information (Local)
  !
      Integer, Allocatable :: EGID (:)!!GLOBAL ID OF LOCAL ELEMENT
      Integer, Allocatable :: NGID (:)!!GLOBAL ID OF LOCAL NODE
  !
  !--Boundary Information: Halo Elements
  !
      Integer :: NHE !!NUMBER OF HALO ELEMENTS
      Integer, Allocatable :: HE_LST (:)!!GLOBAL IDENTITIES OF HALO ELEMENTS
      Integer, Allocatable :: HE_OWN (:)!!OWNER OF HALO ELEMENTS
  !
  !
  !--Boundary Information: Internal Boundary Nodes
  !
      Integer :: NBN !!NUMBER OF BOUNDARY NODES
      Integer :: MX_MLT !!MAX MULTIPLICITY OF BOUNDARY NODES
      Integer, Allocatable, Dimension (:) :: BN_LST !!GLOBAL IDENTITY OF BOUNDARY NODES
      Integer, Allocatable, Dimension (:) :: BN_LOC !!LOCAL IDENTITY OF BOUNDARY NODES
      Integer, Allocatable, Dimension (:) :: BN_MLT !!MULTIPLICITY OF BOUNDARY NODES
      Integer, Allocatable, Dimension (:, :) :: BN_NEY !!NODE OWNER LIST
      Integer, Allocatable, Dimension (:) :: NDE_ID !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY
!
  !
  !--Boundary Information: Halo Nodes
  !
      Integer :: NHN !!NUMBER OF HALO NODES
      Integer, Allocatable, Dimension (:) :: HN_LST !!LIST OF HALO NODES
      Integer, Allocatable, Dimension (:) :: HN_OWN !!PRIMARY OWNER OF HALO NODES
!
  !
  !--Communication Objects    [SIZE: NPROCS]
  !
      Type (COMM), Allocatable, Dimension (:) :: EC, NC, BNC
  !
  !--Maps for Global Array Reconstruction   [SIZE: NPROCS]
  !
      Type (GMAP), Allocatable, Dimension (:) :: EMAP, NMAP
  !
  !--Statistics Calculation (Valid in MASTER Only)   [SIZE: NPROCS]
  !
      Integer, Allocatable :: PNE (:)!!NUMBER OF ELEMENTS IN EACH PROC
      Integer, Allocatable :: PNN (:)!!NUMBER OF NODES IN EACH PROC
      Integer, Allocatable :: PNHE (:)!!NUMBER OF HALO ELEMENTS IN EACH PROC
      Integer, Allocatable :: PNBN (:)!!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
      Integer, Allocatable :: PMBM (:)!!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
      Integer, Allocatable :: PNHN (:)!!NUMBER OF HALO NODES IN EACH PROC
!
  !
  !--Added FILELOCK array to lock a file from being accessed by multiple
  !  processors at the same time
  !
      Integer, Allocatable :: FILELOCK (:)
      Integer :: MYFILELOCK = 0
      Logical :: FILELOCKED = .False.
!
  !===================================================================================|
      CONTAINS !!INCLUDED SUBROUTINES FOLLOW
  !===================================================================================|
  !subroutines:
  !
  ! subroutine INIT_MPI_ENV()
  !	subroutine EXCHANGE()
  !	subroutine EXCHANGE_ALL()
  !	subroutine NODE_MATCH()
  !	subroutine NODE_MATCH_ALL()
  !	subroutine GATHER()
  !	subroutine IGATHER()
  !	subroutine SORT()
  !	subroutine PREAD()
  !	subroutine PWRITE()
  !	subroutine PPRINT()
  !	subroutine LOCK()
  !	subroutine UNLOCK()
  !	subroutine ROOT_OBTAINLOCK()
  !
  !functions:
  !	function GETLOC()
!
      Subroutine INIT_MPI_ENV (MYID, NPROCS, SERIAL, PAR, MSR)
    !===================================================================================|
    !  INITIALIZE MPI ENVIRONMENT                                                       |
    !===================================================================================|
         Integer, Intent (Out) :: MYID, NPROCS
         Logical, Intent (Out) :: SERIAL, PAR, MSR
         Integer IERR
         SERIAL = .True.
         PAR = .False.
         MSR = .True.
         Call MPI_INIT (IERR)
         Call MPI_COMM_RANK (MPI_COMM_WORLD, MYID, IERR)
         Call MPI_COMM_SIZE (MPI_COMM_WORLD, NPROCS, IERR)
         MYID = MYID + 1
         If (NPROCS > 1) SERIAL = .False.
         If (NPROCS > 1) PAR = .True.
         If (MYID /= 1) MSR = .False.
!
         If ( .Not. ALLOCATED(FILELOCK)) Then !Allocate FILELOCKS and also initialize as zero
            Allocate (FILELOCK(1:NPROCS))
            FILELOCK = 0 !zero for unlocked status
         End If
         Return
      End Subroutine INIT_MPI_ENV
!
      Subroutine EXCHANGE (CM, NT, KT, MYID, NPROCS, A, B, C)
    !===================================================================================|
    !    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
    !===================================================================================|
!
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: NT, KT, MYID, NPROCS
         Type (COMM), Intent (In) :: CM (NPROCS)
         Real (SP), Intent (Inout) :: A (0:NT, KT)
         Real (SP), Intent (Inout), Optional :: B (0:NT, KT)
         Real (SP), Intent (Inout), Optional :: C (0:NT, KT)
    !------------------------------------------------------------------------------
         Logical :: BYES, CYES
         Integer :: IREQR (NPROCS), IREQS (NPROCS)
         Real (SP), Allocatable :: RBUF (:), SBUF (:)
         Integer STAT (MPI_STATUS_SIZE), ISTATR (MPI_STATUS_SIZE, &
        & NPROCS), IERR, J, N1, N2, NCNT
         Integer I, IFROM, ITO, ISTAG, IRTAG, TRCV, TSND, NVARS, LBUF, &
        & LP, NMSG, INDX, LPROC, NSZE
    !------------------------------------------------------------------------------
!
         NVARS = 1
         BYES = .False.
         CYES = .False.
         If (PRESENT(B)) Then
            NVARS = NVARS + 1
            BYES = .True.
         End If
!
         If (PRESENT(C)) Then
            NVARS = NVARS + 1
            CYES = .True.
         End If
!
         Allocate (RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
         Allocate (SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))
!
    !===================================================================================|
    !    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
    !===================================================================================|
         TRCV = 0
         Do I = 1, NPROCS
!
            If (CM(I)%NRCV > 0) Then
               IFROM = I - 1
               IRTAG = I * 1000
               TRCV = TRCV + 1
               LP = CM(I)%RCPT * NVARS * KT + 1
               LBUF = NVARS * CM(I)%NRCV * KT
               Call MPI_IRECV (RBUF(LP), LBUF, MPI_F, IFROM, IRTAG, &
              & MPI_COMM_WORLD, IREQR(TRCV), IERR)
            End If
!
         End Do
!
    !===================================================================================|
    !    SEND DATA TO NEIGHBORS                                                         |
    !===================================================================================|
         TSND = 0
         NCNT = 0
         Do I = 1, NPROCS
            LBUF = CM(I)%NSND
            If (LBUF > 0) Then
               NSZE = LBUF * KT * NVARS
          !       ALLOCATE(SBUF(NSZE))
               N2 = NCNT
               Do J = 1, KT
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  SBUF (N1:N2) = A (CM(I)%SNDP(:), J)
                  If (BYES) Then
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = B (CM(I)%SNDP(:), J)
                  End If
                  If (CYES) Then
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = C (CM(I)%SNDP(:), J)
                  End If
               End Do
               TSND = TSND + 1
               ITO = I - 1
               ISTAG = MYID * 1000
               Call MPI_ISEND (SBUF(NCNT+1), NSZE, MPI_F, ITO, ISTAG, &
              & MPI_COMM_WORLD, IREQS(TSND), IERR)
               NCNT = NCNT + LBUF * KT * NVARS
          !       DEALLOCATE(SBUF)
            End If
         End Do
!
    !===================================================================================|
    !    LOOP OVER PROCS UNTIL A MESSAGE IS RECIEVED AND UNPACK                         |
    !===================================================================================|
         Do NMSG = 1, TRCV
            Call MPI_WAITANY (TRCV, IREQR, INDX, STAT, IERR)!,RBUF,SBUF) ! 
            LPROC = STAT (MPI_SOURCE) + 1
            LP = CM(LPROC)%RCPT * NVARS * KT
            LBUF = CM(LPROC)%NRCV
            N2 = LP
            Do J = 1, KT
               N1 = N2 + 1
               N2 = N1 + LBUF - 1
               A (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2)
               If (BYES) Then
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  B (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2)
               End If
               If (CYES) Then
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  C (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2)
               End If
            End Do
         End Do
!
    !===================================================================================|
    !    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
    !===================================================================================|
!
         Call MPI_WAITALL (TSND, IREQS, ISTATR, IERR)!,RBUF,SBUF) ! 
         Deallocate (RBUF, SBUF)
!
         Return
      End Subroutine EXCHANGE
!
  !------------------------------------------------------------------------------
  !
      Subroutine EXCHANGE_ALL (CM, NT, KT, TT, MYID, NPROCS, A)
    !===================================================================================|
    !    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
	!    Note that local array A must be from (0:NT,1:KT,1:TT)							|
    !===================================================================================|
!
         Use MOD_WQM, Only: NAC, AC
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: NT, KT, TT, MYID, NPROCS
         Type (COMM), Intent (In) :: CM (NPROCS)
         Real (SP), Intent (Inout) :: A (0:NT, KT, TT)
    !------------------------------------------------------------------------------
         Integer :: IREQR (NPROCS), IREQS (NPROCS)
         Real (SP), Allocatable :: RBUF (:), SBUF (:)
         Integer STAT (MPI_STATUS_SIZE), ISTATR (MPI_STATUS_SIZE, &
        & NPROCS), IERR, J, N1, N2, NCNT
         Integer I, IFROM, ITO, ISTAG, IRTAG, TRCV, TSND, NVARS, LBUF, &
        & LP, NMSG, INDX, LPROC, NSZE
         Integer II, JCON
    !------------------------------------------------------------------------------
!
         NVARS = NAC
!
         Allocate (RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
         Allocate (SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))
!
    !===================================================================================|
    !    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
    !===================================================================================|
         TRCV = 0
         Do I = 1, NPROCS
            If (CM(I)%NRCV > 0) Then
               IFROM = I - 1
               IRTAG = I * 1000
               TRCV = TRCV + 1
               LP = CM(I)%RCPT * NVARS * KT + 1
               LBUF = NVARS * CM(I)%NRCV * KT
               Call MPI_IRECV (RBUF(LP), LBUF, MPI_F, IFROM, IRTAG, &
              & MPI_COMM_WORLD, IREQR(TRCV), IERR)
            End If
         End Do
!
    !===================================================================================|
    !    SEND DATA TO NEIGHBORS                                                         |
    !===================================================================================|
         TSND = 0
         NCNT = 0
         Do I = 1, NPROCS
            LBUF = CM(I)%NSND
            If (LBUF > 0) Then
               NSZE = LBUF * KT * NVARS
               N2 = NCNT
               Do JCON = 1, NAC
                  II = AC (JCON)
                  Do J = 1, KT
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = A (CM(I)%SNDP(:), J, II)
                  End Do
               End Do
               TSND = TSND + 1
               ITO = I - 1
               ISTAG = MYID * 1000
               Call MPI_ISEND (SBUF(NCNT+1), NSZE, MPI_F, ITO, ISTAG, &
              & MPI_COMM_WORLD, IREQS(TSND), IERR)
               NCNT = NCNT + LBUF * KT * NVARS
            End If
         End Do
!
    !===================================================================================|
    !    LOOP OVER PROCS UNTIL A MESSAGE IS RECIEVED AND UNPACK                         |
    !===================================================================================|
         Do NMSG = 1, TRCV
            Call MPI_WAITANY (TRCV, IREQR, INDX, STAT, IERR, RBUF, &
           & SBUF)! 
            LPROC = STAT (MPI_SOURCE) + 1
            LP = CM(LPROC)%RCPT * NVARS * KT
            LBUF = CM(LPROC)%NRCV
            N2 = LP
            Do JCON = 1, NAC
               II = AC (JCON)
               Do J = 1, KT
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  A (CM(LPROC)%RCVP(:), J, II) = RBUF (N1:N2)
               End Do
            End Do
         End Do
!
    !===================================================================================|
    !    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
    !===================================================================================|
!
         Call MPI_WAITALL (TSND, IREQS, ISTATR, IERR, RBUF, SBUF)! 
         Deallocate (RBUF, SBUF)
!
         Return
      End Subroutine EXCHANGE_ALL
  !===================================================================================|
!
  !===================================================================================|
!
      Subroutine NODE_MATCH (IMATCH, NBN, BN_MLT, BN_LOC, CM, NT, KT, &
     & MYID, NPROCS, A, B, C)
    !===================================================================================|
    ! IMATCH=1:   ENFORCE AGREEMENT OF A,B,C ON BOUNDARY NODES                          |
    ! IMATCH=0:   ACCUMULATE VALUES OF A,B,C AT BOUNDARY NODES                          |
	! Note that A,B,C must of shape (0:NT,1:KT)
    !===================================================================================|
!
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: IMATCH
         Integer, Intent (In) :: NBN
         Integer, Intent (In) :: BN_MLT (NBN)
         Integer, Intent (In) :: BN_LOC (NBN)
         Integer, Intent (In) :: NT, KT, MYID, NPROCS
         Type (COMM), Intent (In) :: CM (NPROCS)
         Real (SP), Intent (Inout) :: A (0:NT, KT)
         Real (SP), Intent (Inout), Optional :: B (0:NT, KT)
         Real (SP), Intent (Inout), Optional :: C (0:NT, KT)
    !------------------------------------------------------------------------------
         Logical :: BYES, CYES
         Integer :: IREQR (NPROCS), IREQS (NPROCS)
         Real (SP), Allocatable :: RBUF (:), SBUF (:)
         Integer STAT (MPI_STATUS_SIZE), ISTATR (MPI_STATUS_SIZE, &
        & NPROCS), IERR, J, N1, N2, NCNT
         Integer I, IFROM, ITO, ISTAG, IRTAG, TRCV, TSND, NVARS, LBUF, &
        & LP, NMSG, INDX, LPROC, NSZE
    !------------------------------------------------------------------------------
!
         NVARS = 1
         BYES = .False.
         CYES = .False.
         If (PRESENT(B)) Then
            NVARS = NVARS + 1
            BYES = .True.
         End If
!
         If (PRESENT(C)) Then
            NVARS = NVARS + 1
            CYES = .True.
         End If
!
         Allocate (RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
         Allocate (SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))
!
    !===================================================================================|
    !    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
    !===================================================================================|
         TRCV = 0
         Do I = 1, NPROCS
!
            If (CM(I)%NRCV > 0) Then
               IFROM = I - 1
               IRTAG = I * 1000
               TRCV = TRCV + 1
               LP = CM(I)%RCPT * NVARS * KT + 1
               LBUF = NVARS * CM(I)%NRCV * KT
               Call MPI_IRECV (RBUF(LP), LBUF, MPI_F, IFROM, IRTAG, &
              & MPI_COMM_WORLD, IREQR(TRCV), IERR)
            End If
!
         End Do
!
    !===================================================================================|
    !    SEND DATA TO NEIGHBORS                                                         |
    !===================================================================================|
         TSND = 0
         NCNT = 0
         Do I = 1, NPROCS
            LBUF = CM(I)%NSND
            If (LBUF > 0) Then
               NSZE = LBUF * KT * NVARS
          !       ALLOCATE(SBUF(NSZE))
               N2 = NCNT
               Do J = 1, KT
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  SBUF (N1:N2) = A (CM(I)%SNDP(:), J)
                  If (BYES) Then
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = B (CM(I)%SNDP(:), J)
                  End If
                  If (CYES) Then
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = C (CM(I)%SNDP(:), J)
                  End If
               End Do
               TSND = TSND + 1
               ITO = I - 1
               ISTAG = MYID * 1000
               Call MPI_ISEND (SBUF(NCNT+1), NSZE, MPI_F, ITO, ISTAG, &
              & MPI_COMM_WORLD, IREQS(TSND), IERR)
               NCNT = NCNT + LBUF * KT * NVARS
          !       DEALLOCATE(SBUF)
            End If
         End Do
!
    !===================================================================================|
    !    LOOP OVER PROCS UNTIL A MESSAGE IS RECIEVED AND UNPACK                         |
    !===================================================================================|
         Do NMSG = 1, TRCV
            Call MPI_WAITANY (TRCV, IREQR, INDX, STAT, IERR, RBUF, &
           & SBUF)! 
            LPROC = STAT (MPI_SOURCE) + 1
            LP = CM(LPROC)%RCPT * NVARS * KT
            LBUF = CM(LPROC)%NRCV
            N2 = LP
            Do J = 1, KT
               N1 = N2 + 1
               N2 = N1 + LBUF - 1
               A (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2) + A &
              & (CM(LPROC)%RCVP(:), J)
               If (BYES) Then
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  B (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2) + B &
                 & (CM(LPROC)%RCVP(:), J)
               End If
               If (CYES) Then
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  C (CM(LPROC)%RCVP(:), J) = RBUF (N1:N2) + C &
                 & (CM(LPROC)%RCVP(:), J)
               End If
            End Do
         End Do
!
    !===================================================================================|
    !    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
    !===================================================================================|
!
         Call MPI_WAITALL (TSND, IREQS, ISTATR, IERR, RBUF, SBUF)! 
         Deallocate (RBUF, SBUF)
!
    !===================================================================================|
    !  USE MULTIPLICITY OF NODES TO COMPUTE TRUE AVERAGE VALUE                          |
    !===================================================================================|
         If (IMATCH /= 1) Return
!
         Do J = 1, KT
            Do I = 1, NBN
               A (BN_LOC(I), J) = A (BN_LOC(I), J) / FLOAT (BN_MLT(I))
               If (BYES) B (BN_LOC(I), J) = B (BN_LOC(I), J) / FLOAT &
              & (BN_MLT(I))
               If (CYES) C (BN_LOC(I), J) = C (BN_LOC(I), J) / FLOAT &
              & (BN_MLT(I))
            End Do
         End Do
!
         Return
      End Subroutine NODE_MATCH
!
  !===================================================================================|
!
!
  !===================================================================================|
!
      Subroutine NODE_MATCH_ALL (IMATCH, NBN, BN_MLT, BN_LOC, CM, NT, &
     & KT, TT, MYID, NPROCS, A)
    !===================================================================================|
    ! IMATCH=1:   ENFORCE AGREEMENT OF A,B,C ON BOUNDARY NODES                          |
    ! IMATCH=0:   ACCUMULATE VALUES OF A,B,C AT BOUNDARY NODES                          |
	! NOTE that A must be of size (0:NT,KT,TT)
    !===================================================================================|
!
         Use MOD_WQM, Only: NAC, AC
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: IMATCH
         Integer, Intent (In) :: NBN
         Integer, Intent (In) :: BN_MLT (NBN)
         Integer, Intent (In) :: BN_LOC (NBN)
         Integer, Intent (In) :: NT, KT, TT, MYID, NPROCS
         Type (COMM), Intent (In) :: CM (NPROCS)
         Real (SP), Intent (Inout) :: A (0:NT, KT, TT)
    !------------------------------------------------------------------------------
         Integer :: IREQR (NPROCS), IREQS (NPROCS)
         Real (SP), Allocatable :: RBUF (:), SBUF (:)
         Integer STAT (MPI_STATUS_SIZE), ISTATR (MPI_STATUS_SIZE, &
        & NPROCS), IERR, J, N1, N2, NCNT
         Integer I, IFROM, ITO, ISTAG, IRTAG, TRCV, TSND, NVARS, LBUF, &
        & LP, NMSG, INDX, LPROC, NSZE
         Integer II, JCON
    !------------------------------------------------------------------------------
!
         NVARS = NAC
         Allocate (RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
         Allocate (SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))
!
    !===================================================================================|
    !    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
    !===================================================================================|
         TRCV = 0
         Do I = 1, NPROCS
            If (CM(I)%NRCV > 0) Then
               IFROM = I - 1
               IRTAG = I * 1000
               TRCV = TRCV + 1
               LP = CM(I)%RCPT * NVARS * KT + 1
               LBUF = NVARS * CM(I)%NRCV * KT
               Call MPI_IRECV (RBUF(LP), LBUF, MPI_F, IFROM, IRTAG, &
              & MPI_COMM_WORLD, IREQR(TRCV), IERR)
            End If
         End Do
!
    !===================================================================================|
    !    SEND DATA TO NEIGHBORS                                                         |
    !===================================================================================|
         TSND = 0
         NCNT = 0
         Do I = 1, NPROCS
            LBUF = CM(I)%NSND
            If (LBUF > 0) Then
               NSZE = LBUF * KT * NVARS
               N2 = NCNT
               Do JCON = 1, NAC
                  II = AC (JCON)
                  Do J = 1, KT
                     N1 = N2 + 1
                     N2 = N1 + LBUF - 1
                     SBUF (N1:N2) = A (CM(I)%SNDP(:), J, II)
                  End Do
               End Do
               TSND = TSND + 1
               ITO = I - 1
               ISTAG = MYID * 1000
               Call MPI_ISEND (SBUF(NCNT+1), NSZE, MPI_F, ITO, ISTAG, &
              & MPI_COMM_WORLD, IREQS(TSND), IERR)
               NCNT = NCNT + LBUF * KT * NVARS
            End If
         End Do
!
!
    !===================================================================================|
    !    LOOP OVER PROCS UNTIL A MESSAGE IS RECIEVED AND UNPACK                         |
    !===================================================================================|
         Do NMSG = 1, TRCV
            Call MPI_WAITANY (TRCV, IREQR, INDX, STAT, IERR, RBUF, &
           & SBUF)! 
            LPROC = STAT (MPI_SOURCE) + 1
            LP = CM(LPROC)%RCPT * NVARS * KT
            LBUF = CM(LPROC)%NRCV
            N2 = LP
            Do JCON = 1, NAC
               II = AC (JCON)
               Do J = 1, KT
                  N1 = N2 + 1
                  N2 = N1 + LBUF - 1
                  A (CM(LPROC)%RCVP(:), J, II) = RBUF (N1:N2) + A &
                 & (CM(LPROC)%RCVP(:), J, II)
               End Do
            End Do
         End Do
!
    !===================================================================================|
    !    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
    !===================================================================================|
!
         Call MPI_WAITALL (TSND, IREQS, ISTATR, IERR, RBUF, SBUF)! 
         Deallocate (RBUF, SBUF)
!
    !===================================================================================|
    !  USE MULTIPLICITY OF NODES TO COMPUTE TRUE AVERAGE VALUE                          |
    !===================================================================================|
         If (IMATCH /= 1) Return
!
         Do JCON = 1, NAC
            II = AC (JCON)
            Do J = 1, KT
               Do I = 1, NBN
                  A (BN_LOC(I), J, II) = A (BN_LOC(I), J, II) / FLOAT &
                 & (BN_MLT(I))
               End Do
            End Do
         End Do
!
         Return
      End Subroutine NODE_MATCH_ALL

  !===================================================================================|
      Subroutine GATHER (N1, N2, NT, NTG, KT, MYID, NPROCS, GM, A, AG)
    !
    !===================================================================================|
    !    GATHER ARRAY/VECTOR A INTO GLOBAL ARRAY/VECTOR                                 |
    !    LOCAL ARRAYS A(N1:N2,KT) --> AG(1:NTG,KT) BY MAPPING GM                        |
	!                  with actual section of A(1:NT,1:KT) gathered  into (AG(1:NTG,KT) |
    !    UPON COMPLETION ONLY MASTER PROCESSOR HAS AG                                   |
    !===================================================================================|
    !   Make this an allocatable array-valued function GWC
    !
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: N1, N2, NT, NTG, KT, MYID, NPROCS
         Type (GMAP), Intent (In) :: GM (NPROCS)
         Real (SP), Intent (In) :: A (N1:N2, KT)!source array bounds
         Real (SP), Intent (Out) :: AG (NTG, KT)
    !------------------------------------------------------------------------------
         Real (SP), Allocatable :: SBUF (:, :), RBUF (:, :)
         Integer STAT (MPI_STATUS_SIZE), IERR, I, IP, J, ISND, IRCV, &
        & NSZE !, INC  !
    !------------------------------------------------------------------------------
!
    !===================================================================================|
    !     CONFINE ARRAY TO INTERIOR ELEMENTS/NODES                                      |
    !===================================================================================|
         If (MYID /= 1) Then
            Allocate (SBUF(NT, KT))
            Do I = 1, NT
               Do J = 1, KT
                  SBUF (I, J) = A (I, J)
               End Do
            End Do
         End If
!
    !===================================================================================|
    !     SEND INFORMATION TO MASTER PROCESSOR                                          |
    !===================================================================================|
         If (MYID /= 1) Then
            ISND = MYID + 100
            Call MPI_SEND (SBUF, NT*KT, MPI_F, 0, ISND, MPI_COMM_WORLD, &
           & IERR)
            Deallocate (SBUF)
            Return
         End If
    !===================================================================================|
    !     COPY MASTER PROC DATA INTO GLOBAL ARRAY                                       |
    !===================================================================================|
!
         Do I = 1, GM(MYID)%NSIZE
            Do J = 1, KT
               AG (GM(MYID)%LOC_2_GL(I), J) = A (I, J)
            End Do
         End Do
    !===================================================================================|
    !     RECEIVE SLAVE DATA AND UNPACK INTO GLOBAL ARRAY                               |
    !===================================================================================|
!
         Do IP = 2, NPROCS
            NSZE = (GM(IP)%NSIZE)
            Allocate (RBUF(NSZE, KT))
            IRCV = IP + 100
            Call MPI_RECV (RBUF, NSZE*KT, MPI_F, IP-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
            Do I = 1, GM(IP)%NSIZE
               Do J = 1, KT
                  AG (GM(IP)%LOC_2_GL(I), J) = RBUF (I, J)
               End Do
            End Do
            Deallocate (RBUF)
         End Do
!
         Return
      End Subroutine GATHER
  !===================================================================================|
!
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  !===================================================================================|
      Subroutine IGATHER (N1, N2, NT, NTG, KT, MYID, NPROCS, GM, A, AG)
    !===================================================================================|
    !    INTEGER GATHER							            |
    !    GATHER ARRAY/VECTOR A INTO GLOBAL ARRAY/VECTOR                                 |
    !    LOCAL ARRAYS A(0:NT,KT) --> AG(1:NTG,KT) BY MAPPING GM                         |
    !    UPON COMPLETION ONLY MASTER PROCESSOR HAS AG                                   |
    !===================================================================================|
    !   Make this an allocatable array-valued function GWC
!
         Implicit None
    !   include "mpif.h"
    !------------------------------------------------------------------------------
         Integer, Intent (In) :: N1, N2, NT, NTG, KT, MYID, NPROCS
         Type (GMAP), Intent (In) :: GM (NPROCS)
         Integer, Intent (In) :: A (N1:N2, KT)
         Integer, Intent (Out) :: AG (NTG, KT)
    !------------------------------------------------------------------------------
         Integer, Allocatable :: SBUF (:, :), RBUF (:, :)
         Integer STAT (MPI_STATUS_SIZE), IERR, I, IP, J, ISND, IRCV, &
        & NSZE !, INC !
    !------------------------------------------------------------------------------
!
    !===================================================================================|
    !     CONFINE ARRAY TO INTERIOR ELEMENTS/NODES                                      |
    !===================================================================================|
         If (MYID /= 1) Then
            Allocate (SBUF(NT, KT))
            Do I = 1, NT
               Do J = 1, KT
                  SBUF (I, J) = A (I, J)
               End Do
            End Do
         End If
!
    !===================================================================================|
    !     SEND INFORMATION TO MASTER PROCESSOR                                          |
    !===================================================================================|
         If (MYID /= 1) Then
            ISND = MYID + 100
            Call MPI_SEND (SBUF, NT*KT, MPI_INTEGER, 0, ISND, &
           & MPI_COMM_WORLD, IERR)
            Deallocate (SBUF)
            Return
         End If
    !===================================================================================|
    !     COPY MASTER PROC DATA INTO GLOBAL ARRAY                                       |
    !===================================================================================|
!
         Do I = 1, GM(MYID)%NSIZE
            Do J = 1, KT
               AG (GM(MYID)%LOC_2_GL(I), J) = A (I, J)
            End Do
         End Do
    !===================================================================================|
    !     RECEIVE SLAVE DATA AND UNPACK INTO GLOBAL ARRAY                               |
    !===================================================================================|
!
         Do IP = 2, NPROCS
            NSZE = (GM(IP)%NSIZE)
            Allocate (RBUF(NSZE, KT))
            IRCV = IP + 100
            Call MPI_RECV (RBUF, NSZE*KT, MPI_INTEGER, IP-1, IRCV, &
           & MPI_COMM_WORLD, STAT, IERR)
            Do I = 1, GM(IP)%NSIZE
               Do J = 1, KT
                  AG (GM(IP)%LOC_2_GL(I), J) = RBUF (I, J)
               End Do
            End Do
            Deallocate (RBUF)
         End Do
         Return
      End Subroutine IGATHER
  !===================================================================================|
!
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  !===================================================================================|
!
      Subroutine SORT (N, M, SZE)
    !===================================================================================|
    !    SORT ELEMENTS IN N AND RETURN ORDER IN M                                       |
    !===================================================================================|
!
    !------------------------------------------------------------------------------
         Implicit None
         Integer, Intent (In) :: SZE
         Integer, Intent (Inout) :: N (SZE)
         Integer, Intent (Out) :: M (SZE)
    !------------------------------------------------------------------------------
         Integer I, LAST, J (1)
    !Real (SP) VALM !
         Intrinsic MINLOC
    !------------------------------------------------------------------------------
!
         M = 0
         LAST = 0
         Do I = 1, SZE
            J = MINLOC (N, MASK=N > LAST)
            M (I) = J (1)
            LAST = N (J(1))
         End Do
         Return
      End Subroutine SORT
  !===================================================================================|
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      Integer Function GETLOC (GLOC, GMAP, N)
    !===================================================================================|
    !    DETERMINE LOCAL IDENTITY OF ELEMENT/NODE I USING GMAP                           |
    !===================================================================================|
!
    !------------------------------------------------------------------------------
         Implicit None
         Integer, Intent (In) :: GLOC, N
         Integer, Intent (In) :: GMAP (N)
         Integer I, FOUND
    !------------------------------------------------------------------------------
         FOUND = 0
         Do I = 1, N
            If (GMAP(I) == GLOC) FOUND = I
         End Do
!
         GETLOC = FOUND
         Return
      End Function GETLOC
!
  !==============================================================================|
  !     READ IN REAL GLOBAL ARRAY AND DECOMPOSE TO LOCAL                         |
  !==============================================================================|
!
      Subroutine PREAD (IUNIT, VAR, I1, I2, IT, GI2, K2, GID, IRS, &
     & VNAME)
!
    !==============================================================================|
         Use MOD_UTILS, Only: PERROR

         Implicit None
         Integer, Intent (In) :: IUNIT, I1, I2, IT, GI2, K2, IRS
         Integer, Intent (In) :: GID (IT)
         Character (Len=*), Intent (In) :: VNAME
         Real (SP), Intent (Out), Dimension (I1:I2, K2) :: VAR
         Real (SP), Dimension (0:GI2, K2) :: TEMP
         Integer :: I, J
         Integer IERR
    !==============================================================================|
!
    !
    !--READ IN GLOBAL ARRAY
    !
         Read (IUNIT, IoStat=IERR) ((TEMP(I, J), J=1, K2), I=IRS, GI2)
    !
    !--ERROR CHECK
    !
         If (IERR /= 0) Call PERROR (IREPORT, "ERROR READING IN HOT START: ", &
        & VNAME)
!
    !
    !--DECOMPOSE
    !
         If (I1 == 0) VAR (0, :) = TEMP (0, :)
         Do I = 1, IT
            Do J = 1, K2
               VAR (I, J) = TEMP (GID(I), J)
            End Do
         End Do
         Return
      End Subroutine PREAD
  !==============================================================================|
!
  !==============================================================================|
  !     GATHER LOCAL ARRAY TO GLOBAL AND WRITE TO FILE                           |
  !==============================================================================|
!
      Subroutine PWRITE (IUNIT, MYID, NPROCS, VAR, I1, I2, IT, GI2, K2, &
     & MP, IRS)!, VNAME) !
!
    !==============================================================================|
         Implicit None
         Integer, Intent (In) :: IUNIT, I1, I2, IT, GI2, K2, IRS, MYID, &
        & NPROCS
         Type (GMAP), Intent (In) :: MP (NPROCS)
    !Character (Len=*), Intent (In) :: VNAME  !
         Real (SP), Intent (In), Dimension (I1:I2, K2) :: VAR
         Real (SP), Dimension (GI2, K2) :: TEMP
         Real (SP), Dimension (0:GI2, K2) :: TEMP2
         Integer :: I, J
    !==============================================================================|
!
    !
    !  GATHER
    !
         Call GATHER (I1, I2, IT, GI2, K2, MYID, NPROCS, MP, VAR, TEMP)
!
    !
    !  WRITE TO FILE
    !
         If (MYID == 1) Then
            TEMP2 (0, :) = 0.0_SP
            TEMP2 (1:GI2, 1:K2) = TEMP (1:GI2, 1:K2)
            Write (IUNIT) ((TEMP2(I, J), J=1, K2), I=IRS, GI2)
         End If
!
         Return
      End Subroutine PWRITE
  !==============================================================================|
  !==============================================================================|
!
      Subroutine PPRINT (IUNIT, LB1, UB1, UB2, VARP, VART, ILOC, K1, &
     & K2, REF, IRPT)
!
    !==============================================================================|
    !  WRITE OUT VARIABLE INFORMATION TO LOCAL FILES                               |
    !									       |
    !  USAGE EXAMPLES                                                              |
    !									       |
    !  write u velocity at surface in triangle 256 to file fort.306 with iteration |
    !  I1 = LBOUND(U,1) ; I2 = UBOUND(U,1)                                         |
    !  CALL PPRINT(306,I1,I2,KB,U,"element",256,1,1,FLOAT(IINT))                   |
    !									       |
    !  I1 = LBOUND(EL,1) ; I2 = UBOUND(EL,1)                                       |
    !  write surface elevation at node 233 to file fort.409 with time in hours     |
    !  CALL PPRINT(406,I1,I2,1,EL,"node",233,1,1,THOUR)                            |
    !									       |
    !  I1 = LBOUND(T1,1) ; I2 = UBOUND(T1,1)                                       |
    !  write vertical distribution of salinity at node 422 to file fort.433        |
    !  CALL PPRINT(433,I1,I2,KB,T1,"node",422,1,KBM1,THOUR)                        |
    !									       |
    !  ARGUMENT LIST                                                               |
    !       PPRINT(IUNIT,LB1,UB1,UB2,VARP,VART,ILOC,K1,K2,REF)                     |
    !    1.) IUNIT - UNIT NUMBER FOR OUTPUT FILE (MUST BE >= 300)                  |
    !    2.) LB1   - LOWER BOUND OF 1ST ARGUMENT OF ARRAY TO PRINT (USUAlLY 0)     |
    !    3.) LB2   - UPPER BOUND OF 1ST ARRAY DIMENSION (USUALLY NT OR MT)         |
    !        NOTE: LB1/LB2 CAN BE DETERMINE AUTOMATICALLY WITH LBOUND/UBOUND       |
    !    4.) UB2   - UPPER BOUND OF SECOND ARRAY DIMENSION                         |
    !        UB2   = 1 FOR SURFACE ARRAYS LIKE EL,UA                               |
    !        UB2   = KB FOR 3D ARRAYS LIKE U/V                                     |
    !    5.) VARP  = VARIABLE TO PRINT (ARRAY NAME = U,V,WW,EL,T1,RHO1, etc)       |
    !    6.) VART  = VARIABLE LOCATION ("element" or "node")                       |
    !    7.) ILOC  = INDEX OF ELEMENT/NODE TO PRINT                                |
    !    8.) K1    = LOWER RANGE OF SIGMA LEVEL TO PRINT                           |
    !    9.) K2    = UPPER RANGE OF SIGMA LEVEL TO PRINT                           |
    !        K1 = 1,K2 = 1 FOR SURFACE VALUES ONLY                                 |
    !        K1 = 1,K2 = KBM1 FOR ALL LEVELS                                       |
    !   10.) REF   = REFERENCE VALUE FOR DATA (MUST BE FLOAT)                      |
    !        REF = THOURS FOR CALCULATION TIME IN HOURS                            |
    !        REF = FLOAT(IINT) FOR ITERATION NUMBER                                |
    !   11.) IRPT = UNIT TO WRITE ERRORS TO (USE IREPORT)                          |
    !==============================================================================|
!
         Use MOD_LIMS, Only: NPROCS, MYID

         Implicit None
         Integer, Intent (In) :: IUNIT, LB1, UB1, UB2, ILOC, K1, K2, &
        & IRPT
         Real (SP), Intent (In) :: VARP (LB1:UB1, UB2)
         Real (SP), Intent (In) :: REF
         Character (Len=*), Intent (In) :: VART
!
         Integer :: I, J, K, PROCMAX, II, IBND
         Logical :: PRINT_PROC
!
    !==============================================================================|
!
    !------------------------------------------------------------------------------|
    !  Process Iunit for Errors                                                    |
    !------------------------------------------------------------------------------|
         If (IUNIT < 300) Then
            Write (IRPT,*) 'ERROR IN PPRINT'
            Write (IRPT,*) 'FILE UNIT < 300 RESERVED FOR FVCOM I/O'
            Write (IRPT,*) 'PLEASE INCREASE IUNIT TO 300+'
            Call PSTOP
         End If
!
    !------------------------------------------------------------------------------|
    !  Process Vartype for Errors                                                  |
    !------------------------------------------------------------------------------|
         If (VART /= "element" .And. VART /= "node") Then
            Write (IRPT,*) 'VART IN PPRINT NOT CORRECT', TRIM (VART)
            Write (IRPT,*) 'SHOULD BE "element" or "node"'
            Call PSTOP
         End If
!
    !------------------------------------------------------------------------------|
    !  Single Processor Case                                                       |
    !------------------------------------------------------------------------------|
         If (NPROCS == 1) Then
            Write (IUNIT,*) REF, (VARP(ILOC, K), K=K1, K2)
            Return
         End If
!
    !------------------------------------------------------------------------------|
    !  Multi Processor Case with Element Based Variable (u,v,ww, etc)              |
    !      Transform to Local Element ID with "ELID"                               |
    !------------------------------------------------------------------------------|
!
#  if defined (MULTIPROCESSOR)
         If (NPROCS /= 1 .And. VART == "element") Then
            If (ELID(ILOC) /= 0) Then
               Write (IUNIT,*) REF, (VARP(ELID(ILOC), K), K=K1, K2)
            End If
            Return
         End If
#  endif
!
    !------------------------------------------------------------------------------|
    !  Multi Processor Case with Node Based Variable (s1,t1,rho1,e1, etc)          |
    !      Transform to Local Node ID with "NLID"                                  |
    !      If Node is Interprocessor Boundary Node, Choose Processor with Highest  |
    !      ID Number to Write Values to File                                       |
    !------------------------------------------------------------------------------|
!
#  if defined (MULTIPROCESSOR)
         If (NPROCS /= 1 .And. VART == "node" .And. NLID(ILOC) > 0) &
        & Then
!
            If (NDE_ID(NLID(ILOC)) == 1) Then !!BOUNDARY NODE
               PRINT_PROC = .True.
!
               Do II = 1, NBN
                  If (BN_LST(II) == ILOC) IBND = II
               End Do
!
               PROCMAX = 10000
               Do J = 1, NPROCS
                  If (BN_NEY(IBND, J) == 1) Then
                     If (J < PROCMAX) PROCMAX = J
                  End If
               End Do
!
               If (PROCMAX /= MYID) PRINT_PROC = .False. !!NOT RESPONSIBLE FOR OUTPUT
            End If
!
            If (PRINT_PROC) Then
               Write (IUNIT,*) REF, (VARP(NLID(ILOC), K), K=K1, K2)
            End If
!
         End If
         Return
#  endif
!
         Return
      End Subroutine PPRINT
  !===================================================================================|
!
      Subroutine LOCK (MYID)
         Integer, Intent (In) :: MYID
         If (MYID >= 1) Then
       !        IF(ALLOCATED(FILELOCK))FILELOCK(MYID)=1  !set to one to lock
            MYFILELOCK = 1
         End If
         Return
      End Subroutine LOCK
!
      Subroutine UNLOCK (MYID)
         Integer, Intent (In) :: MYID
         If (MYID >= 1) Then
       !        IF(ALLOCATED(FILELOCK))FILELOCK(MYID)=0   !set to zero to unlock
            MYFILELOCK = 0
         End If
         Return
      End Subroutine UNLOCK
!
      Subroutine ROOT_OBTAINLOCK (MYID, ROOTID, NPROCS)
    !Passed variables
         Integer, Intent (In) :: MYID, ROOTID
         Integer, Intent (In) :: NPROCS
    !Local variables:
         Integer :: IERR, I
    !Basically to find out whether any process is locking a file
         Write (*,*) 'MYID is gathering :', MYID
         Call MPI_GATHER (MYFILELOCK, 1, MPI_INTEGER, FILELOCK, 1, &
        & MPI_INTEGER, ROOTID-1, MPI_COMM_WORLD, IERR)
         If (MYID .Eq. ROOTID) Then
            Do I = 1, NPROCS
               Write (*,*) 'MYID is:', MYID, 'FILELOCK', '(', I, ')=', &
              & FILELOCK (I)
            End Do
            If (MAXVAL(FILELOCK, 1) >= 1) Then
               FILELOCKED = .True. !file is being used by someone
            Else
               FILELOCKED = .False. !file is not being used by someone
            End If
         End If
         Return
      End Subroutine ROOT_OBTAINLOCK
!
      Subroutine PAR_DEALLOC
!
         If (ALLOCATED(NMAP)) DEALLOCATE (NMAP)
         If (ALLOCATED(EMAP)) DEALLOCATE (EMAP)
         If (ALLOCATED(EC)) DEALLOCATE (EC)
         If (ALLOCATED(BNC)) DEALLOCATE (BNC)
         If (ALLOCATED(NC)) DEALLOCATE (NC)
         If (ALLOCATED(PNHN)) DEALLOCATE (PNHN)
         If (ALLOCATED(EGID)) DEALLOCATE (EGID)
         If (ALLOCATED(NGID)) DEALLOCATE (NGID)
         If (ALLOCATED(HE_LST)) DEALLOCATE (HE_LST)
         If (ALLOCATED(ELID)) DEALLOCATE (ELID)
         If (ALLOCATED(ELID_X)) DEALLOCATE (ELID_X)
         If (ALLOCATED(NDE_ID)) DEALLOCATE (NDE_ID)
         If (ALLOCATED(BN_LST)) DEALLOCATE (BN_LST)
         If (ALLOCATED(BN_MLT)) DEALLOCATE (BN_MLT)
         If (ALLOCATED(BN_LOC)) DEALLOCATE (BN_LOC)
         If (ALLOCATED(BN_NEY)) DEALLOCATE (BN_NEY)
         If (ALLOCATED(HN_LST)) DEALLOCATE (HN_LST)
         If (ALLOCATED(HN_OWN)) DEALLOCATE (HN_OWN)
         If (ALLOCATED(NLID)) DEALLOCATE (NLID)
         If (ALLOCATED(NLID_X)) DEALLOCATE (NLID_X)
         If (ALLOCATED(PNE)) DEALLOCATE (PNE)
         If (ALLOCATED(PNN)) DEALLOCATE (PNN)
         If (ALLOCATED(PNHE)) DEALLOCATE (PNHE)
         If (ALLOCATED(PNBN)) DEALLOCATE (PNBN)
         If (ALLOCATED(PMBM)) DEALLOCATE (PMBM)
!
         If (ALLOCATED(EL_PID)) DEALLOCATE (EL_PID)
         If (ALLOCATED(FILELOCK)) DEALLOCATE (FILELOCK)
!
!
      End Subroutine PAR_DEALLOC
!
# endif
End Module MOD_PAR
