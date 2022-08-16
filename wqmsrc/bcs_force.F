!bcs_force.F
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
!subroutine BCS_FORCE()
!
!==============================================================================|
!   Set Up the Following Boundary Conditions:                                  |
!     Freshwater River Water Quality Input                                     |
!==============================================================================|
!
Subroutine BCS_FORCE
  !
  !------------------------------------------------------------------------------|
      Use MOD_SIZES, Only: NCP
  !
      Use MOD_LIMS, Only: NUMQBC_GL, NUMQBC, NUMPNT_GL, NUMPNT, KBM1, &
     & MLOC
      Use MOD_PREC, Only: SP
      Use MOD_HYDROVARS, Only: DZ, D

      Use MOD_CONTROL, Only: SERIAL, MSR, PAR, CASENAME, &
     & RIV_FILENUMBER, INFLOW_TYPE, POINT_ST_TYPE, PNT_SOURCE
!
      Use MOD_BCS, Only: INODEQ, ICELLQ, RIV_GL2LOC, INOPNT, QBC_TM, &
     & PNT_TM, WVQDIST, DWDIS, PDWQDIS, PDQDIS, N_ICELLQ
!
      Use MOD_UTILS, Only: PERROR, FOPEN
      Use MOD_WQM, Only: JDAY, IJDAY
	  USE MOD_FILEINFO, ONLY: INRIV, INPT,IREPORT
# if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: ELID, NLID
      Use MOD_PREC, Only: MPI_F
# endif
  !
      Implicit None
# if defined (MULTIPROCESSOR)
      Include "mpif.h"
# endif
      Character (Len=80) :: COMT
      Real (SP) :: TTIME
      Real (SP) :: FTEMP1, FTEMP2, FTEMP3, RBUF1, RBUF2, RBUF3
      Real (SP), Allocatable :: RTEMP (:, :), RTEMP1 (:, :, :), RTEMP2 &
     & (:, :)
      Integer, Allocatable :: RIV_2_LOC_PNT (:)
      Integer, Allocatable :: TEMP (:), TEMP2 (:), TEMP3 (:)
      Integer I, J, K, NQTIME, NBFWTIME, ISBCN1, INMAX, IOS, NCNT, IGL, &
     & IERR, JN
      Integer PNQTIME
      Character (Len=13) :: TSTRING
      Character (Len=100) :: ISTR, ISTPT
      Integer ::  NJ
  !
  !------------------------------------------------------------------------------|
  !
  !
  !----------------------------REPORT--------------------------------------------!
      If (MSR) WRITE (*,*) '!'
      If (MSR) WRITE (*,*) '!           SETTING UP PRESCRIBED BOUNDARY &
     &CONDITIONS   '
      If (MSR) WRITE (*,*) '!'
  !
  !==============================================================================|
  !   Input River/Dam/Intake/Outfall Boundary Values                             |
  !==============================================================================|
  !
  !
  !-------Check Selected Combination for Validity--------------------------------!
  !
      IJDAY = AINT (JDAY)
      Write (RIV_FILENUMBER, '(I4.4)') IJDAY
      ISTR = "./inputs/" // trim (CASENAME) // trim (RIV_FILENUMBER)
  !
      Call FOPEN (INRIV, trim(ISTR)//'_riv.dat', "cfr")
  !
      Rewind (INRIV)
      Read (INRIV, '(A4,2X,A10)') INFLOW_TYPE, POINT_ST_TYPE
  !
      If (INFLOW_TYPE /= 'edge' .And. INFLOW_TYPE /= 'node') Then
         Call PERROR (IREPORT, "INFLOW TYPE NOT CORRECT", "SHOULD BE edge or &
        &node")
      End If
  !
      If (POINT_ST_TYPE /= 'calculated' .And. POINT_ST_TYPE /= 'specifi&
     &ed') Then
         Call PERROR (IREPORT, "POINT_ST TYPE NOT CORRECT", "SHOULD BE calcul&
        &ated or specified")
      End If
  !
  !
  !--Read in Number of Discharge Nodes/Edges-------------------------------------!
  !
      If (MSR) Then
         Read (INRIV,*) NUMQBC_GL
      End If
  !
#  if defined (MULTIPROCESSOR)
      If (PAR) Call MPI_BCAST (NUMQBC_GL, 1, MPI_INTEGER, 0, &
     & MPI_COMM_WORLD, IERR)
#  endif
  !
      NUMQBC = 0
      If (NUMQBC_GL > 0) Then
     !
         Allocate (TEMP(NUMQBC_GL), TEMP2(NUMQBC_GL), TEMP3(NUMQBC_GL))
         If (MSR) Then
            Do I = 1, NUMQBC_GL
               Read (INRIV,*) TEMP (I)
            End Do
         End If
     !
#    if defined (MULTIPROCESSOR)
         If (PAR) Call MPI_BCAST (TEMP, NUMQBC_GL, MPI_INTEGER, 0, &
        & MPI_COMM_WORLD, IERR)
#    endif
     !
     !
     !--Determine Global--> Local Mapping of Freshwater Discharge Nodes
     !
         If (SERIAL) Then
            NUMQBC = NUMQBC_GL
            If (INFLOW_TYPE == 'node') Then
               Allocate (INODEQ(NUMQBC))
               INODEQ = TEMP
            Else If (INFLOW_TYPE == 'edge') Then
               Allocate (ICELLQ(NUMQBC))
               ICELLQ = TEMP (1:NUMQBC)
               Allocate (N_ICELLQ(NUMQBC, 2))!Added this to record the node numbes of the edge
           !which should be calculated using set_bndry.F
               N_ICELLQ = 0
            End If
         End If
     !
#   if defined (MULTIPROCESSOR)
         If (PAR) Then
            NCNT = 0
            If (INFLOW_TYPE == 'node') Then
               Do I = 1, NUMQBC_GL
                  If (NLID(TEMP(I)) /= 0) Then
                     NCNT = NCNT + 1
                     TEMP2 (NCNT) = NLID (TEMP(I))
                     TEMP3 (NCNT) = I
                  End If
               End Do
               NUMQBC = NCNT
               Allocate (INODEQ(NUMQBC), RIV_GL2LOC(NUMQBC))
               INODEQ = TEMP2 (1:NCNT)
               RIV_GL2LOC = TEMP3 (1:NCNT)
            Else If (INFLOW_TYPE == 'edge') Then
               Do I = 1, NUMQBC_GL
                  If (ELID(TEMP(I)) /= 0) Then
                     NCNT = NCNT + 1
                     TEMP2 (NCNT) = ELID (TEMP(I))
                     TEMP3 (NCNT) = I
                  End If
               End Do
               NUMQBC = NCNT
               Allocate (ICELLQ(NUMQBC), RIV_GL2LOC(NUMQBC))
               Allocate (N_ICELLQ(NUMQBC, 2))!Added this to record the node numbes of the edge
           !which should be calculated using set_bndry.F
               N_ICELLQ = 0
           !
               ICELLQ = TEMP2 (1:NCNT)
               RIV_GL2LOC = TEMP3 (1:NCNT)
            End If
         End If
#   endif
     !
         Deallocate (TEMP, TEMP2, TEMP3)
     !
     !
     !----Read in Time Dependent DataSets (DQDIS,DSDIS,DTDIS)------------------------!
     !
         If (MSR) READ (INRIV,*) NQTIME
     !
#    if defined (MULTIPROCESSOR)
         If (PAR) Call MPI_BCAST (NQTIME, 1, MPI_INTEGER, 0, &
        & MPI_COMM_WORLD, IERR)
#    endif
     !
         QBC_TM%NTIMES = NQTIME
         QBC_TM%LABEL = "Freshwater Discharge"
         Allocate (QBC_TM%TIMES(NQTIME))
         Allocate (RTEMP1(NUMQBC_GL, NCP, NQTIME))
     !    Added MSR to print out
     !     IF(MSR)WRITE(*,*)'NQTIME', NQTIME
         If (MSR) Then
            Do I = 1, NQTIME
               Read (INRIV,*) TTIME
               QBC_TM%TIMES (I) = TTIME
           !This is where wq data read in for rivers
           !need to create software to generate input file
               Do NJ = 1, NCP
                  Read (INRIV,*) (RTEMP1(J, NJ, I), J=1, NUMQBC_GL)
               End Do
           !
            End Do
         End If
     !
#    if defined (MULTIPROCESSOR)
         If (PAR) Call MPI_BCAST (QBC_TM%TIMES, NQTIME, MPI_F, 0, &
        & MPI_COMM_WORLD, IERR)
         If (PAR) Call MPI_BCAST (RTEMP1(:, :, :), &
        & NUMQBC_GL*NQTIME*NCP, MPI_F, 0, MPI_COMM_WORLD, IERR)
#    endif
     !
     !
     !----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
     !
     !
         If (NUMQBC > 0) Then
            Allocate (DWDIS(NUMQBC, NCP, NQTIME))
        !
            If (SERIAL) Then
               DWDIS (1:NUMQBC_GL, :, :) = RTEMP1 (1:NUMQBC_GL, :, :)
            End If
        !
#     if defined (MULTIPROCESSOR)
            If (PAR) Then
               Do I = 1, NQTIME
                  Do NJ = 1, NCP
                     DWDIS (1:NUMQBC, NJ, I) = RTEMP1 &
                    & (RIV_GL2LOC(1:NUMQBC), NJ, I)
                  End Do
               End Do
            End If
#     endif
        !
         End If
     !
         Deallocate (RTEMP1)
     !
         Close (INRIV)
     !
      End If !! NUMQBC_GL > 0
  !
  !==============================================================================
  ! NON-POINT source Values                             |
  !==============================================================================
  !
  !-------Check Selected Combination for Validity--------------------------------
  !
      ISTPT = "./inputs/" // trim (CASENAME)
      Call FOPEN (INPT, trim(ISTPT)//'_pnt_wq.dat', "cfr")!
      Rewind (INPT)
  !
      Read (INPT, '(A5)') PNT_SOURCE
  !
      If (PNT_SOURCE /= 'point') Then
         Call PERROR (IREPORT, "POINT_SOURCE INPUT FILE IS  NOT CORRECT")
      End If
  !
  !--Read in Number of Discharge Nodes/Edges-------------------------------------
  !
      Read (INPT,*) NUMPNT_GL
  !
      NUMPNT = 0
      If (NUMPNT_GL > 0) Then
     !
         Allocate (TEMP(NUMPNT_GL), TEMP2(NUMPNT_GL), TEMP3(NUMPNT_GL))
         Do I = 1, NUMPNT_GL
            Read (INPT,*) TEMP (I)
         End Do
     !
     !--Determine Global--> Local Mapping of Freshwater Discharge Nodes
     !
         If (SERIAL) Then
            NUMPNT = NUMPNT_GL
            Allocate (INOPNT(NUMPNT))
            INOPNT (:) = TEMP (:)
         End If
     !
#   if defined (MULTIPROCESSOR)
         If (PAR) Then
            NCNT = 0
            Do I = 1, NUMPNT_GL
               If (NLID(TEMP(I)) /= 0) Then
                  NCNT = NCNT + 1
                  TEMP2 (NCNT) = NLID (TEMP(I))
                  TEMP3 (NCNT) = I
               End If
            End Do
            NUMPNT = NCNT
            Allocate (RIV_2_LOC_PNT(NUMPNT), INOPNT(NUMPNT))
            INOPNT (1:NCNT) = TEMP2 (1:NCNT)! The point locally
            RIV_2_LOC_PNT (1:NCNT) = TEMP3 (1:NCNT)! The point globally
         End If
#   endif
     !
         Deallocate (TEMP, TEMP2, TEMP3)
     !
     !
     !----Read in Time Dependent Data and Flow Distribution Rate!
     !
         Allocate (RTEMP(NUMPNT_GL, KBM1))
         Do I = 1, NUMPNT_GL
            Read (INPT,*) J, (RTEMP(I, K), K=1, KBM1)
         End Do
     !
         Allocate (WVQDIST(NUMPNT, KBM1))
         If (SERIAL) WVQDIST (1:NUMPNT, :) = RTEMP (1:NUMPNT_GL, :)
#   if defined (MULTIPROCESSOR)
         If (PAR) Then
            WVQDIST (1:NUMPNT, :) = RTEMP (RIV_2_LOC_PNT(1:NUMPNT), :)
         End If
#   endif
     !
     !    READ DATA
         Read (INPT,*) PNQTIME
     !
         PNT_TM%NTIMES = PNQTIME
         PNT_TM%LABEL = "Point source"
         Allocate (PNT_TM%TIMES(PNQTIME))
         Allocate (RTEMP1(NUMPNT_GL, NCP, PNQTIME))
         Allocate (RTEMP2(NUMPNT_GL, PNQTIME))
     !
         Do I = 1, PNQTIME
            Read (INPT,*) TTIME
            PNT_TM%TIMES (I) = TTIME
            Read (INPT,*) (RTEMP2(J, I), J=1, NUMPNT_GL)
            Do NJ = 1, NCP
               Read (INPT,*) (RTEMP1(J, NJ, I), J=1, NUMPNT_GL)
            End Do
         End Do
     !
     !----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
     !
         Allocate (PDQDIS(NUMPNT, PNQTIME))
         Allocate (PDWQDIS(NUMPNT, NCP, PNQTIME))
     !
         If (SERIAL) Then
            PDQDIS (1:NUMPNT, :) = RTEMP2 (1:NUMPNT_GL, :)
            PDWQDIS (1:NUMPNT, :, :) = RTEMP1 (1:NUMPNT_GL, :, :)
         End If
     !
#     if defined (MULTIPROCESSOR)
         If (PAR) Then
            PDQDIS (1:NUMPNT, :) = RTEMP2 (RIV_2_LOC_PNT(1:NUMPNT), :)
            PDWQDIS (1:NUMPNT, :, :) = RTEMP1 (RIV_2_LOC_PNT(1:NUMPNT), &
           & :, :)
            Deallocate (RIV_2_LOC_PNT)
         End If
#     endif

         Deallocate (RTEMP, RTEMP1, RTEMP2)
     !
         Close (INPT)
     !
      End If !! NUMQBC_GL > 0
  !
      Return
End Subroutine BCS_FORCE
!==============================================================================|
