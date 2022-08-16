!mod_wd.F
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
MODULE MOD_WD  !Wetting and Drying Module for FVCOM-ICM
!9/12/2016
!PNNL/MSL/BSRC, Seattle, WA

#if defined (WET_DRY)
   USE MOD_PREC, ONLY: SP
   USE MOD_HYDROVARS
   IMPLICIT NONE
   SAVE
!
!--Parameters for Wet/Dry Treatment
!
   LOGICAL  :: WET_DRY_ON

!--variables controlling porosities through wet/dry determination-------------------!

   INTEGER , ALLOCATABLE :: ISWETN(:)  !!NODE POROSITY AT NODES FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETC(:)  !!CELL POROSITY AT CELLS FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETNT(:) !!NODE POROSITY AT NODES FOR TIME N-1 INTERNAL
   INTEGER , ALLOCATABLE :: ISWETCT(:) !!CELL POROSITY AT CELLS FOR TIME N-1 INTERNAL
!   INTEGER , ALLOCATABLE :: ISWETCE(:) !!CELL POROSITY AT CELLS FOR TIME N-1 EXTERNAL
   REAL(SP), ALLOCATABLE :: FWET_N_N(:)!!NODE POROSITY AT NODES FOR TIME N
   REAL(SP), ALLOCATABLE :: FWET_C_C(:)!!CELL POROSITY AT CELLS FOR TIME N
   REAL (SP) :: MIN_DEPTH

   CONTAINS

	!------------------------------------------------------------------------------!
	! WD_SET_PARAM        :   READ WET/DRY PARAMETERS FROM INPUT                   !
	! WD_ALLOC		      :   ALLOCATE AND INITIALIZE WET/DRY ARRAYS               !
	! WD_SET_DATA         :   SET VALUES IN WET/DRY ARRAYS                         !
	! WD_WET_JUDGE        :   DETERMINE IF NODES/ELEMENTS ARE WET/DRY              !
	! WD_UPDATE           :   SWAP WET/DRY VARIABLES BETWEEN TIME LEVS             !
	! WD_WRITE_RESTART    :   DUMP WET/DRY FLAGS FOR RESTART                       !
	! WD_READ_RESTART     :   READ WET/DRY FLAGS FOR RESTART                       !
	! -----------------------------------------------------------------------------!


   SUBROUTINE WD_SET_PARAM

!----------------------------------------------------------------------------------|
!  READ IN PARAMETERS CONTROLLING WET/DRY TREATMENT                                |
!----------------------------------------------------------------------------------|


   USE MOD_INPUT, ONLY: SCAN_FILE
   USE MOD_FILEINFO, ONLY:IREPORT,INCTR
   USE MOD_CONTROL, ONLY: CASENAME, MSR,SERIAL
   IMPLICIT NONE
   INTEGER  :: ISCAN
   CHARACTER(LEN=120) :: FNAME


!------------------------------------------------------------------------------|
!   READ IN VARIABLES AND SET VALUES                                           |
!------------------------------------------------------------------------------|

   FNAME = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!   WET/DRY ACTIVE/INACTIVE FLAG
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"WET_DRY_ON",LVAL = WET_DRY_ON, IUNIT=INCTR)
   IF(ISCAN /= 0)THEN
     WRITE(IREPORT,*)'ERROR READING WET_DRY_ON: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IREPORT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!   WET/DRY minimum depth
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"MIN_DEPTH",FSCAL = MIN_DEPTH,IUNIT=INCTR)
   IF(ISCAN /= 0)THEN
     WRITE(IREPORT,*)'ERROR READING MIN_DEPTH: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IREPORT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABlES                                    !
!------------------------------------------------------------------------------|
   IF(MSR)THEN
     WRITE(IREPORT,*)''
     WRITE(IREPORT,*)'!        WET/DRY PARAMETERS                 '
     IF(WET_DRY_ON)THEN
       WRITE(IREPORT,*)'!  # WET/DRY TREATMENT   :  ACTIVE'
	   WRITE(IREPORT,*)'!  # WET/DRY MIN DEPTH   :  ', MIN_DEPTH
     ELSE
       WRITE(IREPORT,*)'!  # WET/DRY TREATMENT   :  NOT ACTIVE'
     END IF
   END IF

   RETURN
   END SUBROUTINE WD_SET_PARAM


!==============================================================================|
!==============================================================================|

   SUBROUTINE WD_SET_DATA

!------------------------------------------------------------------------------|
!  INITIALIZE ARRAYS USED FOR WET/DRY TREATMENT                                |
!------------------------------------------------------------------------------|



   USE MOD_LIMS, ONLY: NPROCS, MYID, MTLOC
   USE MOD_CONTROL, ONLY: PAR
   USE MOD_WQMINIT, ONLY: RES_ICON_IN
   USE MOD_TGE, ONLY: NV

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:NODE_MATCH,NBN,BN_MLT, BN_LOC, BNC
#  endif
   IMPLICIT NONE
   INTEGER :: I

!   IF(.NOT.RES_ICON_IN)THEN  !if cold start then do the very first wet_dry judging
   !IF(RESTART == 'cold_start') THEN

!-------- SET WET/DRY FLAGS AND MODIFY WATER SURFACE ELEVATION-----------------!

     CALL WD_WET_JUDGE

!-------- EXCHANGE MODIFIED FREE SURFACE ELEVATION ACROSS PROCESSOR BOUNDS-----!

#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MTLOC,1,MYID,NPROCS,EL)
#  endif

!-------- TRANSFER ELEVATION FIELD TO DEPTH AND OLD TIME LEVELS----------------!

	 !Get EL1 from EL
     !Do I = 1, NTLOC
	 !    EL1 (I) = EL (NV(I, 1)) + EL (NV(I, 2)) + EL (NV(I, 3))
     !    EL1 (I) = EL1 (I) / 3.0
     !End Do

	 !get D from ELF
     D  = EL + H
	 !Do I = 1, NTLOC
	 !    D1 (I) = D (NV(I, 1)) + D (NV(I, 2)) + D (NV(I, 3))
     !    D1 (I) = D1 (I) / 3.0
     !End Do

	 !set EL to ET
     ET = EL

	 !set EL1 to ET1
	 !ET1 = EL1

	 !set D to DT and DT1
     DT = D

     Do I = 1, NTLOC
            DT1 (I) = DT (NV(I, 1)) + DT (NV(I, 2)) + DT (NV(I, 3))
            DT1 (I) = DT1 (I) / 3.0
     End Do

	 !DT1=D1

     DTFA = D
!   END IF

   RETURN
   END SUBROUTINE WD_SET_DATA

!==============================================================================|

   SUBROUTINE WD_ALLOC

!------------------------------------------------------------------------------|
!  ALLOCATE AND INITIALIZE WET/DRY TREATMENT ARRAYS                            |
!------------------------------------------------------------------------------|


   USE MOD_LIMS, ONLY: MTLOC, NTLOC
   IMPLICIT NONE

!-----variables controlling porosities through wet/dry determination----------------!
   !IF(WET_DRY_ON)THEN
		ALLOCATE(ISWETN(0:MTLOC))        ; ISWETN     = 1
		ALLOCATE(ISWETC(0:NTLOC))        ; ISWETC     = 1
		ALLOCATE(ISWETNT(0:MTLOC))       ; ISWETNT    = 1
		ALLOCATE(ISWETCT(0:NTLOC))       ; ISWETCT    = 1
!		ALLOCATE(ISWETCE(0:NTLOC))       ; ISWETCE    = 1
		ALLOCATE(FWET_N_N(0:MTLOC))      ; FWET_N_N   = 1.0_SP
		ALLOCATE(FWET_C_C(0:NTLOC))      ; FWET_C_C   = 1.0_SP
   !ENDIF

   RETURN
   END SUBROUTINE WD_ALLOC


   SUBROUTINE WD_DEALLOC

!------------------------------------------------------------------------------|
!  DEALLOCATE all wet dry arrays
!------------------------------------------------------------------------------|

   IMPLICIT NONE

		IF(ALLOCATED(ISWETN))DEALLOCATE(ISWETN)
		IF(ALLOCATED(ISWETC))DEALLOCATE(ISWETC)
		IF(ALLOCATED(ISWETNT))DEALLOCATE(ISWETNT)
		IF(ALLOCATED(ISWETCT))DEALLOCATE(ISWETCT)
!		IF(ALLOCATED(ISWETCE))DEALLOCATE(ISWETCE)
		IF(ALLOCATED(FWET_N_N))DEALLOCATE(FWET_N_N)
		IF(ALLOCATED(FWET_C_C))DEALLOCATE(FWET_C_C)


   RETURN
   END SUBROUTINE WD_DEALLOC

!==============================================================================|

   SUBROUTINE WD_WET_JUDGE

!------------------------------------------------------------------------------|
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC, ONLY: SP
   USE MOD_LIMS, ONLY: NPROCS, MYID,MLOC, NLOC, MTLOC, NTLOC
   USE MOD_CONTROL, ONLY: PAR
   USE MOD_TGE, ONLY: NV, NBVE,NTVE

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:EXCHANGE,EC,NC
#  endif
   IMPLICIT NONE

   Real (SP), Parameter :: ONE_THIRD = 1.0 / 3.0

   REAL (SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6


!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWETN = 1
   DO I = 1, MLOC
     DTMP = H(I) + EL(I)
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETN(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWETC = 1
   DO I = 1, NLOC
     DTMP =  MAX(EL(NV(I,1)),EL(NV(I,2)),EL(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETC(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, MLOC
     IF(SUM(ISWETC(NBVE(I,1:NTVE(I)))) == 0)  ISWETN(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   EL = MAX(EL,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   !DO I = 1, NLOC
   !  EL1(I) = ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
   !END DO

!
!--Extend Element/Node Based Wet/Dry Flags to  Halo Domain----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_N_N = FLOAT(ISWETN)
     FWET_C_C = FLOAT(ISWETC)
     CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,MTLOC,1,MYID,NPROCS,FWET_N_N)
     ISWETN = INT(FWET_N_N+.5)
     ISWETC = INT(FWET_C_C+.5)
   END IF
#  endif

   RETURN
   END SUBROUTINE WD_WET_JUDGE

!==============================================================================|


   SUBROUTINE WD_UPDATE(INCASE)

!------------------------------------------------------------------------------|
!  SHIFT WET/DRY VARIABLES TO NEW TIME LEVELS                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC, ONLY: SP
   USE MOD_LIMS, ONLY: NPROCS, MYID, MLOC
   USE MOD_CONTROL, ONLY: PAR

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: INCASE
   INTEGER :: I


   SELECT CASE(INCASE)

!------------------------------------------------------------------------------!
   CASE(1)    !! SHIFT AT END OF EXTERNAL MODE
!------------------------------------------------------------------------------!
!		ISWETCE=ISWETC
!------------------------------------------------------------------------------!
   CASE(2)    !! UPDATE NODE WET/DRY AFTER DEPTH ADJUSTMENT
!------------------------------------------------------------------------------!
		DO I = 1,MLOC
			IF(DTFA(I)-MIN_DEPTH <= 1.0E-5_SP) THEN
				ISWETN(I) = 0
			END IF
		END DO
!------------------------------------------------------------------------------!
   CASE(3)    !! SHIFT VARIABLES AT END OF INTERNAL MODE
!------------------------------------------------------------------------------!
		ISWETCT=ISWETC
		ISWETNT=ISWETN

   END SELECT
   RETURN
   END SUBROUTINE WD_UPDATE

!===============================================================================|
!   DUMP WET/DRY FLAG DATA FOR RESTART                                          |
!===============================================================================|

   SUBROUTINE WD_WRITE_RESTART(FNAME)

!-------------------------------------------------------------------------------|

   USE MOD_PREC, ONLY: SP
   USE MOD_LIMS, ONLY: NPROCS, MYID, MLOC, NLOC,IINT
   USE MOD_CONTROL, ONLY: PAR, MSR, SERIAL

   USE MOD_SIZES, ONLY: MGL, NGL
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY: IGATHER, EMAP, NMAP
#  endif
   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I
   CHARACTER(LEN=120) :: FNAME
!==============================================================================|

   IF(MSR)THEN
     OPEN(1,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
     REWIND(1)
     WRITE(1,*) IINT
     WRITE(1,*) NGL,MGL
   END IF

   IF(SERIAL)THEN
     WRITE(1,*) (ISWETC(I), I=1,NLOC)
     WRITE(1,*) (ISWETN(I), I=1,MLOC)
   ELSE
   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
#  if defined (MULTIPROCESSOR)
   CALL IGATHER(LBOUND(ISWETC,1),UBOUND(ISWETC,1),NLOC,NGL,1,MYID,NPROCS,EMAP,ISWETC,NTEMP1)
   CALL IGATHER(LBOUND(ISWETN,1),UBOUND(ISWETN,1),MLOC,MGL,1,MYID,NPROCS,NMAP,ISWETN,NTEMP2)
   IF(MSR)THEN
     WRITE(1,*) (NTEMP1(I), I=1,NGL)
     WRITE(1,*) (NTEMP2(I), I=1,MGL)
   END IF
   DEALLOCATE(NTEMP1,NTEMP2)
#  endif
   END IF

   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE WD_WRITE_RESTART

!==============================================================================|
!   READ WET/DRY FLAG DATA FOR RESTART                                         |
!==============================================================================|

   SUBROUTINE WD_READ_RESTART(FNAME)

   USE MOD_FILEINFO, ONLY:IREPORT
   USE MOD_LIMS, ONLY: NPROCS, MYID, MLOC, NLOC, MTLOC, NTLOC,IINT

   USE MOD_CONTROL, ONLY: PAR, MSR, SERIAL
   USE MOD_SIZES, ONLY: MGL, NGL
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY:EXCHANGE,EC, NC,EGID,NGID
#  endif
   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I,IINT_TMP
   CHARACTER(LEN=120) :: FNAME
   REAL(SP) :: DTI_RATIO
   LOGICAL :: FEXIST

   DTI_RATIO=1.0_SP

!   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
!
!--Make Sure Wet/Dry Flag Data Exists------------------------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IREPORT,*)'WET/DRY RESTART FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IREPORT,*)'HALTING.....'
     CALL PSTOP
   END IF

!
!--Ensure File Header is Consistent with Main Flow Variable Restart File-------!
!
   OPEN(1,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
   REWIND(1)

!C. Lee, 12.12 2006, restart T.W.
   READ(1,*) IINT_TMP
   IINT_TMP=DTI_RATIO*IINT_TMP
   IINT=DTI_RATIO*IINT
!c. Lee end T.W.

   READ(1,*)
   IF(IINT_TMP /= IINT .AND. MSR)THEN
     WRITE(IREPORT,*)'IINT IN ',FNAME,' NOT EQUAL TO IINT'
     WRITE(IREPORT,*)'FROM MAIN RESTART FILE',IINT,IINT_TMP
     CALL PSTOP
   END IF

!
!--Read Variables--------------------------------------------------------------!
!
   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
   READ(1,*) (NTEMP1(I), I=1,NGL)
   READ(1,*) (NTEMP2(I), I=1,MGL)

!
!--Transfer Variables to Local Domains-----------------------------------------!
!
   IF(SERIAL)THEN
     ISWETC(1:NLOC) = NTEMP1(1:NLOC)
     ISWETN(1:MLOC) = NTEMP2(1:MLOC)
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,NLOC
       ISWETC(I) = NTEMP1(EGID(I))
     END DO
     DO I=1,MLOC
       ISWETN(I) = NTEMP2(NGID(I))
     END DO
   END IF
#  endif

   DEALLOCATE(NTEMP1,NTEMP2)
   CLOSE(1)
!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_C_C = ISWETC
     FWET_N_N = ISWETN
     CALL EXCHANGE(EC,NTLOC,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,MTLOC,1,MYID,NPROCS,FWET_N_N)
     ISWETN = INT(FWET_N_N+.5)
     ISWETC = INT(FWET_C_C+.5)
   END IF
#  endif

   ISWETNT = ISWETN
   ISWETCT = ISWETC
!   ISWETCE = ISWETC
!
!--Extend Element/Node Based Wet/Dry Flags to Previous Time Levels-------------!
!
   RETURN
   END SUBROUTINE WD_READ_RESTART

# endif
END MODULE MOD_WD
