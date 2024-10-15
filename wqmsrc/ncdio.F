!ncdio.F
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
!Subroutines :
!	Subroutine NCD_READ_GRID()
!	Subroutine NCD_READ_SHAPE()
!	Subroutine NCD_READ_OPEN()
!	Subroutine NCD_READ()
!
Subroutine NCD_READ_GRID (INFILE)
  !---------------------------------------------------------------------
  ! READ DIMENSIONS IN A NETCDF FILES
  !---------------------------------------------------------------------
      Use MOD_NCD, Only: NC_FID, GETDIM, GETSVAR, GETDVAR
      Use NETCDF
      Use MOD_LIMS, Only: MLOC, NLOC, NTLOC, KBM2, KB, KBM1, MTLOC
      Use MOD_SIZES, Only: NGL, MGL
      Use MOD_PREC, Only: SP
  !Took CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: MSR
  !
      Implicit None
  !----------------------------------------------------------------------------!
      Character (Len=1024), Intent (In) :: INFILE
  !----------------------------------------------------------------------------!
      Integer :: IERR
      Integer :: N_ELEMS, N_NODES, N_SIG_M1, N_SIG
  !Real (SP), Allocatable, Dimension (:, :) :: TEMP !cleanup
  !----------------------------------------------------------------------------!
  !
  !--Open NetCDF DATA FILE
      IERR = NF90_OPEN (TRIM(INFILE), NF90_NOWRITE, NC_FID)
      If (MSR) WRITE (*,*) 'opening netcdf file, INFILE, NC_FID, NGL ' &
     & // TRIM (INFILE), NC_FID, NGL
      If (IERR /= NF90_NOERR) Then
         Write (*,*) 'ERROR READING ', TRIM (INFILE)
         Write (*,*) TRIM (NF90_STRERROR(IERR))
         Call PSTOP
      End If
  !
  !--Get Model Dimensions
      N_ELEMS = GETDIM (NC_FID, LEN_TRIM('nele'), 'nele')
      N_NODES = GETDIM (NC_FID, LEN_TRIM('node'), 'node')
      N_SIG_M1 = GETDIM (NC_FID, LEN_TRIM('siglay'), 'siglay')
      N_SIG = GETDIM (NC_FID, LEN_TRIM('siglev'), 'siglev')
  !
      If (N_ELEMS /= NGL) Then
         If (MSR) Then
            Write (*,*) 'NGL is inconsistent with the # of ELEMENTS in &
           &NetCDF file'
            Write (*,*) 'NGL should be equal to', N_ELEMS, NGL
            Call PSTOP
         End If
      End If
      If (N_NODES /= MGL) Then
         If (MSR) Then
            Write (*,*) 'MGL IS inconsistent with the # of NODES in Net&
           &CDF file'
            Write (*,*) 'MGL should be equal to', N_NODES
            Call PSTOP
         End If
      End If
      If (N_SIG /= KB) Then
         If (MSR) Then
            Write (*,*) 'KB IS inconsistent with the # of SIGMA LAYERS &
           &in NetCDF file'
            Write (*,*) 'KB should be equal to', N_SIG
            Call PSTOP
         End If
      End If
  !
      MLOC = N_NODES
      NLOC = N_ELEMS
      MTLOC = MLOC
      NTLOC = NLOC
  !
      KB = N_SIG
      KBM1 = N_SIG_M1
      KBM2 = KB - 2
  !
  !--close file
      IERR = NF90_CLOSE (NC_FID)
  !
      Return
End Subroutine NCD_READ_GRID
!
!==============================================================================|
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!==============================================================================|
!
Subroutine NCD_READ_SHAPE (INFILE)
  !---------------------------------------------------------------------
  ! READ BATHYMETRY, SIGMA LEVELS AND GRID COEFFICIENTS IN A NETCDF FILES
  !---------------------------------------------------------------------
      Use MOD_NCD, Only: NC_FID, GETSVAR
      Use NETCDF
      Use MOD_LIMS, Only: KBM2
      Use MOD_HYDROVARS, Only: NVG, XG, YG, HG, Z, Z2DG, &
     & num_hyd_ints,SIGVAR
  !Took CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR
  !
      Use MOD_LIMS, Only: KB, KBM1
      Use MOD_SIZES, Only: NGL, MGL
      Use MOD_PREC, Only: SP
  !
      Implicit None
  !----------------------------------------------------------------------------!
      Character (Len=1024), Intent (In) :: INFILE
  !----------------------------------------------------------------------------!
      Integer :: IERR
      Integer :: K, I !cleanup
      Real (SP), Allocatable, Dimension (:, :) :: TEMP
  !
  !--OPEN NETCDF DATA FILE
      IERR = NF90_OPEN (TRIM(INFILE), NF90_NOWRITE, NC_FID)
      If (MSR) WRITE (*,*) 'opening netcdf file, INFILE, NC_FID' // &
     & TRIM (INFILE), NC_FID
      If (IERR /= NF90_NOERR) Then
         Write (*,*) 'ERROR READING ', TRIM (INFILE)
         Write (*,*) TRIM (NF90_STRERROR(IERR))
         Stop
      End If
  !
  !--Get Node Coordinates
  !Moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F
  !ALLOCATE(XG(0:MGL),YG(0:MGL)) ; XG = 0.0_SP ; YG = 0.0_SP
  !
      Allocate (TEMP(MGL, 1))
      Call GETSVAR (NC_FID, LEN_TRIM('x'), 'x', MGL, 1, TEMP)
      XG (1:MGL) = TEMP (1:MGL, 1)
      Deallocate (TEMP)
  !
      Allocate (TEMP(MGL, 1))
      Call GETSVAR (NC_FID, LEN_TRIM('y'), 'y', MGL, 1, TEMP)
      YG (1:MGL) = TEMP (1:MGL, 1)
      Deallocate (TEMP)
  !
  !--Get Node Numbering
  !Moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F
  !ALLOCATE(NVG(0:NGL,4)); NVG = 0
  !
      Allocate (TEMP(NGL, 3))
      Call GETSVAR (NC_FID, LEN_TRIM('nv'), 'nv', NGL, 3, TEMP)
      NVG (1:NGL, 1:3) = TEMP (1:NGL, 1:3)
      Deallocate (TEMP)
      NVG (:, 4) = NVG (:, 1)
  !
  !--Get Bathymetry
  !!Moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F
  !ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP
  !
      Allocate (TEMP(MGL, 1))
      Call GETSVAR (NC_FID, LEN_TRIM('h'), 'h', MGL, 1, TEMP)
      HG (1:MGL) = TEMP (1:MGL, 1)
      Deallocate (TEMP)
  !
  !Moved this to HYDRO_GEOM_ALLOC() in mod_hydrovars.F
  ! ALLOCATE(Z(KB))  ; Z   = 0.0_SP    !!SIGMA COORDINATE VALUE
  ! ALLOCATE(ZZ(KB)) ; ZZ  = 0.0_SP    !!INTRA LEVEL SIGMA VALUE
  ! ALLOCATE(DZ(KB)) ; DZ  = 0.0_SP    !!DELTA-SIGMA VALUE
  ! ALLOCATE(DZZ(KB)); DZZ = 0.0_SP    !!DELTA OF INTRA LEVEL SIGMA
  !
  !--Get Sigma levels
      Allocate (TEMP(KB, 1))
      Call GETSVAR (NC_FID, LEN_TRIM('siglev'), 'siglev', KB, 1, TEMP)
      Z (1:KB) = TEMP (1:KB, 1)
      Deallocate (TEMP)
      Do K = 1, KB
     !tykim00
     !         Z(K)=0.0-(K-1)*(1.0/(KB*1.0-1.0))
     !
     !
     !: this vertical sigma coordinate formulation is hardired
     !       we should have this information in the netcdf file itself
     !       or we should get this **1.5 as a parameter from input control
     !
         Z (K) = - ((K-1)/FLOAT(KB-1)) ** 1.5
     !
      End Do

      call set_sigma

#if defined (NEWSIGMA)
	  IF(SIGVAR)THEN
		Allocate (TEMP(MGL, KB))
		Call GETSVAR (NC_FID, LEN_TRIM('siglev3dn'), 'siglev3dn', MGL, KB, TEMP)
		Z2DG (1:MGL,1:KB) = TEMP (1:MGL,1:KB)
		Deallocate (TEMP)
	  ENDIF
#endif

  !
  !--Close file
      IERR = NF90_CLOSE (NC_FID)
  !
      Return
End Subroutine NCD_READ_SHAPE
!
!==============================================================================|
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!==============================================================================|
!REMOVE ALL THE CODE FROM THIS ROUTINE AND CALL DUPLICATE CODE IN NCD_READ
Subroutine NCD_READ_OPEN (INFILE, UL, VL, WTSL, UARD_OBCNL, XFLUX_OBCL, &
& DTFAL, KHL, ELL, TL, SL, HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------
  !
      Use MOD_NCD, Only: NF90_OPEN, NF90_NOWRITE, NF90_STRERROR, &
     & NF90_NOERR, NC_FID
      Use MOD_SIZES, Only: NOBTY
      Use MOD_LIMS, Only: MTLOC, NTLOC, KBM1, KB
      Use MOD_PREC, Only: SP
  !
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR !	&           !!TRUE IF MULTIPROCESSOR RUN
  !,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
  !,CASETITLE  	&  			!!CASE TITLE
       !,HMAX       	&  			!!GLOBAL MAXIMUM DEPTH
  !,HMIN       	&  			!!GLOBAL MINIMUM DEPTH
       !,UMOL       	&  			!!VERTICAL DIFFUSION COEFFICIENT
  !,HORCON     	&  			!!HORIZONTAL DIFFUSION COEFFICIENT
       !,DTI        	&  			!!internal time step
  !,HORZMIX    	&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
       !,FILENUMBER	&			!!
  !,PREFF			&			!!
       !,INPDIR		&			!!
  !,GEOAREA		&			!!
       !,RIV_FILENUMBER	&			!!
  !,INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node)
       !,POINT_ST_TYPE 	&			!!(calculated/specified)
  !,PNT_SOURCE    	&			!!point_source
       !,DAY				&
  !,in_jday
      Implicit None
  !----------------------------------------------------------------------------!
  !
      Real (SP), Dimension (0:NTLOC, KB), Intent (Out) :: UL, VL !,WWL
      Real (SP), Dimension (0:MTLOC, KB), Intent (Out) :: KHL, WTSL
      Real (SP), Dimension (0:MTLOC, KBM1), Intent (Out) :: TL, SL
      Real (SP), Dimension (0:MTLOC), Intent (Out) :: ELL, DTFAL
      Real (SP), Dimension (0:NOBTY+1), Intent (Out) :: UARD_OBCNL
      Real (SP), Dimension (0:NOBTY, KBM1), Intent (Out) :: XFLUX_OBCL
      Integer, Intent (In) :: HO
      Character (Len=1024), Intent (In) :: INFILE
  !----------------------------------------------------------------------------!
      Integer :: IERR
  !  CHARACTER(LEN=300) TASK
  !
  !--Open NetCDF Datafile
      IERR = NF90_OPEN (TRIM(INFILE), NF90_NOWRITE, NC_FID)
  ! Added MSR to print out
      If (MSR) WRITE (*,*) 'opening netcdf file, INFILE, NC_FID' // &
     & TRIM (INFILE), NC_FID
  !
      If (IERR /= NF90_NOERR) Then
         Write (*,*) 'ERROR READING ', TRIM (INFILE)
         Write (*,*) TRIM (NF90_STRERROR(IERR))
         Stop
      End If
  !
      Call NCD_READ (INFILE, UL, VL, WTSL, UARD_OBCNL, XFLUX_OBCL, &
     & DTFAL, KHL, ELL, TL, SL, HO)
  !
      Return
End Subroutine NCD_READ_OPEN
!
!
!==============================================================================|
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!==============================================================================|
!
Subroutine NCD_READ (INFILE, UL, VL, WTSL, UARD_OBCNL, XFLUX_OBCL, &
& DTFAL, KHL, ELL, TL, SL, HO)
  !---------------------------------------------------------------------
  ! READ DATA FROM DAILY NETCDF FILES
  !---------------------------------------------------------------------
      Use MOD_PREC, Only: SP
      Use MOD_NCD, Only: NC_FID, GETDVAR
      Use NETCDF
  !
      Use MOD_WQM, Only: DTFAL_GL, ELL_GL, KHL_GL, SL_GL, TL_GL, &
     & UARD_OBCN_GL, UL_GL, VL_GL, WTSL_GL, XFLUX_OBC_GL
  !
      Use MOD_LIMS, Only: MTLOC, NTLOC, MLOC, NLOC, KB, KBM1
      Use MOD_SIZES, Only: MGL, NGL, NOBTY
  !
  !
  !
      Use MOD_HYDROVARS, Only: num_hyd_ints
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, MSR, PAR !	&           !!TRUE IF MULTIPROCESSOR RUN
  !,CASENAME  	&   		!!LETTER ACRONYM SPECIFYING CASE IDENTITY (MAX 80 CHARS)
  !,CASETITLE  	&  			!!CASE TITLE
       !,HMAX       	&  			!!GLOBAL MAXIMUM DEPTH
  !,HMIN       	&  			!!GLOBAL MINIMUM DEPTH
       !,UMOL       	&  			!!VERTICAL DIFFUSION COEFFICIENT
  !,HORCON     	&  			!!HORIZONTAL DIFFUSION COEFFICIENT
       !,DTI        	&  			!!internal time step
  !,HORZMIX    	&   		!!CONTROLS HORIZONTAL DIFFUSION COEF CALC (constant/closure)
       !,FILENUMBER	&			!!
  !,PREFF			&			!!
       !,INPDIR		&			!!
  !,GEOAREA		&			!!
       !,RIV_FILENUMBER	&			!!
  !,INFLOW_TYPE   	&			!!SPECIFIED RIVER INFLOW TYPE (edge/node)
       !,POINT_ST_TYPE 	&			!!(calculated/specified)
  !,PNT_SOURCE    	&			!!point_source
       !,DAY				&
  !,in_jday
  !
      Use MOD_BCMAP, Only: IOBCN, I_OBC_GL
  !
  !
# if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NLID, EGID, NGID, NHE, HE_LST, NBN, BN_LOC, &
     & BN_MLT, NHN, HN_LST, BNC
# endif
      Implicit None
  !
  !----------------------------------------------------------------------------!
  !
  !  REAL(CDF_PREC), DIMENSION(NGL,KBM1)   :: UL_GL,VL_GL                 !u and v
  !  REAL(CDF_PREC), DIMENSION(MGL,KB)     :: WTSL_GL,KHL_GL              !wts, kh
  !  REAL(CDF_PREC), DIMENSION(MGL,KBM1)   :: SL_GL,TL_GL                 !salinity,temp
  !  REAL(CDF_PREC), DIMENSION(MGL)        :: ELL_GL,DTFAL_GL             !zeta,dtfa
  !  REAL(CDF_PREC), DIMENSION(NOBTY+1)    :: UARD_OBCN_GL                !uard_obcn
  !  REAL(CDF_PREC), DIMENSION(NOBTY,KBM1) :: XFLUX_OBC_GL                !xflux_obc
  !
  !
      Real (SP), Dimension (0:NTLOC, KB), Intent (Out) :: UL, VL !,WWL
      Real (SP), Dimension (0:MTLOC, KB), Intent (Out) :: KHL, WTSL
      Real (SP), Dimension (0:MTLOC, KBM1), Intent (Out) :: TL, SL
      Real (SP), Dimension (0:MTLOC), Intent (Out) :: ELL, DTFAL
      Real (SP), Dimension (0:NOBTY+1), Intent (Out) :: UARD_OBCNL
      Real (SP), Dimension (0:NOBTY, KBM1), Intent (Out) :: XFLUX_OBCL
      Integer, Intent (In) :: HO
      Character (Len=1024), Intent (In) :: INFILE
  !----------------------------------------------------------------------------!
      Integer :: IERR
      Integer :: HT
      Integer :: I, K, J
  !
  !--Adjustement to read in Netcdf file
  !  HT=HO+1
      HT = HO
  !
  !---------------------------------------------------------------------
  ! Read Data from file INFILE at time level ht
  !---------------------------------------------------------------------
  !
  ! REORDER DATA READS TO MATCH ORDER IN FILE - FASTER
  !--U velocity
      Call GETDVAR (NC_FID, LEN_TRIM('u'), 'u', NGL, KBM1, UL_GL, HT)
  !
  !--V velocity
      Call GETDVAR (NC_FID, LEN_TRIM('v'), 'v', NGL, KBM1, VL_GL, HT)
  !
  !
  !--WTS velocity
      Call GETDVAR (NC_FID, LEN_TRIM('wts'), 'wts', MGL, KB, WTSL_GL, &
     & HT)
  !
  !--UARD_OBCN
      Call GETDVAR (NC_FID, LEN_TRIM('uard_obcn'), 'uard_obcn', NOBTY, &
     & 1, UARD_OBCN_GL, HT)
  !
  !--XFLUX_OBCN
      Call GETDVAR (NC_FID, LEN_TRIM('xflux_obc'), 'xflux_obc', NOBTY, &
     & KBM1, XFLUX_OBC_GL, HT)
  !
  !--free surface elevation
      Call GETDVAR (NC_FID, LEN_TRIM('dtfa'), 'dtfa', MGL, 1, DTFAL_GL, &
     & HT)
  !
  ! Fix dimensions of KH (KB, not KBM1)- here and above and below
  !--KH
      Call GETDVAR (NC_FID, LEN_TRIM('kh'), 'kh', MGL, KB, KHL_GL, HT)
  !
  !--free surface elevation
      Call GETDVAR (NC_FID, LEN_TRIM('zeta'), 'zeta', MGL, 1, ELL_GL, &
     & HT)
  !
  !--salinity
      Call GETDVAR (NC_FID, LEN_TRIM('salinity'), 'salinity', MGL, &
     & KBM1, SL_GL, HT)
  !
  !--temperature
      Call GETDVAR (NC_FID, LEN_TRIM('temp'), 'temp', MGL, KBM1, TL_GL, &
     & HT)
  !
  !  MOVED ALL COPIES TOGETHER
      If (SERIAL) Then
         ELL (1:MGL) = ELL_GL (1:MGL)
         SL (1:MGL, 1:KBM1) = SL_GL (1:MGL, 1:KBM1)
         TL (1:MGL, 1:KBM1) = TL_GL (1:MGL, 1:KBM1)
         UL (1:NGL, 1:KBM1) = UL_GL (1:NGL, 1:KBM1)
         VL (1:NGL, 1:KBM1) = VL_GL (1:NGL, 1:KBM1)
         WTSL (1:MGL, 1:KB) = WTSL_GL (1:MGL, 1:KB)
         UARD_OBCNL (1:NOBTY) = UARD_OBCN_GL (1:NOBTY)
         XFLUX_OBCL (1:NOBTY, 1:KBM1) = XFLUX_OBC_GL (1:NOBTY, 1:KBM1)
         DTFAL (1:MGL) = DTFAL_GL (1:MGL)
         KHL (1:MGL, 1:KB) = KHL_GL (1:MGL, 1:KB)
      End If
  !
# if defined (MULTIPROCESSOR)
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
            Call PSTOP
         End If
         XFLUX_OBCL = - 99999
         K = 0
         Do I = 1, NOBTY
            J = NLID (I_OBC_GL(I))
            If (J .Ne. 0) Then
               K = K + 1
               XFLUX_OBCL (K, :) = XFLUX_OBC_GL (I, :)
            End If
         End Do
         If (K .Ne. IOBCN) Then
            Write (*,*) "WRONG NUMBER OF BOUNDARIES ON A CPU"
            Call PSTOP
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
      End If
# endif
  !
      If (HT .Eq. num_hyd_ints) Then
         IERR = NF90_CLOSE (NC_FID)
      End If
  !
      Return
End Subroutine NCD_READ
