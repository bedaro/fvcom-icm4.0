!pdomdec.F
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
!Subroutine PDOMDEC()
!
!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|
!
Subroutine PDOMDEC
  !
  !==============================================================================!
      Use MOD_LIMS, Only: NLOC, NTLOC, MLOC, MTLOC, KBM1,KBM2, NPROCS,KB
      Use MOD_PREC, Only: SP
      Use MOD_TGE, Only: NV
      Use MOD_HYDROVARS, Only: ZERO, NVG, XG, YG, HG, XCG, YCG, VXMIN, &
     & VYMIN, VXMAX, VYMAX, XC, YC, VX, VY, H1, H,Z2DG,Z2D,DZ2D, ZZ2D, DZZ2D,SIGVAR
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: SERIAL, PAR, HMAX, HMIN
  !
#  if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NLID, NLID_X, EGID, NGID, NHE, HE_LST, NHN, &
     & HN_LST
  !
#  endif
  !
      Use MOD_SIZES, Only: MGL, NGL
  !
      Implicit None
      Integer I, K, EGL !, J, IERR, I1, I2 !cleanup
  !
  !==============================================================================|
  !  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
  !  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
  !  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
  !==============================================================================|
  !Moved this to HYDRO_ALLOC()
  !ALLOCATE(NV(0:NTLOC,4));    NV = 0  !!NODE NUMBERING FOR ELEMENTS
      If (SERIAL) NV = NVG
  !
#  if defined (MULTIPROCESSOR)
      If (PAR) Then
         Do I = 1, NLOC
            EGL = EGID (I)
            NV (I, 1:4) = NLID (NVG(EGID(I), 1:4))
         End Do
         Do I = 1, NHE
            EGL = HE_LST (I)
            NV (I+NLOC, 1:4) = NLID_X (NVG(EGL, 1:4))
         End Do
      End If
#  endif
  !
  !==============================================================================|
  !   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
  !==============================================================================|
  !
  !--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!
  !
      VXMIN = MINVAL (XG(1:MGL))
      VXMAX = MAXVAL (XG(1:MGL))
      VYMIN = MINVAL (YG(1:MGL))
      VYMAX = MAXVAL (YG(1:MGL))
  !
  !--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!
  !
      XG = XG - VXMIN !DO NOT shift me!!! (removed the shift)
      YG = YG - VYMIN !DO NOT shift me!!! (removed the shift)
      XG (0) = 0.0_SP
      YG (0) = 0.0_SP
  !
  !--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!
  !Moved this to HYDRO_GEOM_ALLOC
  !ALLOCATE(XCG(0:NGL),YCG(0:NGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
  !
      Do I = 1, NGL
         XCG (I) = (XG(NVG(I, 1))+XG(NVG(I, 2))+XG(NVG(I, 3))) / 3.0_SP
         YCG (I) = (YG(NVG(I, 1))+YG(NVG(I, 2))+YG(NVG(I, 3))) / 3.0_SP
      End Do
  !
      XCG (0) = 0.0_SP
      YCG (0) = 0.0_SP
  !
  !
  !--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
  !Moved this to HYDO_ALLOC()
  !ALLOCATE(VX(0:MTLOC));        VX   = ZERO   !!X-COORD AT GRID POINT
  !ALLOCATE(VY(0:MTLOC));        VY   = ZERO   !!X-COORD AT GRID POINT
      If (SERIAL) Then
         VX = XG
         VY = YG
      End If
  !
#  if defined (MULTIPROCESSOR)
      If (PAR) Then
         Do I = 1, MLOC
            VX (I) = XG (NGID(I))
            VY (I) = YG (NGID(I))
         End Do
     !
         Do I = 1, NHN
            VX (I+MLOC) = XG (HN_LST(I))
            VY (I+MLOC) = YG (HN_LST(I))
         End Do
      End If
# endif
  !
  !==============================================================================|
  !   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
  !==============================================================================|
  !
  !--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
  !Moved this to mod_hydrovars.F HYDRO_ALLOC
  !ALLOCATE(H(0:MTLOC));       H = ZERO       !!BATHYMETRIC DEPTH
      If (SERIAL) H = HG
  !
#  if defined (MULTIPROCESSOR)
      If (PAR) Then
         Do I = 1, MLOC
            H (I) = HG (NGID(I))
         End Do
         Do I = 1, NHN
            H (I+MLOC) = HG (HN_LST(I))
         End Do
      End If
#  endif


  !
  !--------------CALCULATE EXTREMUMS---------------------------------------------!
  !
      HMAX = MAXVAL (Abs(HG(1:MGL)))
      HMIN = MINVAL (HG(1:MGL))
  !
  !==============================================================================|
  !   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
  !==============================================================================|
  !
  !Moved these to HYDRO_ALLOC
  !ALLOCATE(XC(0:NTLOC))       ;XC   = ZERO   !!X-COORD AT FACE CENTER
  !ALLOCATE(YC(0:NTLOC))       ;YC   = ZERO   !!X-COORD AT FACE CENTER
  !ALLOCATE(H1(0:NTLOC))       ;H1   = ZERO   !!BATHYMETRIC DEPTH
  !
  !
      Do I = 1, NTLOC
         XC (I) = (VX(NV(I, 1))+VX(NV(I, 2))+VX(NV(I, 3))) / 3.0_SP
         YC (I) = (VY(NV(I, 1))+VY(NV(I, 2))+VY(NV(I, 3))) / 3.0_SP
         H1 (I) = SUM (H(NV(I, 1:3))) / 3.0_SP
      End Do
  !
	  IF(SERIAL) THEN
			Z2D=Z2DG
			!DO K=1,KB
			!DO I=1,MTLOC
			!	write(*,*)'I,K,Z2D(I,K)=' ,I,K,Z2D(I,K)
			!	write(*,*)'I,K,Z2DG(I,K)=',I,K,Z2DG(I,K)
			!ENDDO
			!ENDDO
	  ENDIF

#     if defined (MULTIPROCESSOR)
	  IF(PAR)THEN
         DO I=1,MLOC
         DO K=1,KB
			Z2D(I,K)=Z2DG(NGID(I),K)
         ENDDO
         ENDDO
         DO I=1,NHN
         DO K=1,KB
            Z2D(I+MLOC,K)=Z2DG(HN_LST(I),K)
         ENDDO
         ENDDO
	  ENDIF
#     endif

	  !--Compute derivative and intra-sigma levels
	  DO I=1,MTLOC
			Do K = 1, KBM1
				ZZ2D (I,K) = 0.5_SP * (Z2D(I,K)+Z2D(I,K+1))
				DZ2D (I,K) = Z2D (I,K) - Z2D (I,K+1)
			End Do

			ZZ2D (I,KB) = 2.0_SP * ZZ2D (I,KBM1) - ZZ2D (I,KBM2)

			Do K = 1, KBM2
				DZZ2D (I,K) = ZZ2D (I,K) - ZZ2D (I,K+1)
			End Do
			DZZ2D (I,KB-1) = 0.0
			DZ2D (I,KB) = 0.0
	  ENDDO

      Return
End Subroutine PDOMDEC
!==============================================================================|
