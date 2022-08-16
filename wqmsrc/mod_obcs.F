!mod_obcs.F
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
Module MOD_OBCS
  !
      Use MOD_PREC, Only: SP !
  !
      Use MOD_BCMAP, Only: IOBCN, IOBCN_GL, I_OBC_N !
  !
      Implicit None
      Save
  !
  !--Open Boundary Types, Lists, Pointers
  !Moved the following to bcmap.F
  !INTEGER               :: IOBCN_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY NODES   !move to bcmap.F
  !INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES    !move to bcmap.F
  !INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES		  !move to bcmap.F
  !INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST				  !move to bcmap.F
  !
      Integer, Allocatable :: NEXT_OBC (:)!!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
      Integer, Allocatable :: NEXT_OBC2 (:)!!INTERIOR NEIGHBOR OF NEXT_OBC
      Integer, Allocatable :: TYPE_OBC (:)!!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION) !never allocated
      Integer, Allocatable :: TYPE_OBC_GL (:)!!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION) !never allocated
      Integer :: IBCN (5)!!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
      Integer :: IBCN_GL (5)!!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
      Integer, Allocatable :: OBC_LST (:, :)!!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
      Integer, Allocatable :: NADJN_OBC (:)!!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
      Integer, Allocatable :: ADJN_OBC (:, :)!!ADJACENT OBNs of OBN
      Integer, Allocatable :: NADJC_OBC (:)!!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
      Integer, Allocatable :: ADJC_OBC (:, :)!!ADJACENT OPEN BOUNDARY CELLS
  !
  !--Open Boundary Metrics
      Integer, Allocatable :: NFLUXF_OBC (:)!!NUMBER OF FLUX SEGMENTS TO OBN   !never used
      Real (SP), Allocatable :: FLUXF_OBC (:, :)!!FLUX FRACTION ON EACH SIDE OF OBN	!never used
      Real (SP), Allocatable :: NXOBC (:)!!INWARD POINTING X-NORMAL OF OBN
      Real (SP), Allocatable :: NYOBC (:)!!INWARD POINTING Y-NORMAL OF OBN
      Real (SP), Allocatable :: DLTN_OBC (:)!!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY
  ! tykim
      Real (SP), Allocatable :: NUT_OBC (:, :)!!NUTRIENT AT OPEN BOUNDARY							!never used
      Real (SP), Allocatable :: NUT_OBC_GL (:, :)!!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY	!never used
  !
  !
Contains
  !
  !subroutines:
  !
  !	subroutine SEPARATE_OBC()
  !	subroutine SETUP_OBC()
  !	subroutine OBCS_DEALLOC()
  !functions:
  !
  !
  !==========================================================================|
  !
  !==========================================================================|
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !
  !==============================================================================|
      Subroutine SEPARATE_OBC !Function never used in this code
    !
    !------------------------------------------------------------------------------|
    ! Accumulate separately the amounts of nodes for 11 types of open boundaries   |
    !------------------------------------------------------------------------------|
    !
         Implicit None
    !
         Integer :: I, I1, I2, I3, I4, I5 !, II, J !!
    !
         IBCN = 0
         IBCN_GL = 0
    !
         Do I = 1, IOBCN_GL
            If (TYPE_OBC_GL(I) == 1 .Or. TYPE_OBC_GL(I) == 2) IBCN_GL &
           & (1) = IBCN_GL (1) + 1
            If (TYPE_OBC_GL(I) == 3 .Or. TYPE_OBC_GL(I) == 4) IBCN_GL &
           & (2) = IBCN_GL (2) + 1
            If (TYPE_OBC_GL(I) == 5 .Or. TYPE_OBC_GL(I) == 6) IBCN_GL &
           & (3) = IBCN_GL (3) + 1
            If (TYPE_OBC_GL(I) == 7 .Or. TYPE_OBC_GL(I) == 8) IBCN_GL &
           & (4) = IBCN_GL (4) + 1
            If (TYPE_OBC_GL(I) == 9 .Or. TYPE_OBC_GL(I) == 10) IBCN_GL &
           & (5) = IBCN_GL (5) + 1
         End Do
    !
         Do I = 1, IOBCN
            If (TYPE_OBC(I) == 1 .Or. TYPE_OBC(I) == 2) IBCN (1) = IBCN &
           & (1) + 1
            If (TYPE_OBC(I) == 3 .Or. TYPE_OBC(I) == 4) IBCN (2) = IBCN &
           & (2) + 1
            If (TYPE_OBC(I) == 5 .Or. TYPE_OBC(I) == 6) IBCN (3) = IBCN &
           & (3) + 1
            If (TYPE_OBC(I) == 7 .Or. TYPE_OBC(I) == 8) IBCN (4) = IBCN &
           & (4) + 1
            If (TYPE_OBC(I) == 9 .Or. TYPE_OBC(I) == 10) IBCN (5) = &
           & IBCN (5) + 1
         End Do
    !
         I1 = 0
         I2 = 0
         I3 = 0
         I4 = 0
         I5 = 0
    !
    !
         Allocate (OBC_LST(5, MAXVAL(IBCN)))
         OBC_LST = 0
    !
         Do I = 1, IOBCN
            If (TYPE_OBC(I) == 1 .Or. TYPE_OBC(I) == 2) Then
               I1 = I1 + 1
               OBC_LST (1, I1) = I
            Else If (TYPE_OBC(I) == 3 .Or. TYPE_OBC(I) == 4) Then
               I2 = I2 + 1
               OBC_LST (2, I2) = I
            Else If (TYPE_OBC(I) == 5 .Or. TYPE_OBC(I) == 6) Then
               I3 = I3 + 1
               OBC_LST (3, I3) = I
            Else If (TYPE_OBC(I) == 7 .Or. TYPE_OBC(I) == 8) Then
               I4 = I4 + 1
               OBC_LST (4, I4) = I
            Else If (TYPE_OBC(I) == 9 .Or. TYPE_OBC(I) == 10) Then
               I5 = I5 + 1
               OBC_LST (5, I5) = I
            End If
         End Do
    !
         Return
      End Subroutine SEPARATE_OBC
  !==============================================================================|
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !
  !==============================================================================|
      Subroutine SETUP_OBC
    !
         Use MOD_LIMS, Only: MYID
    !
         Use MOD_TGE, Only: ISONB, NTVE, NTSN, NBSN, NBVE, NV
    !
    !------------------------------------------------------------------------------!
         Use MOD_HYDROVARS, Only: XC, YC, VX, VY
         Implicit None
    !
         Real (SP) :: DXC, DYC, DXN, DYN, CROSS, DOTMAX, DOT, DX, DY, &
        & DXN_TMP, DYN_TMP !E1, E2, !
         Integer :: I, J, INODE, JNODE, I1, I2, IC, N1, N2, N3 !JJ, !
    !Logical :: DEBUG !
    !
         Real (SP), Allocatable :: NXOBC_TMP (:), NYOBC_TMP (:)
    !
    !------------------------------------------------------------------------------!
    !
    !--Determine Adjacent Open Boundary Points-------------------------------------!
    !
         Allocate (NADJN_OBC(IOBCN))
         NADJN_OBC = 0
         Allocate (ADJN_OBC(IOBCN, 2))
         ADJN_OBC = 0
    !
         Allocate (NEXT_OBC(IOBCN))
         NEXT_OBC = 0
         Allocate (NEXT_OBC2(IOBCN))
         NEXT_OBC2 = 0
    !
    !
         Do I = 1, IOBCN
            INODE = I_OBC_N (I)
            Do J = 1, NTSN (INODE) - 1
               JNODE = NBSN (INODE, J)
               If (ISONB(JNODE) == 2 .And. INODE /= JNODE) Then
                  NADJN_OBC (I) = NADJN_OBC (I) + 1
                  ADJN_OBC (I, NADJN_OBC(I)) = JNODE
               End If
            End Do
         End Do
    !
    !
         Do I = 1, IOBCN
            If (NADJN_OBC(I) == 0) Then
               Write (*,*) 'NO ADJACENT NODE FOUND FOR BOUNDARY NODE', &
              & I
               Write (*,*) 'IN PROCESSOR', MYID
               Call PSTOP
            End If
         End Do
    !
    !--Determine Adjacent Cells-(Nonlinear Only)-----------------------------------!
    !--Simultaneously Determine INWARD Pointing Normal NXOBC,NYOBC                 !
    !
         Allocate (NADJC_OBC(IOBCN))
         NADJC_OBC = 0
         Allocate (ADJC_OBC(IOBCN, 2))
         ADJC_OBC = 0
         Allocate (NXOBC(IOBCN))
         NXOBC = 0
         Allocate (NYOBC(IOBCN))
         NYOBC = 0
         Allocate (NXOBC_TMP(IOBCN))
         NXOBC_TMP = 0
         Allocate (NYOBC_TMP(IOBCN))
         NYOBC_TMP = 0
    !
         Do I = 1, IOBCN
            I1 = I_OBC_N (I)
       !
       !!Mark First Cell on Boundary Edge Adjacent to Node I
            I2 = ADJN_OBC (I, 1)
            Do J = 1, NTVE (I1)
               IC = NBVE (I1, J)
               N1 = NV (IC, 1)
               N2 = NV (IC, 2)
               N3 = NV (IC, 3)
               If (N1-I2 == 0 .Or. N2-I2 == 0 .Or. N3-I2 == 0) Then
                  DXN = VX (I2) - VX (I1)
                  DYN = VY (I2) - VY (I1)
                  DXC = XC (IC) - VX (I1)
                  DYC = YC (IC) - VY (I1)
             !
                  CROSS = SIGN (1.0_SP, DXC*DYN-DYC*DXN)
                  NXOBC_TMP (I) = CROSS * DYN / Sqrt (DXN**2+DYN**2)
                  NYOBC_TMP (I) = - CROSS * DXN / Sqrt (DXN**2+DYN**2)
                  NXOBC (I) = NXOBC_TMP (I)
                  NYOBC (I) = NYOBC_TMP (I)
                  NADJC_OBC (I) = NADJC_OBC (I) + 1
                  ADJC_OBC (I, NADJC_OBC(I)) = IC
               End If
            End Do
       !
            If (NADJN_OBC(I) > 1) Then
               I2 = ADJN_OBC (I, 2)
               Do J = 1, NTVE (I1)
                  IC = NBVE (I1, J)
                  N1 = NV (IC, 1)
                  N2 = NV (IC, 2)
                  N3 = NV (IC, 3)
                  If (N1-I2 == 0 .Or. N2-I2 == 0 .Or. N3-I2 == 0) Then
                     DXN = VX (I2) - VX (I1)
                     DYN = VY (I2) - VY (I1)
                     DXC = XC (IC) - VX (I1)
                     DYC = YC (IC) - VY (I1)
                !
                     CROSS = SIGN (1.0_SP, DXC*DYN-DYC*DXN)
                     NXOBC_TMP (I) = NXOBC_TMP (I) + CROSS * DYN / Sqrt &
                    & (DXN**2+DYN**2)
                     NYOBC_TMP (I) = NYOBC_TMP (I) - CROSS * DXN / Sqrt &
                    & (DXN**2+DYN**2)
                     NXOBC (I) = NXOBC_TMP (I) / Sqrt &
                    & (NXOBC_TMP(I)**2+NYOBC_TMP(I)**2)
                     NYOBC (I) = NYOBC_TMP (I) / Sqrt &
                    & (NXOBC_TMP(I)**2+NYOBC_TMP(I)**2)
                !
                     NADJC_OBC (I) = NADJC_OBC (I) + 1
                     ADJC_OBC (I, NADJC_OBC(I)) = IC
                  End If
               End Do
            End If
         End Do
    !
         Deallocate (NXOBC_TMP, NYOBC_TMP)
    !
    !--Determine 1st Layer Neighbor for Open Boundary Points-----------------------!
    !--Node Chosen is Node That is Connected to OBC Node and is Oriented           !
    !--Most Normal to the Boundary.  It is not Necessarily the Closest Node        !
    !--Determine also DLTN_OBC, the normal component of the distance between       !
    !--Next_obc and the open boundary node                                         !
    !
         Do I = 1, IOBCN
            DOTMAX = - 1.0
            INODE = I_OBC_N (I)
            Do J = 1, NTSN (INODE) - 1
               JNODE = NBSN (INODE, J)
               If (ISONB(JNODE) /= 2 .And. INODE /= JNODE) Then
                  DXN_TMP = VX (JNODE) - VX (INODE)
                  DYN_TMP = VY (JNODE) - VY (INODE)
             !
                  DXN = DXN_TMP / Sqrt (DXN_TMP**2+DYN_TMP**2)
                  DYN = DYN_TMP / Sqrt (DXN_TMP**2+DYN_TMP**2)
                  DOT = DXN * NXOBC (I) + DYN * NYOBC (I)
                  If (DOT > DOTMAX) Then
                     DOTMAX = DOT
                     NEXT_OBC (I) = JNODE
                  End If
               End If
            End Do
         End Do
    !
    !--Determine DLTN_OBC----------------------------------------------------------!
         Allocate (DLTN_OBC(IOBCN))
         Do I = 1, IOBCN
            I1 = I_OBC_N (I)
            I2 = NEXT_OBC (I)
       !
            DX = VX (I2) - VX (I1)
            DY = VY (I2) - VY (I1)
            DLTN_OBC (I) = Abs (DX*NXOBC(I)+DY*NYOBC(I))
         End Do
    !
         Return
      End Subroutine SETUP_OBC
  !
  !
      Subroutine OBCS_DEALLOC
    !
    !
         If (ALLOCATED(TYPE_OBC_GL)) DEALLOCATE (TYPE_OBC_GL)!never allocated
         If (ALLOCATED(NUT_OBC_GL)) DEALLOCATE (NUT_OBC_GL)!never allocated
         If (ALLOCATED(NUT_OBC)) DEALLOCATE (NUT_OBC)!never allocated
    !
         If (ALLOCATED(OBC_LST)) DEALLOCATE (OBC_LST)!mod_obcs
    !
         If (ALLOCATED(NADJN_OBC)) DEALLOCATE (NADJN_OBC)!mod_obcs
         If (ALLOCATED(ADJN_OBC)) DEALLOCATE (ADJN_OBC)!mod_obcs
         If (ALLOCATED(NEXT_OBC)) DEALLOCATE (NEXT_OBC)!mod_obcs
         If (ALLOCATED(NEXT_OBC2)) DEALLOCATE (NEXT_OBC2)!mod_obcs
         If (ALLOCATED(NADJC_OBC)) DEALLOCATE (NADJC_OBC)!mod_obcs
         If (ALLOCATED(ADJC_OBC)) DEALLOCATE (ADJC_OBC)!mod_obcs
         If (ALLOCATED(NXOBC)) DEALLOCATE (NXOBC)!mod_obcs
         If (ALLOCATED(NYOBC)) DEALLOCATE (NYOBC)!mod_obcs
         If (ALLOCATED(DLTN_OBC)) DEALLOCATE (DLTN_OBC)!mod_obcs
    !
    !
      End Subroutine OBCS_DEALLOC
  !
  !==============================================================================|
End Module MOD_OBCS
