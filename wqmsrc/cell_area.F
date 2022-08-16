!cell_area.F
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
!==============================================================================|
!    This subroutine is used to calculate the area of individual               !
!    triangle based on the three vertex coordinates and also calculate         !
!    the sigma-surface area of individual control volume consisted of          !
!    triangles with a common node point                                        !
!                                          !
! calculates: art(ntri)   = area of element (triangle)                         !
! calculates: art1(nnode) = area of interior cv (for node value integration)   !
! calculates: art2(nnode) = sum area of all cells around node                  !
!==============================================================================|
!
!subroutine CELL_AREA()
!
Subroutine CELL_AREA
  !
  !==============================================================================!
      Use MOD_LIMS, Only: NTLOC, MTLOC, NLOC, MLOC, MYID, NPROCS
      Use MOD_PREC, Only: SP
  !
      Use MOD_TGE, Only: ISONB, NTVE, NBVE, NBVT, NV
  !
      Use MOD_HYDROVARS, Only: ZERO, XC, YC, VX, VY, ART, ART1, ART2
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: PAR
# if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NBN, BN_LOC, BN_MLT, BNC, NODE_MATCH
# endif
      Implicit None
      Real (SP), Allocatable :: XX (:), YY (:)
      Real (SP) :: ARTMAX, ARTTOT, ARTMIN
      Integer :: I, J, II, J1, J2, MAX_NBRE
  !==============================================================================!
  !
  !---------------INITIALIZE ARRAYS----------------------------------------------!

      MAX_NBRE = MAXVAL (NTVE) + 1
      Allocate (XX(2*MAX_NBRE+1), YY(2*MAX_NBRE+1))
      XX = 0.0_SP
      YY = 0.0_SP
  !
  !---------------COMPUTE AREA OF TRIANGLES USING CROSS PRODUCT------------------!
  !
      Do I = 1, NTLOC !surface area of all elements
         ART (I) = (VX(NV(I, 2))-VX(NV(I, 1))) * (VY(NV(I, 3))-VY(NV(I, &
        & 1))) - (VX(NV(I, 3))-VX(NV(I, 1))) * (VY(NV(I, 2))-VY(NV(I, &
        & 1)))
      End Do
      ART = Abs (.5_SP*ART)
  !
  !---------------COMPUTE MESH STATISTICS----------------------------------------!
  !
      ARTMIN = MINVAL (ART(1:NLOC))
      ARTMAX = MAXVAL (ART(1:NLOC))
      ARTTOT = SUM (ART(1:NLOC))
  !
  !-------COMPUTE CONTROL VOLUME ART1: CV FOR FLUXES OF NODAL BASED VALUES-------!
  !
      Do I = 1, MLOC
         If (ISONB(I) == 0) Then
            Do J = 1, NTVE (I)
               II = NBVE (I, J)
               J1 = NBVT (I, J)
               J2 = J1 + 1 - Int ((J1+1)/4) * 3
               XX (2*J-1) = (VX(NV(II, J1))+VX(NV(II, J2))) * 0.5_SP - &
              & VX (I)
               YY (2*J-1) = (VY(NV(II, J1))+VY(NV(II, J2))) * 0.5_SP - &
              & VY (I)
               XX (2*J) = XC (II) - VX (I)
               YY (2*J) = YC (II) - VY (I)
            End Do
            XX (2*NTVE(I)+1) = XX (1)
            YY (2*NTVE(I)+1) = YY (1)
        !
            Do J = 1, 2 * NTVE (I)
               ART1 (I) = ART1 (I) + 0.5_SP * &
              & (XX(J+1)*YY(J)-XX(J)*YY(J+1))
            End Do
            ART1 (I) = Abs (ART1(I))
         Else
            Do J = 1, NTVE (I)
               II = NBVE (I, J)
               J1 = NBVT (I, J)
               J2 = J1 + 1 - Int ((J1+1)/4) * 3
               XX (2*J-1) = (VX(NV(II, J1))+VX(NV(II, J2))) * 0.5_SP - &
              & VX (I)
               YY (2*J-1) = (VY(NV(II, J1))+VY(NV(II, J2))) * 0.5_SP - &
              & VY (I)
               XX (2*J) = XC (II) - VX (I)
               YY (2*J) = YC (II) - VY (I)
            End Do
            J = NTVE (I) + 1
            II = NBVE (I, J-1)
            J1 = NBVT (I, NTVE(I))
            J2 = J1 + 2 - Int ((J1+2)/4) * 3
        !
            XX (2*J-1) = (VX(NV(II, J1))+VX(NV(II, J2))) * 0.5_SP - VX &
           & (I)
            YY (2*J-1) = (VY(NV(II, J1))+VY(NV(II, J2))) * 0.5_SP - VY &
           & (I)
        !
            XX (2*J) = VX (I) - VX (I)
            YY (2*J) = VY (I) - VY (I)
        !
            XX (2*J+1) = XX (1)
            YY (2*J+1) = YY (1)
        !
            Do J = 1, 2 * NTVE (I) + 2
               ART1 (I) = ART1 (I) + 0.5_SP * &
              & (XX(J+1)*YY(J)-XX(J)*YY(J+1))
            End Do
            ART1 (I) = Abs (ART1(I))
         End If
      End Do
  !
  !---COMPUTE AREA OF CONTROL VOLUME ART2(I) = SUM(ALL TRIS SURROUNDING NODE I)--!
  !
      Do I = 1, MLOC
         ART2 (I) = SUM (ART(NBVE(I, 1:NTVE(I))))!summation of surrounding triangular elements area
     !
      End Do
  !
      ART (0) = ART (1)
      ART1 (0) = ART1 (1)
  !   IF(NTLOC > NLOC)ART(NLOC+1:NTLOC) = ART(NLOC)
      If (MTLOC > MLOC) ART2 (MLOC+1:MTLOC) = ART2 (MLOC)
      If (MTLOC > MLOC) ART1 (MLOC+1:MTLOC) = ART1 (MLOC)
      Deallocate (XX, YY)
  !
# if defined (MULTIPROCESSOR)
      If (PAR) Call NODE_MATCH (1, NBN, BN_MLT, BN_LOC, BNC, MTLOC, 1, &
     & MYID, NPROCS, ART1)
      If (PAR) Call NODE_MATCH (1, NBN, BN_MLT, BN_LOC, BNC, MTLOC, 1, &
     & MYID, NPROCS, ART2)
# endif
  !
  !
      Return
End Subroutine CELL_AREA
!==============================================================================|
