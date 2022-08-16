!viscofh.F
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
!Subroutine VISCOF_H()
!
!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|
!
Subroutine VISCOF_H
  !
  !------------------------------------------------------------------------------|
      Use MOD_LIMS, Only: KBM1, MTLOC, KB, MLOC
  !
      Use MOD_TGE, Only: ISONB, NTVE, NBVT, NBVE, NV
  !
      Use MOD_HYDROVARS, Only: XC, YC, VX, VY, ART1, UU, VV, VISCOFH
      Use MOD_PREC, Only: SP
  !
      Implicit None
  ! REMOVE MLOC DIMENSION - NOT NEEDED FOR PUPX,PUPY,PVPX,PVPY & VISCOFF
      Real (SP) :: PUPX, PUPY, PVPX, PVPY
      Real (SP) :: VISCOFF
      Real (SP) :: X11, Y11, X22, Y22, X33, Y33, TMP1, TMP2
      Integer :: I, I1, IA, IB, J, J1, J2, K, JTMP
  !
  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
  !
      Do K = 1, KBM1
         Do I = 1, MLOC
            PUPX = 0.0_SP
            PUPY = 0.0_SP
            PVPX = 0.0_SP
            PVPY = 0.0_SP
        !
        ! FOLD IN SPECIAL CASE FOR J=1 AND J=NTVE(I)
        ! HAVE TO LEAVE ISONB(I) /= 0 BELOW
            Do J = 1, NTVE (I)
               I1 = NBVE (I, J)
               JTMP = NBVT (I, J)
               J1 = JTMP + 1 - (JTMP+1) / 4 * 3
               J2 = JTMP + 2 - (JTMP+2) / 4 * 3
               X11 = 0.5_SP * (VX(I)+VX(NV(I1, J1)))
               Y11 = 0.5_SP * (VY(I)+VY(NV(I1, J1)))
               X22 = XC (I1)
               Y22 = YC (I1)
               X33 = 0.5_SP * (VX(I)+VX(NV(I1, J2)))
               Y33 = 0.5_SP * (VY(I)+VY(NV(I1, J2)))
           !
               PUPX = PUPX + UU (I1, K) * (Y11-Y33)
               PUPY = PUPY + UU (I1, K) * (X33-X11)
               PVPX = PVPX + VV (I1, K) * (Y11-Y33)
               PVPY = PVPY + VV (I1, K) * (X33-X11)
               If (ISONB(I) /= 0) Then
                  If (J .Eq. 1) Then
                     PUPX = PUPX + UU (I1, K) * (VY(I)-Y11)
                     PUPY = PUPY + UU (I1, K) * (X11-VX(I))
                     PVPX = PVPX + VV (I1, K) * (VY(I)-Y11)
                     PVPY = PVPY + VV (I1, K) * (X11-VX(I))
                  End If
                  If (J .Eq. NTVE(I)) Then
                     PUPX = PUPX + UU (I1, K) * (Y11-VY(I))
                     PUPY = PUPY + UU (I1, K) * (VX(I)-X11)
                     PVPX = PVPX + VV (I1, K) * (Y11-VY(I))
                     PVPY = PVPY + VV (I1, K) * (VX(I)-X11)
                  End If
               End If
            End Do
        !
            PUPX = PUPX / ART1 (I)
            PUPY = PUPY / ART1 (I)
            PVPX = PVPX / ART1 (I)
            PVPY = PVPY / ART1 (I)
            TMP1 = PUPX ** 2 + PVPY ** 2
            TMP2 = 0.5_SP * (PUPY+PVPX) ** 2
            VISCOFF = Sqrt (TMP1+TMP2) * ART1 (I)
        !
            VISCOFH (I, K) = VISCOFF
        !
         End Do
      End Do
  !
      Return
End Subroutine VISCOF_H
!==============================================================================|
