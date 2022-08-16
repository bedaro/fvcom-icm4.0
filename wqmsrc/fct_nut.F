!fct_nut.F
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
!       subroutine FCT_NUT()
!
!==============================================================================|
!       FLUX CONTROL FOR NUTRIENTS                                             |
!==============================================================================|
!
Subroutine FCT_NUT
  !       #  if defined (WET_DRY)
  !
  !==============================================================================|
      Use MOD_LIMS, Only: MLOC, KBM1, NUMQBC, NUMPNT
  !
      Use MOD_TGE, Only: NTSN, NBSN
  !
      Use MOD_HYDROVARS, Only: DZ,DZ2D, THOUR
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: INFLOW_TYPE
      Use MOD_BCS, Only: INODEQ, INOPNT !!LOCAL NON-POINT SOURCE NODES
  !
      Use MOD_WQM, Only: C2, C2F, NAC, AC
      Use MOD_PREC, Only: SP
  !
      Use MOD_BCMAP, Only: IOBCN, I_OBC_N
  !
      Implicit None
      Real (SP) :: NUTMAX, NUTMIN
      Integer :: I, J, K, JCON, II
  !==============================================================================|
      Do JCON = 1, NAC
         II = AC (JCON)
         Do I = 1, MLOC
            If (IOBCN > 0) Then
               Do J = 1, IOBCN
                  If (I == I_OBC_N(J)) Go To 200
               End Do
            End If
        !
            If (NUMQBC > 0) Then
               Do J = 1, NUMQBC
                  If (INFLOW_TYPE == 'node') Then
                     If (I == INODEQ(J)) Go To 200
                  End If
               End Do
            End If
        !
            If (NUMPNT > 0) Then
               Do J = 1, NUMPNT
                  If (I == INOPNT(J)) Go To 200
               End Do
            End If
        !
            Do K = 1, KBM1
               NUTMAX = MAXVAL (C2(NBSN(I, 1:NTSN(I)), K, II))
               NUTMIN = MINVAL (C2(NBSN(I, 1:NTSN(I)), K, II))
           !
               If (K == 1) Then
                  NUTMAX = Max (NUTMAX, (C2(I, K, II)*DZ2D(I,K+1)+C2(I, &
                 & K+1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K+1)))
                  NUTMIN = Min (NUTMIN, (C2(I, K, II)*DZ2D(I,K+1)+C2(I, &
                 & K+1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K+1)))
               Else If (K == KBM1) Then
                  NUTMAX = Max (NUTMAX, (C2(I, K, II)*DZ2D(I,K-1)+C2(I, &
                 & K-1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K-1)))
                  NUTMIN = Min (NUTMIN, (C2(I, K, II)*DZ2D(I,K-1)+C2(I, &
                 & K-1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K-1)))
               Else
                  NUTMAX = Max (NUTMAX, (C2(I, K, II)*DZ2D(I,K-1)+C2(I, &
                 & K-1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K-1)), (C2(I, K, &
                 & II)*DZ2D(I,K+1)+C2(I, K+1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K+1)))
                  NUTMIN = Min (NUTMIN, (C2(I, K, II)*DZ2D(I,K-1)+C2(I, &
                 & K-1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K-1)), (C2(I, K, &
                 & II)*DZ2D(I,K+1)+C2(I, K+1, II)*DZ2D(I,K))/(DZ2D(I,K)+DZ2D(I,K+1)))
               End If

               If (NUTMIN-C2F(I, K, II) > 0.0_SP) C2F (I, K, II) = &
              & NUTMIN
               If (C2F(I, K, II)-NUTMAX > 0.0_SP) C2F (I, K, II) = &
              & NUTMAX

            End Do
        !
200         Continue
        !
         End Do
      End Do
      Return
  !       #  endif
End Subroutine FCT_NUT
!==============================================================================|
