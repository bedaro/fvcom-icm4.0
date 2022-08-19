!mod_tge.F
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
Module MOD_TGE
  !
      Use MOD_LIMS, Only: MX_NBR_ELEM, MTLOC, NTLOC, MYID, MLOC, NLOC, &
     & NE, NCT, NCV_I, NCV, NISBCE_1, NISBCE_2, NISBCE_3
  !
      Use MOD_HYDROVARS, Only: ZERO, XC, YC, VX, VY
  !Took MOD_CONTROL out of MOD_HYDROVARS and put the used variables here
      Use MOD_CONTROL, Only: MSR
#  if defined (MULTIPROCESSOR)
      Use MOD_PAR, Only: NGID
# endif
  !
      Use MOD_BCMAP, Only: IOBCN, I_OBC_N
  !
      Use MOD_PREC, Only: SP
  !
  !
      Integer, Allocatable :: ISBCE (:)
      Integer, Allocatable :: ISONB (:)!!NODE MARKER = 0,1,2
      Integer, Allocatable :: NBE (:, :)!!INDICES OF ELMNT NEIGHBORS
      Integer, Allocatable :: NTVE (:)
      Integer, Allocatable :: NTSN (:)
  !
      Integer, Allocatable :: NV (:, :)!!NODE NUMBERING FOR ELEMENTS
      Integer, Allocatable :: IEC (:, :)!1
      Integer, Allocatable :: IENODE (:, :)!2
  !
      Real (SP), Allocatable :: XIJC (:)!3
      Real (SP), Allocatable :: YIJC (:)!4
      Real (SP), Allocatable :: DLTXC (:)!5
  !
      Real (SP), Allocatable :: DLTYC (:)!6
      Real (SP), Allocatable :: DLTXYC (:)!7
  !
      Real (SP), Allocatable :: SITAC (:)!8
      Integer, Allocatable :: ISBC (:)!9
  !
      Integer, Allocatable :: NBVE (:, :)!10
      Integer, Allocatable :: NBVT (:, :)!11
      Integer, Allocatable :: NBSN (:, :)!12
  !
      Integer, Allocatable :: NIEC (:, :)!13
      Integer, Allocatable :: NTRG (:)!14
  !
      Real (SP), Allocatable :: XIJE (:, :)!15
      Real (SP), Allocatable :: YIJE (:, :)!16
      Real (SP), Allocatable :: DLTXE (:)!17
      Real (SP), Allocatable :: DLTYE (:)!18
      Real (SP), Allocatable :: DLTXYE (:)!19
      Real (SP), Allocatable :: SITAE (:)!20
  !
      Integer, Allocatable :: LISBCE_1 (:)!21   !!LIST OF ELEMENTS WITH ISBCE=1
      Integer, Allocatable :: LISBCE_2 (:)!22   !!LIST OF ELEMENTS WITH ISBCE=2
      Integer, Allocatable :: LISBCE_3 (:)!23   !!LIST OF ELEMENTS WITH ISBCE=3
  !
      Real (SP), Allocatable :: EPOR (:)!24   !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
  !
  !
Contains
  !
  !Subroutines
  !	Subroutine TRIANGLE_GRID_EDGE()
  !	Subroutine TGE_DEALLOC()
  !
      Subroutine TRIANGLE_GRID_EDGE
    !
    !==============================================================================!
    !  This program is used to define the non-overlapped, unstructured             !
    !  triangular meshes used for flux computations. The mesh could be             !
    !  created using the commerical software called "sms8.0" or other              !
    !  mesh generation programs. The mesh file generated by sms8.0 can             !
    !  be directly used for this subroutine, while the mesh file                   !
    !  generated using other programs must be coverted the data format             !
    !  to meet the required format used here.                                      !
    !==============================================================================!
    !     variable list:							       !
    !  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)	               !
    !  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)	               !
    !  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i                         !
    !  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
    !  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
    !  cor(n)   :: cor(i) = f plane coriolis at element i                          !
    !                                                                              !
    !  nbe(n,3) :: nbe(i,1->3) = element index of 1->3 neighbors of element i      !
    !  isbce(n) :: flag if element is on the boundary, see below for values        !
    !  isonb(m) :: flag is node is on the boundary, see below for values           !
    !                                                                              !
    !  ntve(m)  :: the number of neighboring elements of node m                    !
    !  nbve(m,ntve(m)) :: nbve(i,1->ntve(i)) = ntve elements containing node i     !
    !  nbvt(m,ntve(m)) :: nbvt(i,j) = the node number of node i in element         !
    !                     nbve(i,j) (has a value of 1,2,or 3)                      !
    !                                                                              !
    !   ne       ::   number of unique element edges                               !
    !  iec(ne,2) :: nbe(i,1->2) cell number of cell(s) connected to edge i         !
    !   isbc(ne) :: flag marking edge property                                     !
    !              isbc(i) = 0:  element edge i is in the interior                 !
    !              isbc(i) = 1:  element edge i is on the boundary                 !
    !ienode(ne,2):: ienode(i,1->2) node numbers at each end of element edge i      !
    !  xijc(ne)  :: xijc(i) = x-coordinate of mid point of element edge i          !
    !  yijc(ne)  :: yijc(i) = y-coordinate of mid point of element edge i          !
    !  dltxyc(ne):: dltxyc(i) = length of element edge i                           !
    !  dltxc(ne) :: dltxc(i) = deltax (x-projection) of element edge i             !
    !  dltyc(ne) :: dltyc(i) = deltay (y-projection) of element edge i             !
    !  sitac(ne) :: sitac(i) =  arctg(dltyc,dltxc) (angle of inclination of edge)  !
    !                                                                              !
    !==============================================================================!
    !     classification of the triangles nodes, and edges                         !
    !                                                                              !
    !     isonb(i)=0:  node in the interior computational domain                   !
    !     isonb(i)=1:  node on the solid boundary                                  !
    !     isonb(i)=2:  node on the open boundary                                   !
    !                                                                              !
    !     isbce(i)=0:  element in the interior computational domain                !
    !     isbce(i)=1:  element on the solid boundary                               !
    !     isbce(i)=2:  element on the open boundary                                !
    !     isbce(i)=3:  element with 2 solid boundary edges                         !
    !                                                                              !
    !      isbc(i)=0:  element edge in interior                                    !
    !      isbc(i)=1:  element edge on boundary                                    !
    !==============================================================================!
    !
    !
    !==============================================================================|
    !   FIND NEIGHBORING ELEMENTS, MARK BOUNDARY NODES AND ELEMENTS                |
    !									       |
    !   NBE(N,3) :  NBE(I,J) = IDENTITY OF NBOR ELMNT TO TRI I ON EDGE J           |
    !   IBCE(N)  :  DESCRIBED IN SUBROUTINE HEADING			               |	
    !   ISONB(M):  DESCRIBED IN SUBROUTINE HEADING			               |	
    !==============================================================================|
    !
    !
	  
	  	 USE MOD_FILEINFO, ONLY:IREPORT
         Implicit None
         Integer, Allocatable, Dimension (:, :) :: TEMP, TEMP2, NB_TMP, &
        & ISET
         Integer I, J, II, JJ, NTMP, NCNT, INEY, NFLAG
    !DIF1,DIF2,DIF3, !!
         Integer ITMP1, ITMP2, ITMP3, JN, JJB, IBCETMP, NCTMP, NCETMP, &
        & NPT
         Integer, Allocatable :: CELLS (:, :), CELLCNT (:), NBET (:, :)
         Real (SP) DTMP
         Integer N1, N2, N3, J1, J2, J3
    !
    !----------------------------REPORT--------------------------------------------!
         If (MSR) WRITE (*,*) '!'
         If (MSR) WRITE (*,*) '!           SETTING UP TRIS/ELEMENTS/EDG&
        &ES/CVS          '
         If (MSR) WRITE (*,*) '!'
    !----------------------------INITIALIZE----------------------------------------!
    !
         Allocate (ISBCE(0:NTLOC))
         ISBCE = 0
         Allocate (ISONB(0:MTLOC))
         ISONB = 0 !!NODE MARKER = 0,1,2
         Allocate (NBE(0:NTLOC, 3))
         NBE = 0 !!INDICES OF ELMNTNEIGHBORS
    !
    !
    !----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!
    !
         Allocate (NBET(NTLOC, 3))
         NBET = 0
         Allocate (CELLS(MTLOC, 50))
         CELLS = 0
         Allocate (CELLCNT(MTLOC))
         CELLCNT = 0
         Do I = 1, NTLOC
            N1 = NV (I, 1)
            CELLCNT (N1) = CELLCNT (N1) + 1
            N2 = NV (I, 2)
            CELLCNT (N2) = CELLCNT (N2) + 1
            N3 = NV (I, 3)
            CELLCNT (N3) = CELLCNT (N3) + 1
            CELLS (NV(I, 1), CELLCNT(N1)) = I
            CELLS (NV(I, 2), CELLCNT(N2)) = I
            CELLS (NV(I, 3), CELLCNT(N3)) = I
         End Do
         If (maxval(CELLCNT) > 50) WRITE (*,*) 'bad', maxval (CELLCNT)
         Do I = 1, NTLOC
            N1 = NV (I, 1)
            N2 = NV (I, 2)
            N3 = NV (I, 3)
            Do J1 = 1, CELLCNT (N1)
			IF(CELLS(N1,J1) /= I) THEN	
               Do J2 = 1, CELLCNT (N2)
                  IF(CELLS(N1,J1) == CELLS(N2,J2))NBE(I,3) = CELLS(N1,J1)
               End Do
			ENDIF
            End Do
            Do J2 = 1, CELLCNT (N2)
			IF( CELLS(N2,J2) /= I) THEN
               Do J3 = 1, CELLCNT (N3)
                 IF(CELLS(N2,J2) == CELLS(N3,J3))NBE(I,1) = CELLS(N2,J2)
               End Do
			ENDIF
            End Do
            Do J1 = 1, CELLCNT (N1)
			IF( CELLS(N1,J1) /=I) THEN
               Do J3 = 1, CELLCNT (N3)
                  IF(CELLS(N1,J1) == CELLS(N3,J3))NBE(I,2) = CELLS(N3,J3)
               End Do
			ENDIF
            End Do
         End Do
         Deallocate (CELLS, CELLCNT, NBET)!LB added NBET to the deallocation statement
         If (MSR) WRITE (*,*) '!  NEIGHBOR FINDING      :    COMPLETE'
    !
    !--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
    !
         NFLAG = 0
         Do I = 1, NTLOC
            If (SUM(NBE(I, 1:3)) == 0) Then
               NFLAG = 1
               Write (*,*) 'ELEMENT ', I, ' AT ', XC (I), YC (I), ' HAS&
              & NO NEIGHBORS'
               Call PSTOP
            End If
         End Do
         If (NFLAG == 1) Call PSTOP
    !
    !
    !----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
    !
    !
         Do I = 1, NTLOC
            If (Min(NBE(I, 1), NBE(I, 2), NBE(I, 3)) == 0) Then !!ELEMENT ON BOUNDARY
               ISBCE (I) = 1
               If (NBE(I, 1) == 0) Then
                  ISONB (NV(I, 2)) = 1
                  ISONB (NV(I, 3)) = 1
               End If
               If (NBE(I, 2) == 0) Then
                  ISONB (NV(I, 1)) = 1
                  ISONB (NV(I, 3)) = 1
               End If
               If (NBE(I, 3) == 0) Then
                  ISONB (NV(I, 1)) = 1
                  ISONB (NV(I, 2)) = 1
               End If
            End If
         End Do
         If (MSR) WRITE (*,*) '!  ISONB SETTING         :    COMPLETE'
    !
    !==============================================================================|
    !             DEFINE NTVE, NBVE, NBVT                                          !
    !                                                                              !
    ! ntve(1:m):           total number of the surrounding triangles               !
    !                      connected to the given node                             !
    ! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
    !                      triangles with a common node (counted clockwise)        !
    ! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
    !                      each individual surrounding triangle(counted            !
    !                      clockwise)                                              !
    ! ntsn(1:m):           total number of surrounding nodes                       !
    ! nbsn(1:m, ntsn):     the identification number of surrounding nodes          !
    !                      (counted clockwise)                                     !
    ! nbse(1:m,2*ntsn):    the identification number of control volume s           !
    !                      edges between two neighbor nodes                        !
    !==============================================================================|
    !
    !
    !----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
    !
         MX_NBR_ELEM = 0
         Do I = 1, MLOC
            NCNT = 0
            Do J = 1, NTLOC
          !       IF( (NV(J,1)-I)*(NV(J,2)-I)*(NV(J,3)-I)==0) NCNT = NCNT + 1
          !       IF( (NV(J,1)-I) ==0 .OR. (NV(J,2)-I) ==0 .OR. (NV(J,3)-I)==0) NCNT = NCNT + 1
               If (FLOAT(NV(J, 1)-I)*FLOAT(NV(J, 2)-I)*FLOAT(NV(J, &
              & 3)-I) == 0.0_SP) NCNT = NCNT + 1
            End Do
            MX_NBR_ELEM = Max (MX_NBR_ELEM, NCNT)
         End Do
    !   WRITE(*,*) 'MAXIMUM NUMBER OF NEIGHBOR ELEMENTS',MX_NBR_ELEM
    !
    !
    !----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
    !
    !
         Allocate (NBVE(MLOC, MX_NBR_ELEM+1))
         Allocate (NBVT(MLOC, MX_NBR_ELEM+1))
         Allocate (NBSN(MLOC, MX_NBR_ELEM+3))!!MHB
    !
    !--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
    !--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
    !--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
    !
    !
    !
         Allocate (NTVE(0:MTLOC))
         NTVE = 0
         Do I = 1, MLOC
            NCNT = 0
            Do J = 1, NTLOC
          !       IF( (NV(J,1)-I) == 0 .OR.  (NV(J,2)-I) == 0 .OR. (NV(J,3)-I) == 0)THEN
               If (FLOAT(NV(J, 1)-I)*FLOAT(NV(J, 2)-I)*FLOAT(NV(J, &
              & 3)-I) == 0.0_SP) Then
                  NCNT = NCNT + 1
                  NBVE (I, NCNT) = J
                  If ((NV(J, 1)-I) == 0) NBVT (I, NCNT) = 1
                  If ((NV(J, 2)-I) == 0) NBVT (I, NCNT) = 2
                  If ((NV(J, 3)-I) == 0) NBVT (I, NCNT) = 3
               End If
            End Do
            NTVE (I) = NCNT
         End Do
    !
    !
    !--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
    !--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
    !--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!
    !
    !
         Allocate (NTSN(MTLOC))
         NTSN = 0
    !
         Allocate (NB_TMP(MLOC, MX_NBR_ELEM+1))
    !
         Do I = 1, MLOC
            If (ISONB(I) == 0) Then
               NB_TMP (1, 1) = NBVE (I, 1)
               NB_TMP (1, 2) = NBVT (I, 1)
               Do J = 2, NTVE (I) + 1
                  II = NB_TMP (J-1, 1)
                  JJ = NB_TMP (J-1, 2)
                  NB_TMP (J, 1) = NBE (II, JJ+1-Int((JJ+1)/4)*3)
                  JJ = NB_TMP (J, 1)
                  If ((NV(JJ, 1)-I) == 0) NB_TMP (J, 2) = 1
                  If ((NV(JJ, 2)-I) == 0) NB_TMP (J, 2) = 2
                  If ((NV(JJ, 3)-I) == 0) NB_TMP (J, 2) = 3
               End Do
          !
               Do J = 2, NTVE (I) + 1
                  NBVE (I, J) = NB_TMP (J, 1)
               End Do
          !
               Do J = 2, NTVE (I) + 1
                  NBVT (I, J) = NB_TMP (J, 2)
               End Do
          !
               NTMP = NTVE (I) + 1
               If (NBVE(I, 1) /= NBVE(I, NTMP)) Then
                  Write (*,*) MYID, I, 'NBVE(I) NOT CORRECT!!'
                  Call PSTOP
               End If
               If (NBVT(I, 1) /= NBVT(I, NTMP)) Then
                  Write (*,*) I, 'NBVT(I) NOT CORRECT!!'
                  Call PSTOP
               End If
          !
               NTSN (I) = NTVE (I)
          !
               Do J = 1, NTSN (I)
                  II = NBVE (I, J)
                  JJ = NBVT (I, J)
                  NBSN (I, J) = NV (II, JJ+1-Int((JJ+1)/4)*3)
               End Do
          !
               NTSN (I) = NTSN (I) + 1
               NBSN (I, NTSN(I)) = NBSN (I, 1)
          !
            Else
               JJB = 0
          !
               Do J = 1, NTVE (I)
                  JJ = NBVT (I, J)
                  If (NBE(NBVE(I, J), JJ+2-Int((JJ+2)/4)*3) == 0) Then
                     JJB = JJB + 1
                     NB_TMP (JJB, 1) = NBVE (I, J)
                     NB_TMP (JJB, 2) = NBVT (I, J)
                  End If
               End Do
          !
               If (JJB /= 1) Then
                  Write (*,*) 'ERROR IN ISONB !,I,J', I, J
                  call pstop
               End If
          !
               Do J = 2, NTVE (I)
                  II = NB_TMP (J-1, 1)
                  JJ = NB_TMP (J-1, 2)
                  NB_TMP (J, 1) = NBE (II, JJ+1-Int((JJ+1)/4)*3)
                  JJ = NB_TMP (J, 1)
                  If ((NV(JJ, 1)-I) == 0) NB_TMP (J, 2) = 1
                  If ((NV(JJ, 2)-I) == 0) NB_TMP (J, 2) = 2
                  If ((NV(JJ, 3)-I) == 0) NB_TMP (J, 2) = 3
               End Do
          !
               Do J = 1, NTVE (I)
                  NBVE (I, J) = NB_TMP (J, 1)
                  NBVT (I, J) = NB_TMP (J, 2)
               End Do
          !
               NBVE (I, NTVE(I)+1) = 0
               NTSN (I) = NTVE (I) + 1
               NBSN (I, 1) = I
          !
               Do J = 1, NTSN (I) - 1
                  II = NBVE (I, J)
                  JJ = NBVT (I, J)
                  NBSN (I, J+1) = NV (II, JJ+1-Int((JJ+1)/4)*3)
               End Do
          !
               J = NTSN (I)
               II = NBVE (I, J-1)
               JJ = NBVT (I, J-1)
               NBSN (I, J+1) = NV (II, JJ+2-Int((JJ+2)/4)*3)
               NTSN (I) = NTSN (I) + 2
               NBSN (I, NTSN(I)) = I
            End If
         End Do
         Deallocate (NB_TMP)
         If (MX_NBR_ELEM+3-maxval(NTSN) < 0) Then
            Write (*,*) 'CHECK NTSN/NBSN', maxval (NTSN), MX_NBR_ELEM + &
           & 3
            Call PSTOP
         End If
         If (MSR) WRITE (*,*) '!  NBVE/NBVT             :    COMPLETE'
    !
    !
    !==============================================================================!
    !  Define the parameters of each triangular edge                               !
    !                                                                              !
    !  ne           :    number of unique element edges                            !
    !  iec(1:ne,1:2):    counting number identifying two connected cells           !
    !  isbc(1:ne):       0: triangle s edge in the interior                        !
    !                    1: triangle s edge on the boundary                        !
    !  ienode(1:ne,1:2): the identification number of two end points of a          !
    !                    edge                                                      !
    !  xijc(1:ne):       the x-coordinate location of the middle points            !
    !                    of a edge                                                 !
    !  yijc(1:ne):       the y-coordinate location of the middle points            !
    !                    of a edge                                                 !
    !  dltxyc(1:ne):     length of the edge                                        !
    !  dltxc(1:ne):      vx(ienode(i,2))-vx(idnode(i,1))                           !
    !  dltyc(1:ne):      vy(ienode(i,2))-vy(idnode(i,1))                           !
    !  sitac(1:ne):      arctg(dltyc,dltxc)                                        !
    !==============================================================================!
    !
         Allocate (ISET(NTLOC, 3), TEMP((NTLOC)*3, 2), TEMP2((NTLOC)*3, &
        & 2))
         ISET = 0
         NE = 0
         TEMP = 0
         TEMP2 = 0
         Do I = 1, NTLOC
            Do J = 1, 3
               If (ISET(I, J) == 0) Then
                  NE = NE + 1
                  INEY = NBE (I, J)
                  ISET (I, J) = 1
                  Do JN = 1, 3
                     If (I == NBE(INEY, JN)) ISET (INEY, JN) = 1
                  End Do
                  TEMP (NE, 1) = I
                  TEMP (NE, 2) = INEY
                  TEMP2 (NE, 1) = NV (I, J+1-Int((J+1)/4)*3)
                  TEMP2 (NE, 2) = NV (I, J+2-Int((J+2)/4)*3)
               End If
            End Do
         End Do
         Deallocate (ISET)
    !
    !--ALLOCATE ARRAYS REQUIRING NUMBER OF EDGES-----------------------------------!
    !
         Allocate (IEC(NE, 2))
         Allocate (IENODE(NE, 2))
         Allocate (XIJC(NE))
         Allocate (YIJC(NE))
         Allocate (DLTXC(NE))
         Allocate (DLTYC(NE))
         Allocate (DLTXYC(NE))
         Allocate (SITAC(NE))
         Allocate (ISBC(NE))
    !
    !
         IEC (:, 1) = TEMP (1:NE, 1)
         IEC (:, 2) = TEMP (1:NE, 2)
         IENODE (:, 1) = TEMP2 (1:NE, 1)
         IENODE (:, 2) = TEMP2 (1:NE, 2)
    !
    !
         Deallocate (TEMP, TEMP2)
    !
    !
    !
    !------MARK ELEMENT EDGES THAT ARE ON THE BOUNDARY-----------------------------!
    !
         ISBC = 0
         Do I = 1, NE
            If ((IEC(I, 1) == 0) .Or. (IEC(I, 2) == 0)) ISBC (I) = 1
         End Do
    !
    !
    !------CALCULATE ELEMENT EDGE METRICS------------------------------------------!
    !
         Do I = 1, NE
            DLTXC (I) = VX (IENODE(I, 2)) - VX (IENODE(I, 1))
            DLTYC (I) = VY (IENODE(I, 2)) - VY (IENODE(I, 1))
            XIJC (I) = (VX(IENODE(I, 1))+VX(IENODE(I, 2))) / 2.0
            YIJC (I) = (VY(IENODE(I, 1))+VY(IENODE(I, 2))) / 2.0
            DLTXYC (I) = Sqrt (DLTXC(I)**2+DLTYC(I)**2)
            SITAC (I) = Atan2 (DLTYC(I), DLTXC(I))
         End Do
    !
         If (MSR) WRITE (*,*) '!  EDGE SETUP            :    COMPLETE'
    !
    !
    !==============================================================================!
    !  read triangular mesh parameters on open boundary :                          !
    !  iobce:   number of open boundary cells.                                     !
    !  isbcn:   number of open boundary nodes.                                     !
    !  i_obc_e: counter number of open boundary cells   !never used                !
    !  i_obc_n: counter number of open boundary nodes                              !
    !==============================================================================!
    !
    !
    !----TRAVERSE  BOUNDARY NODE NUMBERS AND SET ISONB(NODE)=2---------------------!
    !
         Do I = 1, IOBCN
            ISONB (I_OBC_N(I)) = 2
         End Do
    !
    !
    !----DETERMINE IF ELEMENT IS ON OPEN BOUNDARY (CONTAINS EDGE ON OPEN BOUNDARY)-!
    !
         IBCETMP = 0
         Do I = 1, NLOC
            ITMP1 = ISONB (NV(I, 1))
            ITMP2 = ISONB (NV(I, 2))
            ITMP3 = ISONB (NV(I, 3))
       !
            If (SUM(ISONB(NV(I, 1:3))) == 4) Then
               ISBCE (I) = 2
               IBCETMP = IBCETMP + 1
            Else If (SUM(ISONB(NV(I, 1:3))) > 4) Then
#  if defined (MULTIPROCESSOR)
               Write (*,*) 'SORRY, THE BOUNDARY CELL', NGID (I), 'IS NO&
              &T GOOD FOR MODEL.'
#  else
               Write (*,*) 'SORRY, THE BOUNDARY CELL', I, 'IS NOT GOOD &
              &FOR MODEL.'
#  endif
               Write (*,*) 'IT HAS EITHER TWO SIDES OF OPEN BOUNDARY OR&
              & ONE OPEN BOUNDARY'
               Write (*,*) 'AND ONE SOLID BOUNDARY. PLEASE CHECK AND MO&
              &DIFIED IT.'
               Write (*,*) 'THIS MESSAGE IS IN SUBROUTINE TRIANGLE_GRID&
              &_EDGE (MOD_TGE.F)'
               Write (*,*) 'STOP RUNNING...'
               Call PSTOP
            End If
         End Do
      !
         Do I = 1, NTLOC
            If ((NBE(I, 1)+NBE(I, 2)+NBE(I, 3) == 0) .And. (ISBCE(I) /= &
           & 2)) ISBCE (I) = 3
            If ((NBE(I, 1)+NBE(I, 2) == 0) .And. (ISBCE(I) /= 2)) ISBCE &
           & (I) = 3
            If ((NBE(I, 2)+NBE(I, 3) == 0) .And. (ISBCE(I) /= 2)) ISBCE &
           & (I) = 3
            If ((NBE(I, 1)+NBE(I, 3) == 0) .And. (ISBCE(I) /= 2)) ISBCE &
           & (I) = 3
         End Do
      !
      !==============================================================================!
      !  xije(1:nc,1:2):  the x coordinate locations of starting and ending          !
      !                   points of the control volumes edges                        !
      !  yije(1:nc,1:2):  the y coordinate locations of starting and ending          !
      !                   points of the control volumes edges                        !
      !  niec(1:nc,1:2):  the counting number of left and right nodes                !
      !                   conected to this control volumes edge from                 !
      !                   starting point to ending point                             !
      !  dltxe(1:nc):     the x distance of individual edges                         !
      !  dltye(1:nc)      the y distance of individual edges                         !
      !  dltxye(1:nc):    the length of individual edges                             !
      !  ntrg(1:nc)  :    element associated with this control volume edge           !
      !==============================================================================!
      !
         NCT = NTLOC * 3
         Allocate (XIJE(NCT, 2))
         XIJE = ZERO
         Allocate (YIJE(NCT, 2))
         YIJE = ZERO
         Allocate (NIEC(NCT, 2))
         NIEC = 0
         Allocate (NTRG(NCT))
         NTRG = 0
         Allocate (DLTXE(NCT))
         DLTXE = ZERO
         Allocate (DLTYE(NCT))
         DLTYE = ZERO
         Allocate (DLTXYE(NCT))
         DLTXYE = ZERO
         Allocate (SITAE(NCT))
         SITAE = ZERO
         NCTMP = 0
         NCETMP = 0
      !
         Do I = 1, NE
            If (ISBC(I) == 0) Then
               If (IEC(I, 1) <= NLOC) Then
                  NCTMP = NCTMP + 1
                  NPT = NCTMP
               Else
                  NCETMP = NCETMP + 1
                  NPT = NCETMP + (3*NLOC)
               End If
               XIJE (NPT, 1) = XC (IEC(I, 1))
               YIJE (NPT, 1) = YC (IEC(I, 1))
               XIJE (NPT, 2) = XIJC (I)
               YIJE (NPT, 2) = YIJC (I)
               NIEC (NPT, 1) = IENODE (I, 1)
               NIEC (NPT, 2) = IENODE (I, 2)
               NTRG (NPT) = IEC (I, 1)
               DLTXE (NPT) = XIJE (NPT, 2) - XIJE (NPT, 1)
               DLTYE (NPT) = YIJE (NPT, 2) - YIJE (NPT, 1)
               DTMP = DLTXE (NPT) * DLTXE (NPT) + DLTYE (NPT) * DLTYE &
              & (NPT)
               DLTXYE (NPT) = Sqrt (DTMP)
               SITAE (NPT) = Atan2 (DLTYE(NPT), DLTXE(NPT))
            !
               If (IEC(I, 2) <= NLOC) Then
                  NCTMP = NCTMP + 1
                  NPT = NCTMP
               Else
                  NCETMP = NCETMP + 1
                  NPT = NCETMP + (3*NLOC)
               End If
               XIJE (NPT, 1) = XC (IEC(I, 2))
               YIJE (NPT, 1) = YC (IEC(I, 2))
               XIJE (NPT, 2) = XIJC (I)
               YIJE (NPT, 2) = YIJC (I)
               NIEC (NPT, 1) = IENODE (I, 2)
               NIEC (NPT, 2) = IENODE (I, 1)
               NTRG (NPT) = IEC (I, 2)
               DLTXE (NPT) = XIJE (NPT, 2) - XIJE (NPT, 1)
               DLTYE (NPT) = YIJE (NPT, 2) - YIJE (NPT, 1)
               DTMP = DLTXE (NPT) * DLTXE (NPT) + DLTYE (NPT) * DLTYE &
              & (NPT)
               DLTXYE (NPT) = Sqrt (DTMP)
               SITAE (NPT) = Atan2 (DLTYE(NPT), DLTXE(NPT))
            Else If (ISBC(I) == 1) Then
               If (IEC(I, 1) <= NLOC) Then
                  NCTMP = NCTMP + 1
                  NPT = NCTMP
               Else
                  NCETMP = NCETMP + 1
                  NPT = NCETMP + (3*NLOC)
               End If
               If (IEC(I, 1) == 0) Then
                  Write (*,*) I, 'IEC(I,1)===0'
                  Call PSTOP
               End If
               XIJE (NPT, 1) = XC (IEC(I, 1))
               YIJE (NPT, 1) = YC (IEC(I, 1))
               XIJE (NPT, 2) = XIJC (I)
               YIJE (NPT, 2) = YIJC (I)
               NIEC (NPT, 1) = IENODE (I, 1)
               NIEC (NPT, 2) = IENODE (I, 2)
               NTRG (NPT) = IEC (I, 1)
               DLTXE (NPT) = XIJE (NPT, 2) - XIJE (NPT, 1)
               DLTYE (NPT) = YIJE (NPT, 2) - YIJE (NPT, 1)
               DTMP = DLTXE (NPT) * DLTXE (NPT) + DLTYE (NPT) * DLTYE &
              & (NPT)
               DLTXYE (NPT) = Sqrt (DTMP)
               SITAE (NPT) = Atan2 (DLTYE(NPT), DLTXE(NPT))
            Else
               Write (*,*) 'ISBC(I) NOT CORRECT, I==', I
               Call PSTOP
            End If
         End Do
      !
         NCV_I = NCTMP
         NCV = NCETMP + NCTMP
      !
         If (NCV /= 3*(NTLOC)) Then
            Write (*,*) 'NCV IS NOT CORRECT, PLEASE CHECK THE SETUP'
            Call PSTOP
         End If
         If (NCV_I /= 3*NLOC) Then
            Write (*,*) 'NCV_I IS NOT CORRECT, PLEASE CHECK THE SETUP'
            Call PSTOP
         End If
      !
         Do I = 1, NCV_I
            If (NIEC(I, 1) > MLOC .Or. NIEC(I, 2) > MLOC) Then
               Write (*,*) 'PROBLEMS:', NIEC (I, 1), NIEC (I, 2), MLOC
               Call PSTOP
            End If
         End Do
      !
      !==============================================================================!
      !  nisbce_1/nisbce_2/nisbce_3:  number of elements with isbce of 1,2,3         !
      !  lisbce_1/lisbce_2/lisbce_3:  list of elements with isbce of 1,2,3           !
      !  epor                      :  element porosity (=0 if isbce = 2)             !
      !==============================================================================!
      !
      !
      !--COUNT NUMBER OF ELEMENTS OF EACH TYPE (ISBCE=1,2,3)-------------------------!
      !
         NISBCE_1 = 0
         NISBCE_2 = 0
         NISBCE_3 = 0
         Do I = 1, NLOC
            If (ISBCE(I) == 1) NISBCE_1 = NISBCE_1 + 1
            If (ISBCE(I) == 2) NISBCE_2 = NISBCE_2 + 1
            If (ISBCE(I) == 3) NISBCE_3 = NISBCE_3 + 1
         End Do
      !
      !
      !--ALLOCATE ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3---------------------!
      !
         If (NISBCE_1 > 0) Then
            Allocate (LISBCE_1(NISBCE_1))
         Else
              WRITE(IREPORT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=1'
         End If
      !
         If (NISBCE_2 > 0) Then
            Allocate (LISBCE_2(NISBCE_2))
         Else
              WRITE(IREPORT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=2'
         End If
      !
         If (NISBCE_3 > 0) Then
            Allocate (LISBCE_3(NISBCE_3))
         Else
              WRITE(IREPORT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=3'
         End If
      !
      !
      !--LOAD ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3--------------------------!
      !
         NISBCE_1 = 0
         NISBCE_2 = 0
         NISBCE_3 = 0
         Do I = 1, NLOC
            If (ISBCE(I) == 1) Then
               NISBCE_1 = NISBCE_1 + 1
               LISBCE_1 (NISBCE_1) = I
            End If
            If (ISBCE(I) == 2) Then
               NISBCE_2 = NISBCE_2 + 1
               LISBCE_2 (NISBCE_2) = I
            End If
            If (ISBCE(I) == 3) Then
               NISBCE_3 = NISBCE_3 + 1
               LISBCE_3 (NISBCE_3) = I
            End If
         End Do
      !
      !
      !--SET ELEMENT POROSITY---------------------------------------------------------!
      !
         Allocate (EPOR(0:NTLOC))
         EPOR = 1.0_SP
         Do I = 1, NLOC
            If (ISBCE(I) == 2) EPOR (I) = 0.0_SP
         End Do
      !
         If (MSR) WRITE (*,*) '!  NISBCE/LISBCE/EPOR    :    COMPLETE'
      !
         If (MSR) WRITE (*,*) '!  TRIS/EDGES/CVOLS      :    COMPLETE'
      !
         Return
      End Subroutine TRIANGLE_GRID_EDGE
    !
      Subroutine TGE_ALLOC
      !Moved here from pdomdec.F
         Allocate (NV(0:NTLOC, 4))
         NV = 0 !!NODE NUMBERING FOR ELEMENTS
      !
      End Subroutine TGE_ALLOC
    !
      Subroutine TGE_DEALLOC
      !
      !Deallocate variables defined in MOD_TGE
      !
         If (ALLOCATED(ISBCE)) DEALLOCATE (ISBCE)
         If (ALLOCATED(ISONB)) DEALLOCATE (ISONB)
         If (ALLOCATED(NBE)) DEALLOCATE (NBE)
         If (ALLOCATED(NTVE)) DEALLOCATE (NTVE)
         If (ALLOCATED(NTSN)) DEALLOCATE (NTSN)
      !
         If (ALLOCATED(NV)) DEALLOCATE (NV)
         If (ALLOCATED(IEC)) DEALLOCATE (IEC)
         If (ALLOCATED(IENODE)) DEALLOCATE (IENODE)
      !
         If (ALLOCATED(XIJC)) DEALLOCATE (XIJC)
         If (ALLOCATED(YIJC)) DEALLOCATE (YIJC)
         If (ALLOCATED(DLTXC)) DEALLOCATE (DLTXC)
         If (ALLOCATED(DLTYC)) DEALLOCATE (DLTYC)
         If (ALLOCATED(DLTXYC)) DEALLOCATE (DLTXYC)
         If (ALLOCATED(SITAC)) DEALLOCATE (SITAC)
         If (ALLOCATED(ISBC)) DEALLOCATE (ISBC)
      !
         If (ALLOCATED(NBVE)) DEALLOCATE (NBVE)
         If (ALLOCATED(NBVT)) DEALLOCATE (NBVT)
         If (ALLOCATED(NBSN)) DEALLOCATE (NBSN)
      !
         If (ALLOCATED(NIEC)) DEALLOCATE (NIEC)
         If (ALLOCATED(NTRG)) DEALLOCATE (NTRG)
      !
         If (ALLOCATED(XIJE)) DEALLOCATE (XIJE)
         If (ALLOCATED(YIJE)) DEALLOCATE (YIJE)
         If (ALLOCATED(DLTXE)) DEALLOCATE (DLTXE)
         If (ALLOCATED(DLTYE)) DEALLOCATE (DLTYE)
         If (ALLOCATED(DLTXYE)) DEALLOCATE (DLTXYE)
         If (ALLOCATED(SITAE)) DEALLOCATE (SITAE)
         If (ALLOCATED(EPOR)) DEALLOCATE (EPOR)
      !
      !LB: these ones were missing
         If (ALLOCATED(LISBCE_1)) DEALLOCATE (LISBCE_1)
         If (ALLOCATED(LISBCE_2)) DEALLOCATE (LISBCE_2)
         If (ALLOCATED(LISBCE_3)) DEALLOCATE (LISBCE_3)
      !
      !
      !
         Return
      !
      End Subroutine TGE_DEALLOC
    !
    !==============================================================================!
    !
End Module MOD_TGE
