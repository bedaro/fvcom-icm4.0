!fixformat.F
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
Program FixFormat
  !Simple progam to fix the format of hand prepared W2 input files
  !So that each column is precisely 8 columns
  !
      Implicit None
      Integer, Parameter :: ISHORT = Selected_Int_Kind (8), ILONG = &
     & Selected_Int_Kind (16), SPREC = Selected_Real_Kind (6, 30), &
     & DPREC = Selected_Real_Kind (12, 300), STRLEN_SHORT = 72, &
     & STRLEN_LONG = 1024 !Long string length
  !
      Integer :: i
  !
      Integer (ISHORT) :: NCOL, NHD !number of header lines
  !
      Character (Len=STRLEN_LONG) :: FPREFIX, FAFTERFIX !name of the after fix (output)
  !
      Character (Len=STRLEN_LONG) :: FMTSTR !string to hold format
  !
      Character (Len=STRLEN_LONG) :: HEADERLINE !header line
  !
      Real (SPREC), Allocatable :: ROWFIELDS (:)
  !
  !
      Write (*, '(A)', Advance='NO') 'Please give the original file nam&
     &e (e.g. Cdt_br1.npt):'
      Read (*,*) FPREFIX
  !
      Write (*, '(A)', Advance='NO') 'Please give the output file name &
     &(e.g. Cdt_br1.npt.fixed):'
      Read (*,*) FAFTERFIX
  !
      Write (*, '(A)', Advance='NO') 'Please give the number of columns&
     & in the file (e.g. 62):'
      Read (*,*) NCOL
  !
      Write (*, '(A)', Advance='NO') 'Please give the number of header &
     &lines in the file (e.g. 3):'
      Read (*,*) NHD
  !
  !Allocate the row
      Allocate (ROWFIELDS(1:NCOL))
      ROWFIELDS = 0.0
  !
  !Open the existing file and
  !
      Open (Unit=1, File=FPREFIX, Status='OLD')
  !
  !Open the new file for write
      Open (Unit=2, File=FAFTERFIX, Status='UNKNOWN')
  !
  !read and output the header lines
  !
  !Loop through the BII file and read all records, then output to the CSVFILE
      Do i = 1, NHD
         Read (1, '(A1024)', END=10) HEADERLINE
         Write (2, '(A1024)') TRIM (HEADERLINE)
      End Do
  !
      Write (FMTSTR, '(I0)') NCOL
      FMTSTR = "(" // TRIM (FMTSTR) // '(F12.5,1X)' // ")" !make up some thing like (15(F12.5,1X)) if NCOL=15
      Write (*,*) TRIM (FMTSTR)
  !
      Do WHILE (.True.)
         Read (1,*, END=10) (ROWFIELDS(i), i=1, NCOL)
         Write (2, TRIM(FMTSTR)) (ROWFIELDS(i), i=1, NCOL)
      End Do
  !
      Deallocate (ROWFIELDS)
  !
10    Write (*,*) 'Done! Enter to quit'
      Read (*,*)
      Close (1)
      Close (2)
  !
End Program FixFormat
