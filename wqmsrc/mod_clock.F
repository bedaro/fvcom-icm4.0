      MODULE MOD_CLOCK
        USE MOD_PREC
        IMPLICIT NONE
        REAL(SP) :: TINIT,TFINAL,TCURRENT
        INTEGER :: COUNT_RATE,COUNT_MAX,COUNT,COUNT0
        real(SP) :: lastreport, lastheader
        integer, parameter :: reportinterval = 2 ! Reporting interval (seconds)

        CONTAINS

        !======================================================================!
        !  Initialize Timer => TINIT                                           !
        !======================================================================!
        SUBROUTINE START_CLOCK
           IMPLICIT NONE

           CALL SYSTEM_CLOCK (COUNT0, COUNT_RATE, COUNT_MAX)
           TINIT = FLOAT(COUNT0) / FLOAT(COUNT_RATE)
           lastreport = 0
           lastheader = 0

        END SUBROUTINE START_CLOCK

        !======================================================================!
        !  Retrieve System Time => TCURRENT                                    !
        !======================================================================!
        SUBROUTINE GET_CLOCK 
           IMPLICIT NONE  

           CALL SYSTEM_CLOCK(COUNT=COUNT)
           IF (COUNT < TINIT*FLOAT(COUNT_RATE)) COUNT = COUNT + COUNT_MAX
              TCURRENT = FLOAT(COUNT) / FLOAT(COUNT_RATE)

        END SUBROUTINE GET_CLOCK


        !======================================================================!
        !   Return a Time String DAYSdHRS:MIN:SEC from Number of Seconds       !
        !======================================================================!

        SUBROUTINE GETTIME(INSTRING,INSECS)

           IMPLICIT NONE
           INTEGER, INTENT(IN) :: INSECS 
           CHARACTER(LEN=13), INTENT(INOUT) :: INSTRING
           CHARACTER(LEN=4)  :: S0
           CHARACTER(LEN=2)  :: S1,S2,S3
           INTEGER :: DTCP,HTCP,MTCP,STCP
 
           DTCP = INSECS/(3600*24)
           HTCP = MOD(INSECS,(3600*24))/3600
           MTCP = MOD(INSECS,(3600))/60
           STCP = INSECS - (DTCP*3600*24 + HTCP*3600 + MTCP*60)

           IF(DTCP < 10000.)THEN
              IF(DTCP >= 1000.)THEN
                 WRITE(S0,"(I4)")INT(DTCP)
              ELSE IF(DTCP >= 100.)THEN
                 WRITE(S0,"(' ',I3)")INT(DTCP)
              ELSE IF(DTCP >= 10)THEN
                 WRITE(S0,"('  ',I2)")INT(DTCP)
              ELSE
                 WRITE(S0,"('   ',I1)")INT(DTCP)
              END IF
              IF(HTCP >= 10.)THEN
                 WRITE(S1,"(I2)")INT(HTCP)
              ELSE
                 WRITE(S1,"(' ',I1)")INT(HTCP)
              END IF
              IF(MTCP >= 10.)THEN
                 WRITE(S2,"(I2)")INT(MTCP)
              ELSE
                 WRITE(S2,"('0',I1)")INT(MTCP)
              END IF
              IF(STCP >= 10.)THEN
                 WRITE(S3,"(I2)")INT(STCP)
              ELSE
                 WRITE(S3,"('0',I1)")INT(STCP)
              END IF
              INSTRING = S0//"d"//S1//":"//S2//":"//S3
           ELSE
              INSTRING = "> 1000 DAYS"
           END IF
 
           RETURN
        END SUBROUTINE GETTIME
        !======================================================================!
        !   Report Calculation Speed and Time to Complete                      !
        !======================================================================!

        SUBROUTINE REPORT_TIME(IINT,ISTART,IEND,STIME)

           IMPLICIT NONE
           INTEGER, INTENT(IN) :: IINT,ISTART,IEND
           REAL(SP),INTENT(IN) :: STIME
           REAL(SP) :: TTCP,TAVE
           CHARACTER(LEN=13)  :: SIMTIME,FINTIME 
           CHARACTER(LEN=22) :: PCOMP
           INTEGER :: I,ICMP,ICMP2
 
           !
           !  CALCULATE COMPLETION PERCENTAGE GRAPHIC
           !
           IF(MOD(IINT,1) /= 0) RETURN
           ICMP = INT( 100.*FLOAT(IINT-ISTART+1)/FLOAT(IEND-ISTART+1)) 
           ICMP2 = ICMP/5
           PCOMP = " "
           PCOMP(1:1)  = "|"
           PCOMP(22:22) = "|"

           DO I=2,ICMP2+1
              PCOMP(I:I) = "="
           END DO
           !
           !  CALCULATE CURRENT TIME
           !
           CALL GET_CLOCK
           !
           !  CALCULATE AVERAGE TIME/ITERATION
           !
           TAVE = (TCURRENT-TINIT)/FLOAT(IINT-ISTART+1)
           CALL GETTIME(SIMTIME,INT(STIME))
           !
           !  CALCULATE TIME TO COMPLETION
           !
           TTCP = TAVE*(IEND-IINT)
           CALL GETTIME(FINTIME,INT(TTCP))
           !
           !  REPORT
           ! BDR: don't spam the screen with updates every timestep, rather
           ! only when a certain time has passed since the last one
           if(tcurrent - lastreport > reportinterval) then
              if(tcurrent - lastheader > reportinterval * 10) then
                 write(*,*)
                 write(*,102)
                 lastheader = tcurrent
              end if
              write(*,101)iint,simtime,fintime,tave,pcomp
              lastreport = tcurrent
           endif

        RETURN
 101 FORMAT(I7,3X,A13,3X,A13,3X,F8.4,2X,A22)  
 102 FORMAT(1X,"  IINT ",3X,"     JDAY    ",3X,"  FINISH IN ",3X,&
        "  SECS/IT ",1X,"  PERCENT COMPLETE   ")  
      END SUBROUTINE REPORT_TIME

      END MODULE MOD_CLOCK
