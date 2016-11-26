C
C           *%%%%%*
C        *%%%%%%%%%%%*
C      *%%%%%/\|/\%%%%%*
C    *%%%%%%|--2--|%%%%%%*   C Y C L O P S 2 ...... STAR data name checker
C      *%%%%%\/|\/%%%%%*
C        *%%%%%%%%%%%*                Version 2.1.5
C           *%%%%%*                  29 November 2009
C
C
C
C CYCLOPS2 is a fortran program for checking STAR data names against data
C -------- name dictionaries written in DDL-STAR format proposed by Tony Cook
C          of ORAC Ltd., Leeds. Data names may be checked in any text file.
C
C           CYCLOPS Version 2
C                  by
C
C                  Sydney R. Hall (syd@crystal.uwa.edu.au)
C                  Crystallography Centre
C                  University of Western Australia
C                  Nedlands 6009, AUSTRALIA
C                  
C                  and
C                  
C                  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
C                  Bernstein + Sons
C                  5 Brewster Lane
C                  Bellport, NY 11713, U.S.A.
C
C          Version 2 handles DDL1 and DDL2 dictionaries, and uses CIFtbx2
C
C The latest program source and information is available from:
C
C Em: syd@crystal.uwa.edu.au       ,-_|\      Sydney R. Hall
C sendcif@crystal.uwa.edu.au      /     \     Crystallography Centre
C Fx: +61 9 380 1118          --> *_,-._/     University of Western Australia
C Ph: +61 9 380 2725                   v      Nedlands 6009, AUSTRALIA
C
C
C CYCLOPS reads the text containing data names from the standard input device
C (normally device 5). One or more STAR data name dictionaries (in DDL format)
C are opened (normally from unit 1).  A report on the data names is output to
C the standard output device (normally device 6).  Messages are output to the
C standard error device (normally device 0).  If an error is encountered in
C attempting to write to the standard error device, further messages output
C is diverted to the standard output device.
C
C NOTE:  An attempt is made to open a file named "STARDICT" to use
C either as a dictionary or from which to obtain a list of names
C of dictionaries.  If any dictionaries are specified on the command
C line or by $CYCLOPS_CHECK_DICTIONARY, STARDICT is not examined.
C
C STARDICT may directly contain a dictionary as in version 1 of CYCLOPS, or
C contain a list of file names of dictionaries, one per line.  In that case
C the list must begin with the line for which the first 5 characters are
C "#DICT".  The list of unreferenced data names is suppressed unless
C the line "#VERBOSE" appears in STARDICT, or the command line argument
C "-v yes" is used.  Both the list of names used and the list of unreferenced
C data are suppressed if the line #SHORT appears in STARDICT, or the command
C line argument "-s yes" is used.
C
C
C In a unix-like environment, the program is run as:
C
C cyclops [-i input_text] [-o validation_output] \
C         [-d dictionary] [-p priority] \
C         [-f command_file] [-v verbose] [-s short]\
C         [[[input_text][[validation_output] [dictionary]]]
C
C where:
C         input_text defaults to $CYCLOPS_INPUT_TEXT or stdin
C         validation_output defaults to $CYCLOPS_VALIDATION_OUT or stdout
C         dictionary defaults to $CYCLOPS_CHECK_DICTIONARY
C           (multiple dictionaries may be specified)
C         input_text of "-" is stdin, validation_output of "-" is stdout
C         -p has values of first, final or nodup 
C           (default first for first dictionary has priority) ',
C         -v has values of t or 1 or y vs. f or 0 or n (default f)
C         -s has values of t or 1 or y vs. f or 0 or n (default f)
C
C
C .............................Installation notes.............................
C
C
C CYCLOPS2 is provided as a part of the CIFtbx2 release kit, as the
C file 'cyclops.f', 'cyclops.cmn' and README.cyclops.  CIFtbx2 is
C needed for the files ciftbx.f, ciftbx.sys, ciftbx.cmn and
C hash_funcs.f.  See the CIFtbx2 documentation for further information.
C For non-unix-like environments, you will have to provide replacements
C for iargc, getarg and getenv.  The following are reasonable 
C possibilities:
C
C         integer function iargc(dummy)
C         iargc=1
C         return
C         end
C
C         subroutine getarg(narg,string)
C         integer narg
C         character*(*) string
C         string=char(0)
C         if (narg.eq.1) string='-vn'
C         return
C         end
C
C         subroutine getenv(evar,string)
C         character*(*) evar,string
C         string=char(0)
C         if(evar.eq."CYCLOPS_INPUT_TEXT")
C        *  string='STARTEXT'//char(0)
C         if(evar.eq."CYCLOPS_VALIDATION_OUT")
C        *  string='STARCHEK'//char(0)
C         if(evar.eq."CYCLOPS_CHECK_DICTIONARY") 
C        *  string=char(0)
C         return
C         end
C
C This combination of substitute routines would "wire-in" CYCLOPS2 to
C read its input text from a file named STARTEXT, write the
C validations output to a file named STARCHEK, and check names 
C against the default dictionary or file listing dictionaries
C STARDICT
C
C Note that the PARAMETER 'MAXBUF' should contain the maximum number 
C of characters contained on a single text line. The default value is 200.
C
C
C ..................................Files Used................................
C
C  dictionary input         input   on device 1
C    Optionally a list of dictionaries, 1 per line after a #DICT line
C  validation output        output  on device 6 ('stdout')
C  Input text               input   on device 1, if a file, 5  if 'stdin'
C  Message device           output  on device 0 ('stderr')
C    may be diverted to the validation output file
C  Direct access            in/out  on device 3
C
C***************************************************************************
C
C **** PROGRAM CYCLOPS
C
      include 'ciftbx.sys'
      include 'cyclops.cmn'
C
      iundac = 3
      iunerr = 0
      iunout = 6
      iuninp = 1
      verbose=.false.
      short  =.false.
      CASE=ICHAR('a')-ICHAR('A')
      XTAB=CHAR(05)
      IF(CASE.LT.0) GOTO 100
      XTAB=CHAR(09)
C
      VERS='2.1.5 (29 Nov 09)'
      WRITE(IUNERR,
     *  '(/35H <CYCLOPS2> STAR data name checker ,A,/1X,59(1H-))'
     *  ,ERR=40) VERS
      GO TO 50
 40   IUNERR=IUNOUT
C
 50   call getfiles
      if (valout(1:lvalout) .ne.' ' .and. 
     *  valout(1:lvalout) .ne. '-') then
        OPEN(IUNOUT,
     *  FILE=valout(1:lvalout),STATUS='UNKNOWN',FORM='FORMATTED')
      endif
      if (iunerr.eq.iunout)
     *  WRITE(IUNERR,
     *  '(/35H <CYCLOPS2> STAR data name checker ,A,/1X,59(1H-))')
     *  VERS

      iuninp=5
      if (inptext(1:linptext) .ne.' ' .and. 
     *  inptext(1:linptext) .ne. '-') then
        iuninp=1
        OPEN(IUNINP,
     *  FILE=inptext(1:linptext),STATUS='OLD',FORM='FORMATTED')
      endif
C
C
100   CONTINUE
      CALL CHEKDN
      CALL REPORT
C
      END
C
C***************************************************************************
C     CHEKDN checks the input text file on the standard input file for data
C     names and checks them against the stored dictionary names in dicnam().
C     If it is not present the data name is added to the end of dicnam(). The
C     line numbers are added to NPOS(). XGETSTR supplies strings as the
C     variable 'stri' of length 'nchr'. NTYP=1 for a data name.
C***************************************************************************
C
      SUBROUTINE CHEKDN
      integer  lastnb
      include 'ciftbx.sys'
      include 'cyclops.cmn'
      INTEGER I,J,N,IL,IR
C
      CHARACTER HEAD*14,TAIL*19,C*1,CBS*1
      CHARACTER*20 XDICNAM(30)
      INTEGER      NAMLEN(30),NNCHR
      DATA HEAD/'  ,.:([{</|"  '/
      DATA TAIL/'  ,.?)]}>!;:/-=|"  '/
      DATA XDICNAM /
     * '_category           ','_definition         ',
     * '_dictionary_history ','_dictionary_name    ',
     * '_dictionary_update  ','_dictionary_version ',
     * '_enumeration        ','_enumeration_default',
     * '_enumeration_detail ','_enumeration_range  ',
     * '_example            ','_example_detail     ',
     * '_list               ','_list_level         ',
     * '_list_link_child    ','_list_link_parent   ',
     * '_list_mandatory     ','_list_reference     ',
     * '_list_uniqueness    ','_name               ',
     * '_related_item       ','_related_function   ',     
     * '_type               ','_type_conditions    ',
     * '_type_construct     ','_units              ',
     * '_units_detail       ','                    ',
     * '                    ','                    '/
C
      DO I = 1,30
        NAMLEN(I)=nblen(XDICNAM(I))
      ENDDO
C
      CBS='\\'
      HEAD(1:1)=''''
      HEAD(2:2)=CBS(1:1)
      HEAD(14:14)=XTAB
      TAIL(1:1)=''''
      TAIL(2:2)=CBS(1:1)
      TAIL(19:19)=XTAB
      MCNT=NDICT+1
      LLINE=0
      IEOF=0
      MAXC=MAXBUF
      ICHR=MAXC
      STRI=' '
C
10    CALL XGETSTR
      IF(IEOF.NE.0)  GOTO 200
      I=1
      J=NCHR
      IF(NTYP.EQ.1)  GOTO 20
      IF(NTYP.NE.4)  GOTO 10
C -------------------------------------- look for data name in string
      IL=3
      DO 15 I=1,NCHR
      IF(STRI(I:I).EQ.'_') GOTO 20
15    CONTINUE
      GOTO 10
C
20    IF(I.EQ.1) GOTO 22
      IL=INDEX(HEAD,STRI(I-1:I-1))
      IF(IL.GE.6.AND.IL.LE.8) THEN
        ILEV=1
        DO J = I,NCHR
          IF(STRI(J:J).EQ.STRI(I-1:I-1)) ILEV=ILEV+1
          IR=INDEX(TAIL,STRI(J:J))
          IF (IR.EQ.IL) THEN
            ILEV=ILEV-1
            IF(ILEV.EQ.0) THEN
              ICHR=ICHR-(NCHR-J)+1
              NCHR=J
              GOTO 21
            ENDIF
          ENDIF
        ENDDO
      ENDIF
21    IF(IL.EQ.0 .AND. STRI(I-1:I-1).NE.'*') GOTO 10
      NNCHR=INDEX(STRI(I:NCHR),' ')
      IF(NNCHR.GT.0 .AND. NNCHR.LT.NCHR-I+1) THEN
        NCHR=NNCHR+I-1
      ENDIF
22    DO 25 J=NCHR,I,-1
      IR=INDEX(TAIL,STRI(J:J))
      IF(IR.EQ.0. .OR.
     *  (IR.EQ.7.AND.IR.NE.IL)) GO TO 30
25    CONTINUE
      J=I
30    IF(J.GT.I+1.AND.STRI(J-1:J).eq.CBS(1:1)//'n') THEN
       NCHR=J-2
       GOTO 22
      ENDIF
      NCHR=J-I+1
C -------------------------------------- identify mode of search
35    IF(NCHR.LT.2) GOTO 10
      IF(NCHR.GT.NUMCHAR)
     * CALL CERR(0,'Data name in text is > NUMCHAR chars',
     *        STRI(I:J))
      NCHR = MIN(NUMCHAR,NCHR)
      XNAME=STRI(I:J)
      UNAME=XNAME
      DO 37 N=1,NCHR
      C=XNAME(N:N)
      IF(C.GE.'A'.AND.C.LE.'Z') XNAME(N:N)=CHAR(ICHAR(C)+CASE)
37    CONTINUE
      IF(STRI(J:J).EQ.'_')     GOTO 50
      IF(I.EQ.1)               GOTO 40
      IF(STRI(I-1:I-1).EQ.'*') GOTO 60
C --------------------------------------- search for full data name
40    call hash_find(XNAME,dicnam,dicchain,NUMDICT,ndict,
     *  dichash,NUMHASH,I)
      IF (I.NE.0) GOTO 170
      GOTO 100
C ---------------------------------------- search for group data name
50    DO 55 I=1,NDICT
      IF(XNAME(1:NCHR).EQ.dicnam(I)(1:NCHR)) GOTO 170
55    CONTINUE
      GOTO 100
C ---------------------------------------- search for end-part data name
60    DO 65 I=1,NDICT
      N=lastnb(dicnam(I))
      IF(NCHR.GT.N) GOTO 65
      IF(XNAME(1:NCHR).EQ.dicnam(I)(N-NCHR+1:N)) GOTO 170
65    CONTINUE
C ---------------------------------------- new data name add to list
100   DO 120 I=1,30
      IF(NAMLEN(I).NE.NCHR)       GOTO 120
      IF(XDICNAM(I).EQ.XNAME(1:20)) GOTO 10
120   CONTINUE
      J = NDICT+1
      call hash_store(XNAME,dicnam,dicchain,NUMDICT,ndict,
     *  dichash,NUMHASH,I)
      IF (I.EQ.0) 
     *  CALL CERR(1,'Data name table exceeded (Current max is NUMDICT)'
     *  ,' ')
      IF (J.NE.I) GOTO 170
      dictag(I) = UNAME
      DO 160 J=1,20
160   NPOS(NDICT,J)=0
      I=NDICT
C --------------------------------------- add line numbers to npos()
170   IF(NPOS(I,1).EQ.19) THEN
        DO J = 12,19
        NPOS(I,J)=NPOS(I,J+1)
        ENDDO
        NPOS(I,12)=-NPOS(I,12)
        NPOS(I,1)=NPOS(I,1)-1
      ENDIF
      NPOS(I,1)=NPOS(I,1)+1
      N=NPOS(I,1)
180   NPOS(I,N+1)=LLINE
      GO TO 10
C
200   RETURN
      END
C
C***************************************************************************
C     REPORT prints a summary of the line occurrences for each data name
C     in the dictionary and the text file itself (if it does not appear
C     in the dictionary). This summary is output as the file STARCHEK.
C***************************************************************************
C
      SUBROUTINE REPORT
      integer  lastnb
      include 'ciftbx.sys'
      include 'cyclops.cmn'
      INTEGER I,II,J,JJ,KK,MAXPAS,N,NNT,NT,NPASS
      CHARACTER*2 BREAK
      CHARACTER*4 RFLAG
      CHARACTER*4 DNUM
      INTEGER NUMCP7
      PARAMETER (NUMCP7 = NUMCHAR+7)
      CHARACTER*(NUMCP7) NAMT
C
      WRITE(IUNOUT,
     *  '(//,20X,18HCYCLOPS Check List,/,20x,18(1H-),//)')
      WRITE(IUNOUT,
     *  '(15X,25H Dictionary data names  =,I5)') MCNT-1
      WRITE(IUNOUT,
     *  '(15X,25H New data names in text =,I5)') NDICT-MCNT+1
      WRITE(IUNERR,
     *  '(25H New data names in text =,I5)') NDICT-MCNT+1
      NT=0
      DO I = 1,KDICT
        IF (I.LT.10) THEN
          WRITE(DNUM,'(''['',I1,''] '')') I
        ELSE
          WRITE(DNUM,'(''['',I2,'']'')') I
        ENDIF
        WRITE(IUNERR,
     *  '(1x,A,I5)')DNUM//' Dictionary '//
     *  DICTNM(I)(1:MAX(1,LASTNB(DICTNM(I))))//
     *  ' '//
     *  DICTVR(I)(1:MAX(1,LASTNB(DICTVR(I))))//
     *  ' data names = ',
     *  MDICT(I)-NT
        IF(KDICT.GT.1) THEN
        WRITE(IUNOUT,'(16X,A,I5)')DNUM//' Dictionary '//
     *  DICTNM(I)(1:MAX(1,LASTNB(DICTNM(I))))//
     *  ' '//
     *  DICTVR(I)(1:MAX(1,LASTNB(DICTVR(I))))//
     *  ' data names = ',
     *  MDICT(I)-NT
        ELSE
        WRITE(IUNOUT,'(16X,A,I5)')'Dictionary '//
     *  DICTNM(I)(1:MAX(1,LASTNB(DICTNM(I))))//
     *  ' '//
     *  DICTVR(I)(1:MAX(1,LASTNB(DICTVR(I))))//
     *  ' data names = ',
     *  MDICT(I)-NT
        ENDIF
        NT=MDICT(I)
      ENDDO
      IF (NDICT.GT.MDICT(KDICT)) THEN
        WRITE(IUNOUT,'(//1X,28HData names NOT in Dictionary,25X,
     *  12HLine Numbers,/)')
        call ssort(DICNAM(MDICT(KDICT)+1),
     *  NDICT-MDICT(KDICT),NUMCHAR,NORDER)
        DO KK=1,-MDICT(KDICT)+NDICT
        I=MDICT(KDICT)+NORDER(KK)
        N=NPOS(I,1)
        BREAK="  "
        NAMT = DICTAG(I)
        NNT = LASTNB(NAMT)
        DO II = 1,NUMCP7,2
          IF(II.GE.NNT+3) NAMT(II:II)='.'
        ENDDO
        IF(NPOS(I,12).LT.0) BREAK=".."
        IF(N.LE.10) THEN
        IF(N.LE.4) THEN
        WRITE(IUNOUT,'(1X,A,4I6)')
     *  NAMT,
     *  (IABS(NPOS(I,J+1)),J=1,N)
        ELSE
        WRITE(IUNOUT,'(1X,A,4I6,/,44X,6I6)')
     *  NAMT,
     *  (IABS(NPOS(I,J+1)),J=1,N)
        ENDIF
        ELSE
        IF(N.LE.16) THEN
        WRITE(IUNOUT,'(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6)')
     *  NAMT,(IABS(NPOS(I,J+1)),J=1,10),BREAK,
     *  (IABS(NPOS(I,J+1)),J=11,N)
        ELSE
        WRITE(IUNOUT,
     *  '(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6,/,44X,6I6)')
     *  NAMT,(IABS(NPOS(I,J+1)),J=1,10),BREAK,
     *  (IABS(NPOS(I,J+1)),J=11,N)
        ENDIF
        ENDIF
        ENDDO
      ENDIF
C
      if (short) go to 900
      NT=1
      DO I=1,MDICT(KDICT)
      IF(I.EQ.MDICT(NT)+1) NT=NT+1
      NNDICT(I)=NT
      ENDDO
      call ssort(DICNAM,MDICT(KDICT),NUMCHAR,NORDER)
      MAXPAS = 1
      IF (verbose) MAXPAS = 2
      DO NPASS = 1,MAXPAS
      WRITE(IUNOUT,'(//)')
      DO II = 1,KDICT
      IF (II.LT.10) THEN
        WRITE(DNUM,'(''['',I1,''] '')') II
      ELSE
        WRITE(DNUM,'(''['',I2,'']'')') II
      ENDIF
      IF(KDICT.GT.1) THEN
      WRITE(IUNOUT,'(1x,A)')DNUM//' Dictionary '//
     *  DICTNM(II)(1:MAX(1,LASTNB(DICTNM(II))))//
     *  ' '//
     *  DICTVR(II)(1:MAX(1,LASTNB(DICTVR(II))))
      ELSE
      WRITE(IUNOUT,'(1x,A)')'Dictionary '//
     *  DICTNM(II)(1:MAX(1,LASTNB(DICTNM(II))))//
     *  ' '//
     *  DICTVR(II)(1:MAX(1,LASTNB(DICTVR(II))))
      ENDIF
      ENDDO
      IF (NPASS.EQ.1) THEN
          WRITE(IUNOUT,'(54X,12HLine Numbers,/)') 
        ELSE
          WRITE(IUNOUT,'(54X,20HNames Not Referenced,/)')
        ENDIF
C
        NT=1
        DO 50 JJ=1,MDICT(KDICT)
        I=NORDER(JJ)
        IF(NBLEN(DICNAM(I)).NE.0) THEN
20      N=NPOS(I,1)
        BREAK="  "
        IF(NPOS(I,12).LT.0) BREAK=".."
        IF((NPASS.EQ.1.AND.N.GT.0).OR.
     *   (NPASS.EQ.2.AND.N.EQ.0)) THEN
C       IF (NPASS.EQ.2.AND.AROOT(I).NE.0.AND.AROOT(I).NE.I) GO TO 50
        IF (NPASS.EQ.2.AND.ALIAS(I).NE.0) THEN
          II = I
30        IF(NPOS(II,1).NE.0) GO TO 50
          II = ALIAS(II)
          IF(II.NE.0) GO TO 30
        ENDIF
        NAMT = DICTAG(I)
        IF (KDICT.GT.1) THEN
        IF (NNDICT(I).LT.10) THEN
          WRITE(RFLAG,'(''['',I1,''] '')') NNDICT(I)
        ELSE
          WRITE(RFLAG,'(''['',I2,'']'')') NNDICT(I)
        ENDIF
        NAMT = RFLAG//DICTAG(I)(1:MAX(1,LASTNB(DICTAG(I))))
        ENDIF
        NNT = LASTNB(NAMT)
        IF (NPASS.EQ.1) THEN
        DO II = 1,NUMCP7,2
        IF(II.GE.NNT+3) NAMT(II:II)='.'
        ENDDO
        ENDIF
        IF(N.LE.10) THEN
          IF(N.LE.4) THEN
          WRITE(IUNOUT,'(1X,A,4I6)')
     *    NAMT,
     *    (IABS(NPOS(I,J+1)),J=1,N)
          ELSE
          WRITE(IUNOUT,'(1X,A,4I6,/,44X,6I6)')
     *    NAMT,
     *    (IABS(NPOS(I,J+1)),J=1,N)
          ENDIF
        ELSE
          IF(N.LE.16) THEN
          WRITE(IUNOUT,'(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6)')
     *    NAMT,(IABS(NPOS(I,J+1)),J=1,10),BREAK,
     *    (IABS(NPOS(I,J+1)),J=11,N)
          ELSE
          WRITE(IUNOUT,
     *    '(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6,/,44X,6I6)')
     *    NAMT,(IABS(NPOS(I,J+1)),J=1,10),BREAK,
     *    (IABS(NPOS(I,J+1)),J=11,N)
          ENDIF
        ENDIF
        IF (AROOT(I).NE.0.OR.ALIAS(I).NE.0) THEN
          II = AROOT(I)
          IF (II.EQ.0) II = I
40        IF (II.NE.I) THEN
            RFLAG=' '
            IF(KDICT.GT.1.AND.NNDICT(II).NE.NNDICT(I)) THEN
            IF (NNDICT(II).LT.10) THEN
              WRITE(RFLAG,'(''['',I1,''] '')') NNDICT(II)
            ELSE
              WRITE(RFLAG,'(''['',I2,'']'')') NNDICT(II)
            ENDIF
            ENDIF
            NAMT = RFLAG//'= '//
     *        DICTAG(II)(1:max(1,lastnb(DICTAG(II))))
            IF(NPOS(II,1).eq.0) THEN
            WRITE(IUNOUT,'(1X,A)') NAMT
            ELSE
            BREAK="  "
            IF(NPOS(II,12).LT.0) BREAK=".."
            IF(NPOS(II,1).LE.10) THEN
              IF(NPOS(II,1).LE.4) THEN
                WRITE(IUNOUT,'(1X,A,4I6)')
     *            NAMT,
     *            (IABS(NPOS(II,J+1)),J=1,NPOS(II,1))
              ELSE
                WRITE(IUNOUT,'(1X,A,4I6,/,44X,6I6)')
     *          NAMT,
     *          (IABS(NPOS(II,J+1)),J=1,NPOS(II,1))
              ENDIF
            ELSE
            IF(NPOS(II,1).LE.16) THEN
            WRITE(IUNOUT,
     *        '(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6)')
     *        NAMT,(IABS(NPOS(II,J+1)),J=1,10),BREAK,
     *        (IABS(NPOS(II,J+1)),J=11,NPOS(II,1))
            ELSE
            WRITE(IUNOUT,
     *        '(1X,A,4I6,/,44X,6I6,/41X,A2,1X,6I6,/,44X,6I6)')
     *        NAMT,(IABS(NPOS(II,J+1)),J=1,10),BREAK,
     *        (IABS(NPOS(II,J+1)),J=11,NPOS(II,1))
            ENDIF
            ENDIF
            ENDIF
          ENDIF
          II=ALIAS(II)
          IF (II.NE.0) GOTO 40
        ENDIF
        ENDIF
        ENDIF
50      CONTINUE
      ENDDO
C
900   CLOSE(7)
      RETURN
      END
C
C***************************************************************************
C     XGETSTR extract a string of characters delimited by blanks, XTABs, quotes.
C     NTYP=1 data name; =3 number data; =4 character data; =5 text data.
C***************************************************************************
      SUBROUTINE XGETSTR
      include 'ciftbx.sys'
      include 'cyclops.cmn'
C
      INTEGER       I
      CHARACTER C*1,NUM*13,CBS*2
      DATA NUM/'0123456789+-.'/
      maxc=maxbuf
      CBS='\\'
      CBS(2:2)='n'
C
C---------------------------------- loop over data items in each line
10    ICHR=ICHR+1
      IF(ICHR.LE.MAXC)     GOTO 30
C---------------------------------- get a new line from device iuninp
20    READ(IUNINP,'(A)',END=95) IBUF(1:MAXC)
      LLINE=LLINE+1
      ICHR=1
      NTYP=0
C---------------------------------- test for text data
      IF(IBUF(1:1).NE.';') GOTO 30
      ICHR=2
C---------------------------------- test for delimiter character
30    C=IBUF(ICHR:ICHR)
      IF(C.EQ.' ')       GOTO 10
      IF(C.EQ.XTAB)      GOTO 10
      IF(ICHR.LT.MAXC) THEN
        IF(IBUF(ICHR:ICHR+1).EQ.CBS) THEN
          IBUF(ICHR:ICHR+1) = '  '
          GOTO 10
        ENDIF
      ENDIF
      IF(C.EQ.'#')       GOTO 20
      IF(C.EQ.'''')      GOTO 60
      IF(C.EQ.'"')       GOTO 60
C---------------------------------- test for data name
      IF(C.NE.'_')       GOTO 40
      NTYP=1
      GOTO 45
C---------------------------------- test for number or character data
40    NTYP=3
      IF(INDEX(NUM,C).EQ.0) NTYP=4
C---------------------------------- get blank-limited char string
45    DO 50 I=ICHR,MAXC
      IF(IBUF(I:I).EQ.' ')       GOTO 90
      IF(IBUF(I:I).EQ.XTAB)       GOTO 90
      IF(I.LT.MAXC) THEN
        IF(IBUF(I:I+1).EQ.CBS) THEN
          IBUF(I:I+1) = '  '
          GOTO 90
        ENDIF
      ENDIF
50    CONTINUE
      I=MAXC+1
      GOTO 90
C---------------------------------- get quote-limited char string
60    NTYP=4
      DO 70 I=ICHR+1,MAXC
      IF(IBUF(I:I).EQ.C)         GOTO 80
70    CONTINUE
      I=MAXC
      IBUF(I:I)=C
80    I=I+1
C---------------------------------- count & store string
90    NCHR=I-ICHR
      ICHR=I-1
      STRI(1:NCHR)=IBUF(I-NCHR:ICHR)
      GOTO 100
95    IEOF=1
100   RETURN
      END
C
C***************************************************************************
C     Error message generator. Output to device 6 and stop.
C***************************************************************************
      SUBROUTINE CERR(I,MESS,STRG)
      include 'ciftbx.sys'
      include 'cyclops.cmn'
      INTEGER I
C
      CHARACTER*(*) MESS,STRG
C
      WRITE(iunerr,'(11H Error >>> ,A,2X,A)') MESS,STRG
      IF(I.EQ.0) RETURN
C
      WRITE(iunerr,'(11X,30HFatal error -- Text file line ,I5/)') LLINE
      CLOSE(4)
50    STOP
      END
       subroutine ssort(strg,nstrg,nchr,norder)
C
C      simple sort routine to establish order of array
C      strg of nstrg string of length nchr in array norder
C      no attempt is made to remove duplicates nor to
C      preserve the initial ordering of duplicates.
C
C      Herbert J. Bernstein, Bernstein+Sons, 13 July 1996
C
       integer nstrg, nchr, ns, ii, jj, ilow, ihigh, imid
       integer norder(nstrg)
       character*(*) strg(nstrg)
       norder(1)=1
       ns=1
       do ii = 2,nstrg
       ilow = 1
       ihigh = ns
 100   imid = (ilow+ihigh+1)/2
       if (ihigh-ilow.le.1) then
       if (strg(ii) .ge. strg(norder(ilow)).and.
     *     strg(ii) .le. strg(norder(ihigh))) then
         do jj = ns,ihigh,-1
           norder(jj+1) = norder(jj)
         enddo
         norder(ihigh) = ii
         ns=ns+1
         go to 200
       endif
       endif
       if (strg(ii) .ge. strg(norder(ihigh))) then
         do jj = ns,ihigh+1,-1
           norder(jj+1) = norder(jj)
         enddo
         norder(ihigh+1) = ii
         ns=ns+1
       else
         if (strg(ii) .le. strg(norder(ilow))) then
           do jj = ns,ilow,-1
             norder(jj+1) = norder(jj)
           enddo
           norder(ilow) = ii
           ns=ns+1
         else
           if (strg(ii) .le. strg(norder(imid))) then
             ihigh = imid
           else
             ilow = imid
           endif
           go to 100
         endif
       endif
 200   continue
       enddo
       return
       end
C
C
      subroutine getfiles
C
C     This routine, derived from the similar routine in cif2pdb
C     by H. J. Bernstein and F. C. Bernstein, is used to control
C     initial processing of command line arguments and envrionment
C     variables
C
C     On systems lacking a unix style iargc and getarg, provide
C     equivalent routines.  At the very least, you will need
C     to have iargc() return 0 (for no command line arguments)
C     and have getarg(karg,strg) return a null-terminated
C     empty string in strg.
C
C     On systems lacking a unix-style getenv, provide a dummy
C     routine returning a null-terminated empty string
C
C     
C
      integer  lastnb
      logical  init_
      logical  dict_
      include 'ciftbx.sys'
      include 'cyclops.cmn'
      character*256 temp,temp2,cline
      character*256 cldict(99)
      integer ncld
      integer nblen,i,j,kk
      character*(MAXBUF) XTEST
      character*5 prior, cats
      integer iargc,ll,karg,kfarg,nfarg,
     *  ifound,iwant,isi,iso,ii
      logical backarg,ffile
      ncld = 0
      numarg = iargc()
      call getenv("CYCLOPS_INPUT_TEXT",inptext)
      call getenv("CYCLOPS_VALIDATION_OUT",valout)
      call getenv("CYCLOPS_CHECK_DICTIONARY",ckdict)
      lckdict = max(1,nblen(ckdict))
      if(ckdict(1:lckdict).ne.' ') then
        ncld = ncld+1
        if(ncld.gt.99)
     *    call cerr(1,'Too many dictionaries',' ')
        cldict(ncld) = ckdict(1:lckdict)
        ckdict = " "
        lckdict = 1
      endif
      IUNINP=1
      IEOF=0
      NDICT=0
      FLAG=0
      LLINE=0
      KDICT=1
      karg = 0
      iwant = 0
      ifound = 0
      isi = 0
      iso = 0       
      backarg = .false.
      ffile = .false.
      prior = ' '
      cats = 'catno'
 100  if(.not.ffile) then
        karg = karg+1
        if (karg.le.numarg) then
          if (.not.backarg) then
            call getarg(karg,temp)
            temp2 = ' '
          else
            temp = temp2
            backarg = .false.
          endif
        else
          go to 500
        endif
      else
        go to 300
 200    close(unit=iuninp)
        ffile = .false.
        go to 100
 300    kfarg = kfarg+1
        if (kfarg.gt.nfarg) then
          read(iuninp,'(a)',end=200) cline
          call splitstr(nfarg,cline,cstr,128,' ')
          kfarg = 0
          go to 300
        else
          if (.not.backarg) then
            temp = cstr(kfarg)
            temp2 = ' '
          else
            temp = temp2
            backarg = .false.
          endif
        endif
      endif
      ll = nblen(temp)
      if (ll.eq.0) then
        temp = ' '
        ll = 1
      endif
      if (iwant.ne.0) then
        if (iwant.eq.1) inptext = temp(1:ll)
        if (iwant.eq.2) valout = temp(1:ll)
        if (iwant.eq.3) then
          ckdict = temp(1:ll)
          lckdict = max(1,nblen(ckdict))
          ncld = ncld+1
          if(ncld.gt.99)
     *    call cerr(1,'Too many dictionaries',' ')
          cldict(ncld) = ckdict(1:lckdict)
          ckdict = " "
          lckdict = 1
        endif
        if (iwant.eq.4) then
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            verbose=.false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              verbose=.true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.5) then
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            short=.false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              short=.true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.6) then
          if(prior.ne.' ') go to 900
          if(temp(1:ll).eq.'NODUP' .or.
     *      temp(1:ll).eq.'nodup') then
            prior='nodup'
          else
            if(temp(1:ll).eq.'FIRST' .or.
     *        temp(1:ll).eq.'first') then
              prior='first'
            else
              if(temp(1:ll).eq.'FINAL' .or.
     *          temp(1:ll).eq.'final') then
                prior='first'
              else
                go to 900
              endif
            endif
          endif
        endif
        if (iwant.eq.7) then
          open(unit=iuninp,file=temp(1:ll),status='OLD',err=900)
          kfarg = 1
          nfarg = 0
          ffile = .true.
        endif
        if (iwant.eq.8) then
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            cats='catno'
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              cats='catck'
            else
              go to 900
            endif
          endif
        endif
        iwant = 0
      else
        if (temp(1:1).eq."-") then
          temp2=temp(3:256)
          if (ll.gt.2) then
            backarg = .true.
            karg = karg-1
          endif
          if (temp(2:2).eq."i") then
            iwant = 1
            if (isi.gt.0) go to 900
            isi = 1
          endif
          if (temp(2:2).eq."o") then
            iwant = 2
            if (iso.gt.0) go to 900
            iso = 1
          endif
          if (temp(2:2).eq."d") iwant = 3
          if (temp(2:2).eq."v") iwant = 4
          if (temp(2:2).eq."s") iwant = 5
          if (temp(2:2).eq."p") iwant = 6
          if (temp(2:2).eq."f") iwant = 7
          if (temp(2:2).eq."c") iwant = 8
          if (iwant.eq.0) go to 900
        else
          ifound = ifound+1
          if (ifound.eq.1) then
            if (isi.gt.0) then
              ifound = ifound+1
            else
              inptext = temp(1:ll)
              isi = 1
            endif
          endif
          if (ifound.eq.2) then
            if (iso.gt.0) then
              ifound = ifound+1
            else
              valout = temp(1:ll)
              iso = 1
            endif
          endif
          if (ifound.eq.3) then
            ckdict = temp(1:ll)
          endif
          if (ifound.gt.3) go to 900
        endif
      endif
      go to 100
 500  linptext = max(1,nblen(inptext))
      lvalout = max(1,nblen(valout))
      lckdict = max(1,nblen(ckdict))
      if(.not.init_(iuninp,iunout,iundac,iunerr))
     *  call cerr(1,' failed init',' ')
      XDICT(1)='STARDICT'
      kdict = 1
      MDICT(1)=0
      if(ncld.gt.0) goto 560
      open(iuninp,FILE='STARDICT',STATUS='OLD',FORM='FORMATTED',
     *  err=560)
      READ(iuninp,'(A)',END=560)XTEST
      if(XTEST.eq.'#DICT') then
        do kk = 1,100
        READ(iuninp,'(A)',END=550) XDICT(kdict)
        if (XDICT(kdict).eq.'#VERBOSE') then
          verbose=.true.
        else
          if (XDICT(kdict).eq.'#SHORT') then
            short=.true.
          else
            do ii = 1,kdict-1
            if (XDICT(ii).eq.XDICT(kdict)) then
              if(prior.eq.'nodup')
     *          call cerr(1,' Duplicate dictionary ',
     *          XDICT(ii))
              go to 540
            endif
            enddo
            kdict=kdict+1
 540        continue
          endif
        endif
        enddo
 550    kdict=kdict-1
      endif
      go to 570
 560  kdict=0
 570  do ii = 1,ncld
        do kk = 1,kdict
          if (XDICT(kk).eq.cldict(ii)) then
            if(prior.eq.'nodup')
     *        call cerr(1,' Duplicate dictionary ',
     *        XDICT(kk))
            go to 580
          endif
        enddo
        kdict=kdict+1
        if (kdict.gt.99) call cerr(1, 'Too many dictionaries', ' ')
        XDICT(kdict)=cldict(ii)
 580    continue
      enddo
      if(kdict.lt.1) call cerr(1,' Dictionary list empty',' ')
      close(1)
      do i = 1,kdict
      if(.not.dict_(XDICT(i)(1:max(1,lastnb(XDICT(i)))),
     *  prior//cats))
     *  call cerr(1,XDICT(i)(1:max(1,lastnb(XDICT(i))))//
     *' not found',' ')
      MDICT(i)=NDICT
      DICTNM(i)=dicname_
      DICTVR(i)=dicver_
      enddo
      close(unit=iundac)
      DO 610 I=1,NDICT
        DO 660 J=1,20
 660    NPOS(I,J)=0
 610  CONTINUE
      RETURN
 900  write(iunerr,'(a)')
     *  ' cyclops [-i input_text] [-o validation_output] '//
     *           '[-d dictionary] [-p priority] [-c catck]',
     *  '         [-f command_file] [-v verbose] [-s short] '//
     *           '[[[input_text][[validation_output]'//
     *           ' [dictionary]]]',
     *  ' input_text defaults to $CYCLOPS_INPUT_TEXT or stdin',
     *  ' validation_output defaults to $CYCLOPS_VALIDATION_OUT'//
     *  ' or stdout',
     *  ' dictionary defaults to $CYCLOPS_CHECK_DICTIONARY',
     *  ' multiple dictionaries may be specified ',
     *  ' input_text of "-" is stdin ',
     *  ' validation_output of "-" is stdout '
      write(iunerr,'(a)')
     *  ' -p has values of first, final or nodup ',
     *  '   (default first for first dictionary has priority) ',
     *  ' -c has values of t or 1 or y vs. f or 0 or n,',
     *  '   (default f, i.e. no category checking) ',
     *  ' -v has values of t or 1 or y vs. f or 0 or n,',
     *  '   (default f, i.e. non-verbose),',
     *  ' -s has values of t or 1 or y vs. f or 0 or n,',
     *  '   (default f, i.e. not short),',
     *  '   short restricts output to items not in dictionaries'
       stop
       end
C
C
       function nblen(str)
C
C      variant of lastnb which also detects a null character to
C      terminate a string
C
       character*(*) str
       integer nblen,ll,jj,ii
       ll = len(str)
       jj = index(str,char(0))
       if (jj.ne.0) ll = jj-1
       if (ll.eq.0) then
         nblen = 0
         return
       endif
       do ii = 1,ll
       nblen = 1+ll-ii
       if(str(nblen:nblen).ne.' ' .and.
     *   str(nblen:nblen).ne. char(9) ) return
       enddo
       nblen = 0
       return
       end
C
C
C
C
       subroutine splitstr(numf,str,sarry,maxf,fs)
C
C      split the string str into a maximum of maxf fields using the
C      field separator fs.  The number of fields found is reported in
C      numf.  The fields are split into the array sarry, with each
C      field terminated by a null character.
C
C
C      The field separator ' ' is treated as a special case, with
C      full blank and tab removal, and only the count of non-blank
C      fields returned.
C
       integer numf, maxf
       integer jj,jjt,is,ll,lfs
       character*(*) str, sarry(maxf), fs
       ll = len(str)
       lfs = len(fs)
       if (fs.eq.' ') lfs = 1
       is = 1
       numf = 0
 100   jj = index(str(is:ll),fs)
       if (fs.eq.' ') then
         jjt = index(str(is:ll),char(9))
         if (jj.ne.0) then
           if (jjt.ne.0) jj = min(jjt,jj)
         else
           jj = jjt
         endif
       endif
       if (jj.eq.0) then
         numf = numf+1
         if (numf.le.maxf) then
           sarry(numf) = str(is:ll)//char(0)
         endif
         return
       else
         if (fs.ne.' ' .or. jj.gt.1) then
         numf = numf+1
           if (numf.le.maxf) then
             if (jj.gt.1) then
               sarry(numf) = str(is:is-1+jj-1)//char(0)
             else
               sarry(numf) = char(0)
             endif
           endif
         endif
         is = is+jj+lfs-1
         if (is.le.ll) go to 100
         if (fs.ne.' ') then 
           numf = numf+1
           if (numf.le.maxf) sarry(numf) = char(0)
         endif
         return
       endif
       end


