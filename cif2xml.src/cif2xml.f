C
C
C    \ | /      \ | /
C     \|/        \|/ 
C    -- -->>>>>>-- --               c i f 2 x m l  ....CIF to XML COPY PROGRAM
C     /|\        /|\
C    / | \      / | \                       Version 0.1.0 - alpha
C                                              18 November 2000
C
C
C
C 
C cif2xml  is a fortran program using CIFtbx2 to copy a CIF on standard
C -------  input to an equivalent XML on standard output, while checking
C          data names against dictionaries and reformating numbers with
C          esd's to conform to the rule of 19.  A quasar-style request
C          list may be specified, otherwise the entire CIF is copied.
C          The XML output may be literally derived from the CIF input,
C          or transformations may be specified in a dictionary.
C          The declarations required for the XML document may either be
C          embedded in the new document, written to an external DTD, or
C          referred to an existing file.
C
C          This program is based on cif2cif, and differs from the program
C          primarily in the format of the output.
C
C           cif2xml
C                  by
C
C                  Copyright (C) 2000
C                  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
C                  Bernstein + Sons
C                  5 Brewster Lane
C                  Bellport, NY 11713, U.S.A.
C
C This is a copyrighted open source program subject to the license
C conditions in the file NOTICE.html
C
C cif2xml reads the input CIF from the standard input device (normally
C device 5). An optional STAR data name dictionary (in DDL format) is opened.
C A reformatted copy of the input CIF is written to standard output (device 6).
C Messages are output to the standard error device (normally device 0).
C Note that the PARAMETER 'MAXBUF' should contain the maximum number of char-
C acters contained on a single text line. The default value is 200.
C If a request list (a file listing data_ block names and tags) is provided
C that list controls the ordering and selection of tags and values to copy.
C Otherwise the entire CIF is copied in the order presented
C
C In a unix-like environment, the program is run as:
C
C cif2xml [-i input_cif] [-o output_xml] [-d dictionary] [-c catck]\
C         [-f command_file] [-e esdlim_] [-a aliaso_] [-p prefix]\
C         [-t tabl_] [-q request_list] [-b {row|col} ] [-x {xfer|keep|zap}]\
C         [-u {drop|insert}]
C         [-s {inline | referto spec_dtd | writeto spec_dtd} ] \
C         [input_cif [output_xml [dictionary [request_list [spec_dtd ]]]]] 
C
C where:
C         input_cif defaults to $CIF2XML_INPUT_CIF or stdin
C         output_xml defaults to $CIF2XML_OUTPUT_XML or stdout
C         dictionary defaults to $CIF2XML_CHECK_DICTIONARY
C           (multiple dictionaries may be specified)
C         request_list defaults to $CIF2XML_REQUEST_LIST
C         input_cif of "-" is stdin, output_xml of "-" is stdout
C         request_list of "-" is stdin
C         -e has integer values (e.g. 9, 19(default) o 29)
C         -a has values of t or 1 or y vs. f or 0 or n
C         -p has string values in which "_" is replaced by blank
C         -t has values of t or 1 or y vs. f or 0 or n, default f
C         -s defaults to inline, -b defaults to col
C         -x defaults to zap, -u to drop
C
C
C .............................Installation notes.............................
C
C
C cif2xml is provided as a part of the CIFtbx2 release kit, as the
C files 'cif2xml.f' and 'cif2xml.cmn'.  CIFtbx2 is needed for the files 
C ciftbx.f, ciftbx.sys and hash_funcs.f.  See the CIFtbx2 documentation 
C for further information.  For non-unix-like environments, you will have
C to provide replacements for iargc, getarg and getenv.  The following
C are reasonable possibilities:
C
C         integer function iargc(dummy)
C         iargc=0
C         return
C         end
C
C         subroutine getarg(narg,string)
C         integer narg
C         character*(*) string
C         string=char(0)
C         return
C         end
C
C         subroutine getenv(evar,string)
C         character*(*) evar,string
C         string=char(0)
C         if(evar.eq.'CIF2XML_INPUT_CIF')
C        *  string='INPCIF.CIF'//char(0)
C         if(evar.eq.'CIF2XML_OUTPUT_XML')
C        *  string='OUTXML.XML'//char(0)
C         if(evar.eq.'CIF2XML_CHECK_DICTIONARY') 
C        *  string='CIF_CORE.DIC'//char(0)
C         if(evar.eq.'CIF2XML_REQUEST_LIST') 
C        *  string='REQLST.DAT'//char(0)
C         return
C         end
C
C This combination of substitute routines would "wire-in" cif2xml to
C read its input cif from a file named INPCIF.CIF, write its output
C cif to a file named OUTXML.XML, check names against CIF_CORE.DIC
C and use the tag names given in REQLST.DAT to selects the ones to copy. 
C
C ..................................Files Used................................
C
C  dictionary input         input   on device 2
C  Reformatted XML          output  on device 6 ('stdout')
C  Input CIF                input   on device 2, if a file, 5  if 'stdin'
C  Message device           output  on device 0 ('stderr')
C  Direct access            in/out  on device 3
C  Request list             input   on device 4, if a file, 5 if 'stdin'
C
C
C***************************************************************************
C
      program cif2xml
C
      character*28 VERS
      integer  lastnb
      include 'cif2xml.cmn'
      include 'ciftbx.sys'
      include 'ciftbx.cmf'
      character*(MAXBUF) tbxxlocs
      character*(MAXBUF) tbxxxsub
      character*(MAXBUF) string,temp,blok
      character*(MAXBUF) topcat
      integer ipcnt(NUMBLOCK), ipdup(NUMBLOCK)
      character*3 xvc,xtc
      character*1 xxtab,c
      integer ixname,iyname,ii,ipdata,ipsave,kkr,kkrl
      integer kreqnam,lreqnam,kpass,markdb,kbstrt,lbstrt
      integer marklh
      integer nblen,kkp,krr,krrl,jpsave
      double precision dpn,dpesd
C
      VERS='0.1.0 - alpha (18 Nov 2000)'
C
C     Initialization of variables
C
      xxtab = char(5)
      if (ichar('a')-ichar('A') .gt.0) xxtab = char(9)
      iuninp = 2
      iunout = 6
      iundac = 3
      iunerr = 0
      iunrql = 4
      nreqnam = 0
C
      xmlout_ = .true.
      tabl_ = .false.
      align_ = .false.
      tabx_ = .false.
      ptabx_ = .false.
      pdecp_ = .false.
      plzero_ = .false.
      kbstrt = -1
      result = init_(iuninp,iunout,iundac,iunerr)
      call getfiles
      if (tabl_) align_ = .true.
      if (ckdict(1:lckdict).ne.' ')
     *  result = dict_(ckdict(1:lckdict),'valid dtype '//cats)
      if (inpcif(1:linpcif).eq.'-' .or. 
     *  inpcif(1:linpcif) .eq. ' ') then
        iuninp = 5
        inpcif = ' '
        xvc=vcheck
        xtc=tcheck
        result = init_(iuninp,iunout,iundac,iunerr)
        vcheck=xvc
        tcheck=xtc
      endif
C
C     Open input CIF
C
      result = ocif_(inpcif(1:linpcif))
      if(.not.result) then
        call c2xerr(
     *  ' failed to open '//inpcif(1:linpcif))
      endif
      if (nblen(reqlst).gt.0) then
        if (reqlst(1:lreqlst).eq.'-') then
          iunrql = 5
        else
          open(unit=iunrql,file=reqlst(1:lreqlst),status='old',
     *      err=150)
          go to 160
150       call c2xerr(' Unable to open request list '//
     *      reqlst(1:lreqlst))
          go to 190
        endif
160     kkr = MAXBUF+1
        kkrl = 0
        kkp = 1
170     if (kkr.gt.kkrl) then
          read(iunrql,'(a)',end=180) string
          kkr = 1
          kkrl = nblen(string)
          go to 170
        endif
        c = string(kkr:kkr)
        kkr = kkr+1
        if (kkp.eq.1.and.(c.eq.' '.or.c.eq.xxtab)) go to 170
        if (kkp.eq.1.and.c.eq.'#') then
          kkr = kkrl+1
          go to 170
        endif
        if (kkp.gt.1.and.(c.eq.' '.or.c.eq.xxtab)) go to 175
        if (ichar(c).ge.ichar('A')
     *    .and.ichar(c).le.ichar('Z'))
     *    c = char(ichar(c)+ichar('a')-ichar('A'))
        if (kkp.le.NUMCHRR) then
          temp(kkp:kkp) = c
          kkp=kkp+1
        endif
        if (kkr.le.kkrl) go to 170
175     nreqnam = nreqnam+1
        if (nreqnam .le. NUMREQ) then
          reqnam(nreqnam) = temp(1:kkp-1)
          irord(nreqnam) = 0
        else
          call c2xerr(' More than NUMREQ items in request list')
        endif
        kkp = 1
        go to 170
180     if (reqlst(1:lreqlst).ne.'-') close(unit=iunrql)
      endif
190   continue
C
C     open output XML
C
      if (outxml(1:loutxml).eq.'-')  outxml = ' '
      result=pfile_(outxml(1:loutxml))
      if(.not.result) then
        call c2xerr(
     *  ' failed to open '//outxml(1:loutxml))
      endif
      call tbxxpstr(char(0))
      result=pcmnt_(' CIF Copied by cif2xml, version '// VERS)
      result=find_(' ','head',string)
      ipsave = 0
      markdb = 0
      marklh = 0
      kreqnam = 1
C
C     This is the main loop.  return to 199 for a full new pass
C
199   kpass = 1
      lreqnam = kreqnam
      lbstrt = kbstrt
      blok = ' '
      if(markdb.ne.0 .and.nreqnam.ne.0) result=bkmrk_(markdb)
C
200   if (nreqnam.eq.0) then
        call cpcmnt
        result=data_(' ')
        if(.not.result) go to 900
        do ii = 1,nname
          ipcnt(ii) = 1
          irord(ii) = ii
        enddo
      else
        if (kreqnam.gt.nreqnam) go to 900
        if (reqnam(kreqnam)(1:5).eq.'data_') then
          blok = reqnam(kreqnam)(6:NUMCHRR)
          temp = blok
          if (blok.eq.'which_contains:') temp = ' '
          result = data_(temp)
          if (.not.result) then
C
C           We have failed to find the requested block
C
C           If we have made a full circle on a blank name or
C           the block was named, all we can do is report
C           the failure
C
            if (kpass.eq.2.or.temp.ne.' ') then
              call c2xwarn(' Did not find '//
     *          reqnam(kreqnam))
              kreqnam = kreqnam+1
              go to 200
            else
C
C           If we failed on the first pass on a blank name
C           make a second pass from the front of the CIF
C
              kpass = kpass+1
              result = find_(' ','head',string)
              go to 200
            endif
          endif
          kbstrt = recn_
          result=bkmrk_(markdb)
        else
          kreqnam = kreqnam+1
          go to 200
        endif
C
C       Process the list of tags for this block, expanding
C       wild-card requests
C
        do ii = 1,nname
          ipcnt(ii) = 0
        enddo
        krr = 1
C
210     kreqnam = kreqnam+1
        if (kreqnam.gt.nreqnam) go to 220
        if (reqnam(kreqnam)(1:5).eq.'data_') go to 220
        krrl = nblen(reqnam(kreqnam))
        if (reqnam(kreqnam)(krrl:krrl).eq.'_' .or.
     *    reqnam(kreqnam)(krrl:krrl).eq.'.') then
          krrin = krr
          do ii = 1,nname
          if (dname(ii)(1:krrl).eq.
     *      reqnam(kreqnam)(1:krrl) ) then
            ipcnt(ii) = ipcnt(ii)+1
            if (krr.gt.NUMCHRR)
     *        call c2xerr(' Number of tags greater than NUMCHRR')
            irord(krr) = ii
            krr = krr+1
          endif
          enddo
          if (krrin.eq.krr) then
            if (unknown.eq.'insert') then
              if (krr.gt.NUMCHRR)
     *          call c2xerr(' Number of tags greater than NUMCHRR')
              irord(krr) = -kreqnam
              krr = krr+1
            else
              if (blok.ne.'which_contains:'.or.
     *          (kpass.eq.2.and.kbstrt.ge.lbstrt)) then
                 call c2xwarn(' No match to requested name '//
     *             reqnam(kreqnam)(1:krrl))
              endif
            endif
          endif
        else
          call hash_find(reqnam(kreqnam)(1:krrl),
     *     dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *     ii)
          if (ii.gt.0) then
            if (krr.gt.NUMCHRR)
     *        call c2xerr(' Number of tags greater than NUMCHRR')
            irord(krr) = ii
            krr = krr+1
          else
            if (unknown.eq.'insert') then
              if (krr.gt.NUMCHRR)
     *          call c2xerr(' Number of tags greater than NUMCHRR')
              irord(krr) = -kreqnam
              krr = krr+1
            else
              if (blok.ne.'which_contains:'.or.
     *          (kpass.eq.2.and.kbstrt.ge.lbstrt)) then
                call c2xwarn(' No match to requested name '//
     *            reqnam(kreqnam)(1:krrl))
              endif
            endif
          endif
        endif
        go to 210
C
220     krr = krr-1
        if (krr.le.0) then
          if (blok.ne.'which_contains:'.or.
     *      (kpass.eq.2.and.kbstrt.ge.lbstrt)) then
            call c2xwarn(' No tags requested from block '//bloc_)
            go to 199
          else
            kreqnam = lreqnam
            go to 200
          endif
        endif
        if (unknown.eq.'insert') then
          krrfnd = 0
          do ii = 1,krr
            if (irord(ii).gt.0) krrfnd = krrfnd+1
          enddo
          if (krrfnd.eq.0) then
            if (blok.eq.'which contains:' .and.
     *        (kpass.ne.2 .or. kbstrt.lt.lbstrt)) then
              kreqnam = lreqnam
              go to 200
            endif
          endif
        endif
      endif
C
      ipdata=posnam_
      jpsave=ipsave
      ipsave=posval_
      string=bloc_
      saveo_=save_
      globo_=glob_
      if (pdblok(1:1).eq.' ') then
      call tbxxpstr(char(0))
      call tbxxpstr('<!DOCTYPE')
      if (globo_) then
        call putxn('global_')
      else
        if (xmdata.eq.0) then      
          call putxn(string)
        else
          call putxn(tbxxxsub(string,xmlate(xmdata)))
        endif
      endif
      if (usedtd(1:3).ne.'inl') then
        call tbxxpstr('SYSTEM '//'"'//spcdtd(1:lastnb(spcdtd))//'"')
      else
        call tbxxpstr('[')
        call c2xwarn(' inline DTD not supported yet ')
        call tbxxpstr(']')
      endif
      call tbxxpstr('>')
      call tbxxpstr(char(0))
      else
        call c2xwarn(' multiple root elements not legal in XML')
      endif
      result=pdata_(string)
      if(.not.result) then
        call c2xwarn(
     *  ' duplicate data block '//bloc_)
        pposnam_ = 40
        result=pcmnt_
     *    ('<---- duplicate data block')
      endif
      ixname = 1
      if (nreqnam.eq.0) krr = nname
250   if(ixname.le.krr) then
        if (irord(ixname).lt.0) then
          pposval_ = 0
          pposnam_ = 0
          krrl = nblen(reqnam(-irord(ixname)))
          if (unknown.eq.'insert')
     *      result = pchar_(reqnam(-irord(ixname)),'?')
          ixname = ixname+1
          go to 250
        endif
        loop_ = .false.
        if(nloop(irord(ixname)).ne.0) then
        do ii = 1,nname
        ipdup(ii) = 0
        enddo
        ipdup(irord(ixname)) = 1
255       if(jpsave.eq.0 .and.nreqnam.eq.0) call cpcmnt
          jpsave = 0
          if (nreqnam.eq.0) then
             result=find_(' ',' ',string)
            if(type_.ne.'loop') then
            call c2xwarn(' Expected loop_, found: '//string)
            go to 255
            endif
          else
            posval_ = 1
          endif
          pposval_=posval_
          result=ploop_(' ')
          iyname=ixname
          do ii = ixname+1,krr
            if (irord(ii).le.0) go to 260
            if (nloop(irord(ixname)).ne.
     *        nloop(irord(ii))) go to 260
            if (ipdup(irord(ii)).gt.0) go to 260
            iyname=ii
            ipdup(irord(ii)) = 1
          enddo
260       continue
          if (nreqnam .ne.0) result = pcmnt_(' ')
          if (by.eq.'col') then
            result = bkmrk_(marklh)
          else
          do ii = ixname,iyname
            if (nreqnam.eq.0) then
            call cpcmnt
            result=find_(' ','name',string)
            if(tbxxlocs(string).ne.dname(irord(ii))) then
              call c2xwarn(' Expected '
     *         //dname(irord(ii))(1:lastnb(dname(irord(ii))))
     *         //', found '
     *         //string(1:lastnb(string)))
            endif
            else
              result=find_(dname(irord(ii)),'name',string)
            endif
            if(.not.tabl_ .and. by.ne.'col') pposnam_=posnam_
            if (by.ne.'col') then
              result=ploop_(string(1:max(1,lastnb(string))))
            endif
          enddo
          endif
          if (nreqnam.ne.0) then
          result = bkmrk_(markdb)
          result = bkmrk_(markdb)
          endif
270       do ii = ixname,iyname
          if (by.eq.'col') then
            result = bkmrk_(marklh)
            result = bkmrk_(marklh)
            result=find_(dname(irord(ii)),'name',string)
            result = bkmrk_(marklh)
            if(ii.lt.iyname) result = bkmrk_(marklh)
            if(.not.tabl_ .and. by.ne.'col') pposnam_=posnam_
            result=ploop_(string(1:max(1,lastnb(string))))
          endif
275         if (nreqnam.eq.0 .and. by.ne.'col') call cpcmnt
            type_ = 'null'
            result=test_(dname(irord(ii)))
            pposval_ = 0
            pposend_ = 0
            if(.not.tabl_ .and. by.ne.'col') then
              pposval_=posval_
              pposend_=posend_
            endif
            if(type_.eq.'null') then
            result=char_(dname(irord(ii)),string)
            if(.not.tabl_ .and. by.ne.'col') then
              pposval_=posval_
              pposend_=posend_
            endif
              if (long_.eq.1.and.strg_(1:1).eq.'?') then
                 result=pchar_(' ','?')
              else
                 result=pchar_(' ','.')
              endif
              type_=' '
            endif
            if(type_.eq.'char') then
              result=char_(dname(irord(ii)),string)
              pquote_=quote_
              if(.not.tabl_ .and. by.ne.'col') then
                pposval_=posval_
                pposend_=posend_
              endif
              if(string.eq.'?'.or.string.eq.'.') then
                if(pquote_.eq.' ') pquote_=''''
              endif
              result=pchar_(' ',string(1:max(1,lastnb(string))))
              type_=' '
            endif
            if(type_.eq.'text') then
280           result=char_(dname(irord(ii)),string)
              result=ptext_(' ',string(1:max(1,lastnb(string))))
              if(text_) go to 280
              call tbxxeot
              type_=' '
            endif
            if(type_.eq.'numb') then
              dpesd = 0.D0
              result=numd_(dname(irord(ii)),dpn,dpesd)
              if(.not.tabl_ .and. by.ne.'col') then
                pposval_=posval_
                pposend_=posend_
                pposdec_=posdec_
              endif
              pesddig_=esddig_
              pdecp_=decp_
              plzero_=lzero_
              result=pnumd_(' ',dpn,dpesd)
            endif
          if (by.eq.'col') then
            if (loop_) go to 275
            result = pcmnt_(char(0))
          endif
          enddo
          if (nreqnam.ne.0) result=pcmnt_(char(0))
          if(loop_ .and. by.ne.'col') go to 270
          ixname=iyname+1
          if (nreqnam.ne.0) result=pcmnt_(' ')
          go to 250
        else
          if(jpsave.eq.0 .and.nreqnam.eq.0) call cpcmnt
          jpsave=0
          if (nreqnam.eq.0) then
          result=name_(temp)
          else
          result =find_(dname(irord(ixname)),'name',temp)
          endif
          if(.not.result)
     *      call c2xerr(' misaligned data item '//
     *      dname(irord(ixname)))
          if (nreqnam.eq.0) call cpcmnt
          result=test_(dname(irord(ixname)))
          if(type_.eq.'null') then
            if (.not.tabl_) then
              pposval_=posval_
              pposend_=posend_
              pposnam_=posnam_
            endif
            if (long_.eq.1.and.strg_(1:1).eq.'?') then
               result=pchar_(temp,'?')
            else
               result=pchar_(temp,'.')
            endif
            type_=' '
          endif  
          if(type_.eq.'char') then
            result=char_(temp,string)
            pquote_=quote_
            if (.not.tabl_) then
              pposval_=posval_
              pposend_=posend_
              pposnam_=posnam_
            endif
            if(string.eq.'?'.or.string.eq.'.') then
              if(pquote_.eq.' ') pquote_=''''
            endif
            result=pchar_(temp,string(1:max(1,lastnb(string))))
            type_=' '
          endif
          if(type_.eq.'text') then
290         result=char_(temp,string)
            if (.not.tabl_) pposnam_=posnam_
            result=ptext_(temp,string(1:max(1,lastnb(string))))
            if(text_) go to 290
            call tbxxeot
            type_=' '
          endif
          if(type_.eq.'numb') then
            dpesd=0.D0
            result=numd_(temp,dpn,dpesd)
            if (.not.tabl_) then
              pposdec_=posdec_
              pposval_=posval_
              pposend_=posend_
              pposnam_=posnam_
            endif
            pdecp_=decp_
            plzero_=lzero_
            pesddig_=esddig_
            result=pnumd_(temp,dpn,dpesd)
            type_=' '
          endif
        endif
        ixname=ixname+1
        go to 250
      endif
      if (nreqnam.eq.0) call cpcmnt
      pposnam_=ipdata
      pposval_=ipsave
      result=pdata_(' ')
      if(nreqnam.gt.0) then
        result=pcmnt_(' ')
        result=pcmnt_('           ---end-of-data-block--- ')
        result=pcmnt_(' ')
        result=pcmnt_(' ')
      endif
      if(ipsave.gt.0) then
        irecd=lrecd+1
        jrecd=-1
      endif
      jchar=MAXBUF
      go to 199
900   call cpcmnt
      call close_
      stop

C
      END
C
C
C >>>>>> Put out the name at the beginning of a string
C
C        Note that the string may have embedded blanks and
C        parameters.  Only the first token will be used.
C
         subroutine putxn(string)
C
         integer lastnb
         include 'ciftbx.sys'
         character sbuf*(MAXBUF)
         character*(*) string
         integer ik
         logical tbxxxcln

         if (string(1:1).eq.' ') return
         do ik = 1,len(string)
           if (string(ik:ik).eq.' ') go to 100
         enddo
         ik = len(string)+1
 100     ik = ik-1
         sbuf(1:ik)=string(1:ik)
         if (.not.tbxxxcln(sbuf(1:ik),ik)) then
           call c2xwarn(' XML required remapping for '//sbuf(1:ik))
         endif
         call tbxxpstr(sbuf(1:ik))
         return
         end 

C
C
      subroutine cpcmnt
C
C     routine to copy a set of comments (if any) from the
C     input cif to the output cif
C
      include 'cif2xml.cmn'
      logical  cmnt_
      logical  pcmnt_
      include 'ciftbx.cmv'
      character*(MAXBUF) string
100   result=cmnt_(string)
      if (result) then
        pposnam_ = posnam_
        result=pcmnt_(string(1:long_))
        goto 100
      endif
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
      integer  nblen
      logical  dict_
      logical  prefx_
      include 'cif2xml.cmn'
      include 'ciftbx.cmv'
      character*(MAXBUF) tbxxlocs
      character*256 temp,temp2,cline
      character*10 bfill
      character*11 temp3
      integer iargc,ll,karg,kfarg,nfarg,
     *  ifound,iwant,isi,iso,isd,isq,ii
      logical backarg,ffile
      data bfill/'          '/
      numarg = iargc()
      call getenv('CIF2XML_INPUT_CIF',inpcif)
      call getenv('CIF2XML_OUTPUT_XML',outxml)
      call getenv('CIF2XML_CHECK_DICTIONARY',ckdict)
      call getenv('CIF2XML_REQUEST_LIST',reqlst)
      spcdtd = ' '
      usedtd = 'inl'
      cats = 'catno'
      by = 'col'
      unknown = 'drop'
      xferesd = 'zap'
      karg = 0
      iwant = 0
      ifound = 0
      isi = 0
      iso = 0
      isd = 0
      isq = 0
      iss = 0
      backarg = .false.
      ffile = .false.
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
        if (iwant.eq.1) inpcif = temp(1:ll)
        if (iwant.eq.2) outxml = temp(1:ll)
        if (iwant.eq.3) then
          ckdict = temp(1:ll)
          lckdict = max(1,nblen(ckdict))
          result = dict_(ckdict(1:lckdict),'valid dtype '//cats)
          ckdict = ' '
          lckdict = 1
        endif
        if (iwant.eq.4) then
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            aliaso_=.false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              aliaso_=.true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.5) then
          open(unit=iuninp,file=temp(1:ll),status='OLD',err=900)
          kfarg = 1
          nfarg = 0
          ffile = .true.
        endif
        if (iwant.eq.6) then
          if(ll.gt.10) go to 900
          temp3=bfill(1:11-ll)//temp(1:ll)
          read(temp3,'(1x,i10)') esdlim_
            if (iabs(esdlim_).lt.9.or.
     *        iabs(esdlim_).gt.99999) go to 900
        endif
        if (iwant.eq.7) then
          do ii = 1,ll
          if (temp(ii:ii).eq.'_') temp(ii:ii)=' '
          enddo
          result=prefx_(temp(1:ll),ll)
        endif
        if (iwant.eq.8) then
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
              tabl_=.false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              tabl_=.true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.9) then
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
        if (iwant.eq.10) reqlst = temp(1:ll)
        if (iwant.eq.11) then
          if (temp(1:ll).eq.'inline') then
            usedtd = 'inl'
            lspcdtd = 0
          else
            if (temp(1:ll).eq.'referto') then
              usedtd = 'ref'
              iwant = 12
              go to 100
            else
              if (temp(1:ll).eq.'write') then
                usedtd = 'wri'
                call c2xwarn(' writing of DTD not supported yet')
                iwant = 12
                go to 100
              else
                go to 900
              endif
            endif
          endif
        endif
        if (iwant.eq.12) spcdtd = temp(1:ll)
        if (iwant.eq.13) then
          if (temp(1:ll).eq.'drop' .or. temp(1:ll).eq.'insert' ) then
            unknown = temp(1:ll)
          else
            go to 900
          endif          
        endif
        if (iwant.eq.14) then
          if (temp(1:ll).eq.'row' .or. temp(1:ll).eq.'col' ) then
            by = temp(1:ll)
          else
            go to 900
          endif
        endif
        if (iwant.eq.15) then
          if (temp(1:ll).eq.'xfer' .or.
     *      temp(1:ll).eq.'keep' .or.
     *      temp(1:ll).eq.'zap' ) then
            xferesd = temp(1:ll)
          else
            go to 900
          endif
        endif
        iwant = 0
      else
        if (temp(1:1).eq.'-') then
          temp2=temp(3:256)
          if (ll.gt.2) then
            backarg = .true.
            karg = karg-1
          endif
          if (temp(2:2).eq.'i') then
            iwant = 1
            if (isi.gt.0) go to 900
            isi = 1
          endif
          if (temp(2:2).eq.'o') then
            iwant = 2
            if (iso.gt.0) go to 900
            iso = 1
          endif
          if (temp(2:2).eq.'d') then
            iwant = 3
            isd = 1
          endif
          if (temp(2:2).eq.'a') iwant = 4
          if (temp(2:2).eq.'f') iwant = 5
          if (temp(2:2).eq.'e') iwant = 6
          if (temp(2:2).eq.'p') iwant = 7
          if (temp(2:2).eq.'t') iwant = 8
          if (temp(2:2).eq.'c') iwant = 9
          if (temp(2:2).eq.'q') then
            iwant = 10
            if (isq.gt.0) go to 900
            isq = 1
          endif
          if (temp(2:2).eq.'s') then
            iwant = 11
            if (iss.gt.0) go to 900
            iss = 1
          endif
          if (temp(2:2).eq.'u') iwant = 13
          if (temp(2:2).eq.'b') iwant = 14
          if (temp(2:2).eq.'x') iwant = 15
          if (iwant.eq.0) go to 900
        else
          ifound = ifound+1
          if (ifound.eq.1) then
            if (isi.gt.0) then
              ifound = ifound+1
            else
              inpcif = temp(1:ll)
              isi = 1
            endif
          endif
          if (ifound.eq.2) then
            if (iso.gt.0) then
              ifound = ifound+1
            else
              outxml = temp(1:ll)
              iso = 1
            endif
          endif
          if (ifound.eq.3) then
            if (isd.gt.0) then
              ifound = ifound+1
            else
              ckdict = temp(1:ll)
              isd = 1
            endif
          endif
          if (ifound.eq.4) then
            if (isq.gt.0) then
              ifound = ifound+1
            else
              reqlst = temp(1:ll)
              isq = 1
            endif
          endif
          if (ifound.eq.5) then
            if (iss.gt.0) then
              ifound = ifound+1
            else
              spcdtd = temp(1:ll)
              usedtd = 'ref'
              iss = 1
            endif
          endif
          if (ifound.gt.4) go to 900
        endif
      endif
      go to 100
 500  linpcif = max(1,nblen(inpcif))
      loutxml = max(1,nblen(outxml))
      lckdict = max(1,nblen(ckdict))
      lreqlst = max(1,nblen(reqlst))
      return
 900  write(iunerr,'(a)')
     *  ' cif2xml [-i input_cif] [-o output_cif] '//
     *           '[-d dictionary] [-c catck]',
     *  '         [-f command_file] [-e esdlim_] [-a aliaso_] '//
     *           '[-p prefix] [-t tabl_]',
     *  '         [-q request_list] '//
     *           '[-s {inline | referto spec_dtd| write spec_dtd}] '//
     *           '[-b {row|col}]',
     *  '         [-x {xfer|keep|zap] '// 
     *           '[input_cif [output_cif '//
     *           '[dictionary',
     *  '         [request_list '//
     *           '[spec_dtd]]]]]'
       write(iunerr,'(a)')
     *  ' input_cif defaults to $CIF2XML_INPUT_CIF or stdin',
     *  ' output_cif defaults to $CIF2XML_OUTPUT_XML or stdout',
     *  ' dictionary defaults to $CIF2XML_CHECK_DICTIONARY',
     *  ' multiple dictionaries may be specified ',
     *  ' request_list defaults to $CIF2XML_REQUEST_LIST',
     *  ' input_cif of "-" is stdin, output_cif of "-" is stdout',
     *  ' request_list of "-" is stdin',
     *  ' -e has integer values (e.g. 9, 19(default) or 29),',
     *  ' -c has values of t or 1 or y vs. f or 0 or n,',
     *  ' -a has values of t or 1 or y vs. f or 0 or n,',
     *  ' -p has a string value in which "_" is replaced by blank',
     *  ' -t has values of t or 1 or y vs. f or 0 or n (default f),',
     *  ' -s defaults to inline, -b defaults to col, -x defaults to zap'
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
       subroutine c2xerr(mess)
C
C      variant of ciftbx err routine for cif2xml
C
       character*(*) mess
       call c2xmsg('error',mess)
       stop
       end
C
C
       subroutine c2xwarn(mess)
C
C      variant of ciftbx warn for cif2xml
C
       character*(*) mess
       call c2xmsg('warning',mess)
       return
       end
C
C
       subroutine c2xmsg(flag,mess)
C
       integer    nblen
       include   'cif2xml.cmn'
       include   'ciftbx.cmv'
       character*(*) flag
       character*(*) mess
       character*(MAXBUF)  tline
       integer       ll,ls,ltry,ii,i
C
       tline = ' cif2xml'
       if (nblen(flag).gt.0) 
     *   tline = ' cif2xml'//' '//flag(1:nblen(flag))
       tline= tline(1:nblen(tline))//': '
     *   //outxml(1:loutxml)//' '
     *   //' line:'
       ll = max(1,nblen(tline))
       write(iunerr,'(a,i7)')tline(1:ll),precn_
       ll=len(mess)
       ls=1
100    if(ll-ls.le.79) then
         write(iunerr,'(1X,a)') mess(ls:ll)
         return
       else
         ltry = min(ll,ls+79)
         do ii = ls+1,ltry
         i = ltry-ii+ls+1
         if(mess(i:i).eq.' ') then
           write(iunerr,'(1X,a)') mess(ls:i-1)
           ls=i+1
           if(ls.le.ll) go to 100
           return
         endif
         enddo
         write(iunerr,'(1X,a)') mess(ls:ltry)
         ls=ltry+1
         if(ls.le.ll) go to 100
         return
       endif  
       end
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


