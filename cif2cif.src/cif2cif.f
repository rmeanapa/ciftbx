C
C
C    \ | /      \ | /
C     \|/        \|/ 
C    -- -->>>>>>-- --               c i f 2 c i f  ...... CIF COPY PROGRAM
C     /|\        /|\
C    / | \      / | \                          Version 2.0.0
C                                               29 Nov 2009
C
C
C
C 
C cif2cif  is a fortran program using CIFtbx4 to copy a CIF on standard
C -------  input to an equivalent CIF on standard output, while checking
C          data names against dictionaries and reformating numbers with
C          esd's to conform to the rule of 19.  A quasar-style request
C          list may be specified, otherwise the entire CIF is copied.
C
C           cif2cif
C                  by
C
C                  Copyright (C) 1997, 1998, 2000, 2009
C                  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
C                  Bernstein + Sons
C                  5 Brewster Lane
C                  Bellport, NY 11713, U.S.A.
C
C           based on a suggestion by
C
C                  Sydney R. Hall (syd@crystal.uwa.edu.au)
C                  Crystallography Centre
C                  University of Western Australia
C                  Nedlands 6009, AUSTRALIA
C
C           with the request list handling suggested by the program
C           QUASAR by Sydney R. Hall and
C
C                  Rolf Sievers (unc411@dbnrhrz1.bitnet)
C                  Institut fuer Anorganische-Chemisches der Universitaet
C                  Gerhard-Domagk-Str.
C                  Bonn, GERMANY
C
C This is a copyrighted open source program subject to the license
C conditions in the file NOTICE.html
C
C cif2cif reads the input CIF from the standard input device (normally
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
C cif2cif [-i input_cif] [-o output_cif] [-d dictionary] [-a aliaso_]\
C         [-c catck] [-e esdlim_] [-f command_file] [-g guess_type]\
C         [-m maxline] [-p prefix] [-q request_list] [-t tabl_]\
C         [-u unfold] [-w wrap] [-B read|cif2read] \
C         [input_cif [output_cif [dictionary [request_list]]]] 
C
C where:
C         input_cif defaults to $CIF2CIF_INPUT_CIF or stdin
C         output_cif defaults to $CIF2CIF_OUTPUT_CIF or stdout
C         dictionary defaults to $CIF2CIF_CHECK_DICTIONARY
C           (multiple dictionaries may be specified)
C         request_list defaults to $CIF2CIF_REQUEST_LIST
C         input_cif of "-" is stdin, output_cif of "-" is stdout
C         request_list of "-" is stdin
C         -a has values of t or 1 or y vs. f or 0 or n
C         -c has values of t or 1 or y vs. f or 0 or n
C         -e has integer values (e.g. 9, 19(default) o 29)
C         -g has values of t or 1 or y vs. f or 0 or n, default t
C         -m has values from 80 to 2048 for a maximum line width
C         -p has string values in which "_" is replaced by blank
C         -t has values of t or 1 or y vs. f or 0 or n, default f
C         -u has values of t or 1 or y vs. f or 0 or n (default f)
C         -w has values of 0 for no wrap or a column to wrap on
C
C
C .............................Installation notes.............................
C
C
C cif2cif is provided as a part of the CIFtbx4 release kit, as the
C files 'cif2cif.f' and 'cif2cif.cmn'.  CIFtbx4 is needed for the files 
C ciftbx.f, ciftbx.sys and hash_funcs.f.  See the CIFtbx4 documentation 
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
C         if(evar.eq.'CIF2CIF_INPUT_CIF')
C        *  string='INPCIF.CIF'//char(0)
C         if(evar.eq.'CIF2CIF_OUTPUT_CIF')
C        *  string='OUTCIF.CIF'//char(0)
C         if(evar.eq.'CIF2CIF_CHECK_DICTIONARY') 
C        *  string='CIF_CORE.DIC'//char(0)
C         if(evar.eq.'CIF2CIF_REQUEST_LIST') 
C        *  string='REQLST.DAT'//char(0)
C         return
C         end
C
C This combination of substitute routines would "wire-in" cif2cif to
C read its input cif from a file named INPCIF.CIF, write its output
C cif to a file named OUTCIF.CIF, check names against CIF_CORE.DIC
C and use the tag names given in REQLST.DAT to selects the ones to copy. 
C
C ..................................Files Used................................
C
C  dictionary input         input   on device 2
C  Reformatted CIF          output  on device 6 ('stdout')
C  Input CIF                input   on device 2, if a file, 5  if 'stdin'
C  Message device           output  on device 0 ('stderr')
C  Direct access            in/out  on device 3
C  Request list             input   on device 4, if a file, 5 if 'stdin'
C
C
C***************************************************************************
C
      program cif2cif
C
      character*28 VERS
      integer  lastnb
      include 'cif2cif.cmn'
      include 'ciftbx.sys'
      include 'ciftbx.cmf'
      character*(MAXBUF) tbxxlocs
      character*(MAXBUF) string,temp,blok
      integer ipcnt(NUMBLOCK), ipdup(NUMBLOCK)
      character*3 xvc,xtc
      character*1 xxtab,c
      character*4 mytype
      character*1 mydelim
      integer ixname,iyname,ii,ipdata,ipsave,kkr,kkrl
      integer kreqnam,lreqnam,kpass,markdb,kbstrt,lbstrt
      integer nblen,kkp,krr,krrl,jpsave
      integer kfold, kprevpos, knextpos, krord
      integer basedepth,ndepth,pposdlm,recdlm
      double precision dpn,dpesd
C
      VERS='2.0.0 (29 Nov 2009)'
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
      tabl_ = .false.
      pfold_ = 0
      unfold_ = .false.
      align_ = .false.
      tabx_ = .false.
      ptabx_ = .false.
      pdecp_ = .false.
      plzero_ = .false.
      guesst = .true.
      kbstrt = -1
      kprevpos = 0
      knextpos = 8
      koffset = 0
      basedepth = 0
      krr = 0
      result = init_(iuninp,iunout,iundac,iunerr)
      call getfiles
      kfold = line_
      if (pfold_.gt.0 .and. pfold_.lt.line_) 
     * kfold = min(line_,max(20,pfold_))
      if (tabl_) align_ = .true.
      isline = line_
      line_ = 2048
      if (ckdict(1:lckdict).ne.' ')
     *  result = dict_(ckdict(1:lckdict),'valid dtype '//cats)
      line_ = isline
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
        call c2cerr(
     *  ' failed to open '//inpcif(1:linpcif))
      endif
      if (nblen(reqlst).gt.0) then
        if (reqlst(1:lreqlst).eq.'-') then
          iunrql = 5
        else
          open(unit=iunrql,file=reqlst(1:lreqlst),status='old',
     *      err=150)
          go to 160
150       call c2cerr(' Unable to open request list '//
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
        else
          call c2cerr(' More than NUMREQ items in request list')
        endif
        kkp = 1
        go to 170
180     if (reqlst(1:lreqlst).ne.'-') close(unit=iunrql)
      endif
190   continue
C
C     open output CIF
C
      if (outcif(1:loutcif).eq.'-')  outcif = ' '
      result=pfile_(outcif(1:loutcif))
      if(.not.result) then
        call c2cerr(
     *  ' failed to open '//outcif(1:loutcif))
      endif
      result=pcmnt_(' CIF Copied by cif2cif, version '// VERS)
      result=find_(' ','head',string)
      ipsave=0
      markdb=0
      kreqnam = 1
C
C     This is the main loop.  return to 199 for a full new pass
C
199   kpass = 1
      lreqnam = kreqnam
      lbstrt = kbstrt
      blok = ' '
      if(markdb.ne.0 .and.nreqnam.ne.0)result=bkmrk_(markdb)
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
              call c2cwarn(' Did not find '//
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
        if (reqnam(kreqnam)(krrl:krrl).eq.'_') then
          do ii = 1,nname
          if (dname(ii)(1:krrl).eq.
     *      reqnam(kreqnam)(1:krrl) ) then
            ipcnt(ii) = ipcnt(ii)+1
            if (krr.gt.NUMCHRR)
     *        call c2cerr(' Number of tags greater than NUMCHRR')
            irord(krr) = ii
            krr = krr+1
          endif
          enddo
        else
          call hash_find(reqnam(kreqnam)(1:krrl),
     *     dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *     ii)
          if (ii.gt.0) then
            if (krr.gt.NUMCHRR)
     *        call c2cerr(' Number of tags greater than NUMCHRR')
            irord(krr) = ii
            krr = krr+1
          else
            call c2cwarn(' No match to requested name '//
     *      reqnam(kreqnam)(1:krrl))
            if (krr.gt.NUMCHRR)
     *        call c2cerr(' Number of tags greater than NUMCHRR')
            irord(krr) = -kreqnam
            krr = krr+1
          endif
        endif
        go to 210
C
220     krr = krr-1
        if (krr.le.0) then
          if (blok.ne.'which_contains:'.or.
     *      (kpass.eq.2.and.kbstrt.ge.lbstrt)) then
            call c2cwarn(' No tags requested from block '//bloc_)
            go to 199
          else
            kreqnam = lreqnam
            go to 200
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
      result=pdata_(string)
      if(.not.result) then
        call c2cwarn(
     *  ' duplicate data block '//bloc_)
        pposnam_ = 40
        result=pcmnt_
     *    ('<---- duplicate data block')
      endif
      ixname=1
      if (nreqnam.eq.0) krr = nname
250   if(ixname.le.krr) then
        loop_ = .false.
        krord = 0
        do ii = 1,nname
        ipdup(ii) = 0
        enddo
        if (irord(ixname).gt.0) krord = nloop(irord(ixname)) 
        if(krord.ne.0) then
        ipdup(irord(ixname)) = 1
255       if(jpsave.eq.0 .and.nreqnam.eq.0) call cpcmnt
          jpsave = 0
          if (nreqnam.eq.0) then
          result=find_(' ',' ',string)
          if(type_.ne.'loop') then
            call c2cwarn(' Expected loop_, found: '//string)
            go to 255
          endif
          else
          posval_ = 1
          endif
          pposval_=posval_
          result=ploop_(' ')
          iyname=ixname
          if (irord(ixname).lt.0) go to 260
          do ii = ixname+1,krr
            if (irord(ii).lt.0) then
              krord=nloop(irord(ixname))
            else
              krord=nloop(irord(ii))
            endif
            if (nloop(irord(ixname)).ne.
     *        krord) go to 260
            if (irord(ii).gt.0) then
              if (ipdup(irord(ii)).gt.0) go to 260
            endif
            iyname=ii
            if (irord(ii).gt.0) then
              ipdup(irord(ii)) = 1
            endif
          enddo
260       continue
          if (nreqnam .ne.0) result = pcmnt_(' ')
          do ii = ixname,iyname
            if (nreqnam.eq.0) then
            call cpcmnt
            result=find_(' ','name',string)
            if(tbxxlocs(string).ne.dname(irord(ii))) then
              call c2cwarn(' Expected '
     *         //dname(irord(ii))(1:lastnb(dname(irord(ii))))
     *         //', found '
     *         //string(1:lastnb(string)))
            endif
            else
              if (irord(ii).gt.0) then
                result=find_(dname(irord(ii)),'name',string)
              endif
            endif
            if (irord(ii).gt.0) then
              if(.not.tabl_) pposnam_=posnam_
              result=ploop_(string(1:max(1,lastnb(string))))
            else
              if(.not.tabl_) pposnam_ = 1
              result = ploop_(reqnam(-irord(ii)))
              pposnam_ = 40
              result=pcmnt_
     *         ('<---- requested item not present')
            endif
          enddo
          if (nreqnam.ne.0) then
          result = bkmrk_(markdb)
          result = bkmrk_(markdb)
          endif
270       do ii = ixname,iyname
            pquote_ = ' '
            if (nreqnam.eq.0) then
              if (unfold_.and.kprevpos.ne.0) then
                call cpcmnt2
              else
                call cpcmnt
              endif
            endif
            if (irord(ii).gt.0) then
              testfl = 'no '
              result=test_(dname(irord(ii)))
              basedepth = 0
CDBG          print *,' Calling test_ ',dname(irord(ii))
              if (type_.eq.'numb'.and.(.not.guesst)) then
                if (.not.dtype_(dname(irord(ii)),mytype)) type_='char'
              endif
CDBG          print *,' post test_ type_ ',type_
            else
              type_ = 'null'
            endif
            if (ii .eq. ixname .and. posval_ .lt. 17) then
              kprevpos=0
              knextpos=8
              koffset=0
            endif
            if(.not.tabl_) then
              pposval_=mod(posval_-koffset-1,kfold)+1
              if (unfold_) then
                pposval_=max(8-mod(kprevpos+7,8)
     *            +kprevpos,posval_-koffset)
                pposval_=max(pposval_,knextpos)
                pposval_=mod(pposval_-1,kfold)+1
              endif
              pposend_=0
              if (posend_.gt.0)
     *          pposend_=pposval_+(posend_-posval_)
            endif
275         continue
            if (depth_.le.basedepth .and.depth_.gt.0) then
              if (delim_(depth_,mydelim,pposdlm,recdlm)) then
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm) 
              endif
              basedepth = depth_
            end if
            if(type_.eq.'null') then
            testfl = 'no '
            if (depth_.gt.basedepth) then
              do ndepth = basedepth+1,depth_
              if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm)
              endif
              enddo
              basedepth = depth_
            endif
            if(.not.tabl_) then
              pposval_=mod(posval_-koffset-1,kfold)+1
              if (unfold_) then
                pposval_=max(8-mod(kprevpos+7,8)
     *            +kprevpos,posval_-koffset)
                pposval_=max(pposval_,knextpos)
                pposval_=mod(pposval_-1,kfold)+1
                if(pposval_+long_.gt.kfold) then
                  kprevpos=0
                  knextpos=1
                  pposval_=1
                  koffset=posval_-1
                endif
                kprevpos=pposval_
                koffset=posval_-pposval_
                knextpos=mod(kprevpos+3,kfold)+1
              endif
              pposend_=0
              if (posend_.gt.0) then
                 pposend_=pposval_+(posend_-posval_)
                 knextpos=mod(max(knextpos,pposend_+2)-1,kfold)+1
              endif
            endif
              if (irord(ii).lt.0) then
                 result=pchar_(' ','?')
              else
              if (long_.eq.1.and.strg_(1:1).eq.'?') then
                 result=pchar_(' ','?')
              else
                 result=pchar_(' ','.')
              endif
              endif
              if (depth_.eq.0) then
                type_=' '
              else
                call cpcmnt
                if (depth_.gt.0) then
                  result=test_(dname(irord(ii)))
                  go to 275
                else
                  type_=' '
                end if
              end if
            endif
            if(type_.eq.'char') then
CDBG          print *,' handling as char at basedepth ',basedepth
              result=char_(dname(irord(ii)),string)
              if (depth_.gt.basedepth) then
                do ndepth = basedepth+1,depth_
                if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                  result = pdelim_(mydelim,.false.,pposdlm)
                  call cpbcmnt(pposdlm,recdlm)
                endif
                enddo
                basedepth = depth_
              endif
              pquote_=quote_
              if(.not.tabl_) then
              pposval_=mod(posval_-koffset-1,kfold)+1
              if (unfold_) then
                pposval_=max(8-mod(kprevpos+7,8)
     *            +kprevpos,posval_-koffset)
                pposval_=max(pposval_,knextpos)
                pposval_=mod(pposval_-1,kfold)+1
                  if(pposval_+long_.gt.kfold) then
                  kprevpos=0
                  knextpos=1
                  pposval_=1
                  koffset=posval_-1
                endif
                kprevpos=pposval_
                koffset=posval_-pposval_
                knextpos=mod(kprevpos+3,kfold)+1
              endif
              pposend_=0
              if (posend_.gt.0) then
                 pposend_=pposval_+(posend_-posval_)
                 knextpos=mod(max(knextpos,pposend_+4)-1,kfold)+1
              endif
              endif
              if(string.eq.'?'.or.string.eq.'.') then
                if(pquote_.eq.' ') pquote_=''''
              endif
              if (depth_.eq.0) then
                result=pchar_(' ',string(1:max(1,lastnb(string))))
                type_=' '
              else
                result=pchar_(char(0),string(1:max(1,lastnb(string))))
                call cpcmnt
                if (depth_.gt.0) then
                  result=test_(dname(irord(ii)))
                  go to 275
                else
                  type_=' '
                end if
              end if
            endif
            if(type_.eq.'text') then
280           result=char_(dname(irord(ii)),string)
              if (depth_.gt.basedepth) then
                do ndepth = basedepth+1,depth_
                if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                  result = pdelim_(mydelim,.false.,pposdlm)
                  call cpbcmnt(pposdlm,recdlm)
                endif
                enddo
                basedepth = depth_
              endif
              if (depth_.eq.0) then
                if (type_.eq.'text') 
     *            result=ptext_(' ',string(1:max(1,lastnb(string))))
                if(text_) go to 280
                call tbxxeot
                type_=' '
                kprevpos = 0
              else
                if (type_.eq.'text') 
     *          result=ptext_(char(0),string(1:max(1,lastnb(string))))
                if(text_) go to 280
                call tbxxeot
                call cpcmnt
                if (depth_.gt.0) then
                  result=test_(dname(irord(ii)))
                  go to 275
                else
                  type_=' '
                end if
              end if
            endif
            if(type_.eq.'numb') then
CDBG          print *,' processing numb', basedepth
              dpesd = 0.D0
              result=numd_(dname(irord(ii)),dpn,dpesd)
              if (depth_.gt.basedepth) then
                do ndepth = basedepth+1,depth_
                if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
CDBG              print *,' numb delim ',mydelim
                  result = pdelim_(mydelim,.false.,pposdlm)
                  call cpbcmnt(pposdlm,recdlm)
                endif
                enddo
                basedepth = depth_
              endif
              if(.not.tabl_) then
              pposval_=mod(posval_-1-koffset,kfold)+1
              if (unfold_) then
                pposval_=max(8-mod(kprevpos+7,8)
     *            +kprevpos,posval_-koffset)
                pposval_=max(pposval_,knextpos)
                pposval_=mod(pposval_-1,kfold)+1
                if(pposval_+long_.gt.kfold) then
                  kprevpos=0
                  knextpos=1
                  pposval_=1
                  koffset=posval_-1
                endif
                kprevpos=pposval_
                koffset=posval_-pposval_
                knextpos=mod(kprevpos+3,kfold)+1
              endif
              pposend_=0
              if (posend_.gt.0) then
                 pposend_=pposval_+(posend_-posval_)
                 knextpos=mod(max(knextpos,pposend_+2)-1,kfold)+1
              endif
              pposdec_=0
              if (posdec_.gt.0)
     *          pposdec_ = pposval_+(posdec_-posval_)
              endif
              pesddig_=esddig_
              pdecp_=decp_
              plzero_=lzero_
              if (depth_.eq.0) then
                result=pnumd_(' ',dpn,dpesd)
                type_=' '
              else
                result=pnumd_(char(0),dpn,dpesd)
                call cpcmnt
                if (depth_.gt.0) then
                  result=test_(dname(irord(ii)))
                  go to 275
                else
                  type_=' '
                end if
              end if
            endif
          enddo
          if (nreqnam.ne.0) result=pcmnt_(char(0))
          if(loop_) go to 270
          ixname=iyname+1
          if (nreqnam.ne.0) result=pcmnt_(' ')
          go to 250
        else
          if(jpsave.eq.0 .and.nreqnam.eq.0) call cpcmnt
          jpsave=0
          if (nreqnam.eq.0) then
          result=name_(temp)
          else
          if (irord(ixname).lt.0) then
            temp = reqnam(-irord(ixname))
            if (.not.tabl_) pposnam_ = 1
            result = pchar_(temp,'?')
            pposnam_ = 40
            result=pcmnt_
     *      ('<---- requested item not present')
            ixname = ixname+1
            go to 250
          endif
          result =find_(dname(irord(ixname)),'name',temp)
          endif
          if(.not.result)
     *      call c2cerr(' misaligned data item '//
     *      dname(irord(ixname)))
          if (nreqnam.eq.0) call cpcmnt
          basedepth = 0
          result=test_(dname(irord(ixname)))
285       if (type_.eq.'numb'.and.(.not.guesst)) then
              if (.not.dtype_(dname(irord(ixname)),mytype)) type_='char'
          endif
          if (depth_.le.basedepth .and.depth_.gt.0) then
            if (delim_(depth_,mydelim,pposdlm,recdlm)) then
              result = pdelim_(mydelim,.false.,pposdlm)
              call cpbcmnt(pposdlm,recdlm)
            endif
            basedepth = depth_
          end if
          if(type_.eq.'null') then
            result=char_(dname(irord(ixname)),string)
            if (depth_.gt.basedepth) then
              if (basedepth.eq.0) then
                pposnam_=0
                if (posnam_.gt.0)
     *            pposnam_ = pposval_+(posnam_-posval_)
                if (pposnam_ .lt. 0) pposnam_ = 1
                result=pchar_(temp,char(0))
              end if
              do ndepth = basedepth+1,depth_
              if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm)
              endif
              enddo
              basedepth = depth_           
            endif
            if (.not.tabl_) then
              pposval_=mod(posval_-1,kfold)+1
              pposend_=0
              if (posend_.gt.0)
     *          pposend_=pposval_+(posend_-posval_)
              pposnam_=0
              if (posnam_.gt.0)
     *          pposnam_ = pposval_+(posnam_-posval_)
              if (pposnam_ .lt. 0) pposnam_ = 1
            endif
            if (long_.eq.1.and.strg_(1:1).eq.'?') then
               if (depth_.eq. 0) then
                 result=pchar_(temp,'?')
               else
                 result=pchar_(char(0),'?')
               end if
            else
               if (depth_.eq. 0) then
                 result=pchar_(temp,'.')
               else
                 result=pchar_(char(0),'.')
               end if
            endif
            if (depth_.eq.0) then 
              type_=' '
            else
              call cpcmnt
              if (depth_.gt.0) then
                result = test_(dname(irord(ixname)))
                if (.not.(type_.eq.'null'.and.depth_.eq.0))
     *          go to 285
              end if
              type_=' '
            endif
          endif  
          if(type_.eq.'char') then
            result=char_(dname(irord(ixname)),string)
CDBG        print *,' item as read: ',string(1:lastnb(string))
            if (depth_.gt.basedepth) then
              if (basedepth.eq.0) then
                pposnam_=0
                if (posnam_.gt.0)
     *            pposnam_ = pposval_+(posnam_-posval_)
                if (pposnam_ .lt. 0) pposnam_ = 1
                result=pchar_(temp,char(0))
              end if
              do ndepth = basedepth+1,depth_
              if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
CDBG            print *, ' ndepth,mydelim,pposdlm ',
CDBG *            ndepth,mydelim,pposdlm
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm)
              endif
              enddo
              basedepth = depth_
            endif
            pquote_=quote_
            if (.not.tabl_) then
              pposval_=mod(posval_-1,kfold)+1
              pposend_=0
              if (posend_.gt.0)
     *          pposend_=pposval_+(posend_-posval_)
              pposnam_=0
              if (posnam_.gt.0)
     *          pposnam_ = pposval_+(posnam_-posval_)
              if (pposnam_ .lt. 0) pposnam_ = 1
            endif
            if(string.eq.'?'.or.string.eq.'.') then
              if(pquote_.eq.' ') pquote_=''''
            endif
            if (depth_.eq.0) then
              result=pchar_(temp,string(1:max(1,lastnb(string))))
              type_=' '
            else
              result=pchar_(char(0),string(1:max(1,lastnb(string))))
CDBG          print *, 'depth before cpcmnt ', depth_
              call cpcmnt
              if (depth_.gt.0) then
                result = test_(dname(irord(ixname)))
                if (.not.(type_.eq.'null'.and.depth_.eq.0))
     *          go to 285
              end if
              type_=' '
            end if
          endif
          if(type_.eq.'text') then
290         result=char_(dname(irord(ixname)),string)
            if (depth_.gt.basedepth) then
              if (basedepth.eq.0) then
                pposnam_=0
                if (posnam_.gt.0)
     *            pposnam_ = pposval_+(posnam_-posval_)
                if (pposnam_ .lt. 0) pposnam_ = 1
                result=pchar_(temp,char(0))
              end if
              do ndepth = basedepth+1,depth_
              if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm) 
              endif
              enddo
              basedepth = depth_
            endif
            if (.not.tabl_) pposnam_=posnam_
            if (depth_.eq.0) then
              if (type_.eq.'text')
     *        result=ptext_(temp,string(1:max(1,lastnb(string))))
              if(text_) go to 290
              call tbxxeot
              type_=' '
            else
              if (type_.eq.'text')
     *        result=ptext_(char(0),string(1:max(1,lastnb(string))))
              if (text_) go to 290
              call tbxxeot
CDBG          print *, 'depth before cpcmnt ', depth_
              call cpcmnt
              if (depth_.gt.0) then
                result = test_(dname(irord(ixname)))
                if (.not.(type_.eq.'null'.and.depth_.eq.0))
     *          go to 285
              end if
              type_=' '
            end if
          endif
          if(type_.eq.'numb') then
            dpesd=0.D0
            result=numd_(dname(irord(ixname)),dpn,dpesd)
            if (depth_.gt.basedepth) then
              if (basedepth.eq.0) then
                pposnam_=0
                if (posnam_.gt.0)
     *            pposnam_ = pposval_+(posnam_-posval_)
                if (pposnam_ .lt. 0) pposnam_ = 1
                result=pchar_(temp,char(0))
              end if
              do ndepth = basedepth+1,depth_
              if (delim_(ndepth,mydelim,pposdlm,recdlm)) then
                result = pdelim_(mydelim,.false.,pposdlm)
                call cpbcmnt(pposdlm,recdlm) 
              endif
              enddo
              basedepth = depth_
            endif
            if (.not.tabl_) then
              pposval_=mod(posval_-1,kfold)+1
              pposend_=0
              if (posend_.gt.0)
     *          pposend_=pposval_+(posend_-posval_)
              pposdec_ = 0
              if (posdec_.gt.0)
     *          pposdec_=pposval_+(posdec_-posval_)
              if (posnam_.gt.0)
     *          pposnam_ = pposval_+(posnam_-posval_)
              if (pposnam_ .lt. 0) pposnam_ = 1
            endif
            pdecp_=decp_
            plzero_=lzero_
            pesddig_=esddig_
            if (depth_.eq.0) then
              result=pnumd_(temp,dpn,dpesd)
              type_=' '
            else
              result=pnumd_(char(0),dpn,dpesd)
CDBG          print *, 'depth before cpcmnt ', depth_
              call cpcmnt
              if (depth_.gt.0) then
                result = test_(dname(irord(ixname)))
                if (.not.(type_.eq.'null'.and.depth_.eq.0))
     *          go to 285
              end if
              type_=' '
            end if
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
      subroutine cpbcmnt(posdlm,recdlm)
C
C     Fudge to back up and copy the comments, if any,
C     after a delimter and prior to the current position
C
      include 'cif2cif.cmn'
      integer  posdlm
      integer  recdlm
      logical  cotdb_
      logical  pcmnt_
      logical  isdlm
      logical  pdelim_
      integer  lstring
      include 'ciftbx.cmv'
      character*(MAXBUF) string
100   result=cotdb_(string,lstring,isdlm,posdlm,recdlm)
      if (result) then
        if (isdlm) then
          result = pdelim_(string(1:1),.false.,posval_)
        else
          pposnam_ = posnam_
          pquote_=quote_
          result=pcmnt_(string(1:lstring))
        end if
        goto 100
      endif
      return
      end
C
C
      subroutine cpcmnt
C
C     routine to copy a set of comments (if any) from the
C     input cif to the output cif
C
      include 'cif2cif.cmn'
      logical  cotd_
      logical  pcmnt_
      logical  isdlm
      logical  pdelim_
      include 'ciftbx.cmv'
      character*(MAXBUF) string
      call tbxxbtab
100   result=cotd_(string,isdlm)
      if (result) then
        if (isdlm) then
          result = pdelim_(string,.false.,posval_)
        else
          pposnam_ = posnam_
          pquote_=quote_
          result=pcmnt_(string(1:long_))
        end if
        goto 100
      endif
      call tbxxetab
      return
      end
C
C
      subroutine cpcmnt2
C
C     routine to copy a set of comments (if any) from the
C     input cif to the output cif within a loop_ row
C     do not copy empty comments
C
      include 'cif2cif.cmn'
      logical  cotd_
      logical  pcmnt_
      logical  isdlm
      logical  pdelim_
      external pdelim_
      include 'ciftbx.cmv'
      character*(MAXBUF) string
      character*1 xxtab
      xxtab = char(5)
      if (ichar('a')-ichar('A') .gt.0) xxtab = char(9)
      call tbxxbtab
100   result=cotd_(string,isdlm)
      if (result) then
        if(isdlm) then
          result = pdelim_(string,.false.,posval_)
          go to 100
        else
        if(quote_.eq.'#' .or.
     *    ( string(1:long_).ne.char(0) .and.
     *      string(1:long_).ne.xxtab.and.
     *      string(1:long_).ne.' ')) then
          pposnam_ = posnam_
          pquote_=quote_
          result=pcmnt_(string(1:long_))
          goto 100
        endif
        endif
      endif
      call tbxxetab
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
      include 'cif2cif.cmn'
      include 'ciftbx.cmv'
      character*256 temp,temp2,cline
      character*10 bfill
      character*11 temp3
      character*2 temp4
      integer iargc,ll,karg,kfarg,nfarg,
     *  ifound,iwant,isi,iso,isd,isq,ii
      logical backarg,ffile
      data bfill/'          '/
      numarg = iargc()
      call getenv('CIF2CIF_INPUT_CIF',inpcif)
      call getenv('CIF2CIF_OUTPUT_CIF',outcif)
      call getenv('CIF2CIF_CHECK_DICTIONARY',ckdict)
      call getenv('CIF2CIF_REQUEST_LIST',reqlst)
      cats = 'catno'
      karg = 0
      kfarg = 0
      iwant = 0
      ifound = 0
      isi = 0
      iso = 0
      isd = 0
      isq = 0
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
        if (iwant.eq.2) outcif = temp(1:ll)
        if (iwant.eq.3) then
          ckdict = temp(1:ll)
          lckdict = max(1,nblen(ckdict))
          isline = line_
          line_ = 2048
          result = dict_(ckdict(1:lckdict),'valid dtype '//cats)
          line_ = isline
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
          if(temp(1:ll).eq.'F.' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            unfold_ = .false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              unfold_ = .true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.12) then
          if(ll.gt.10) go to 900
          temp3=bfill(1:11-ll)//temp(1:ll)
          read(temp3,'(1x,i10)') pfold_
            if (pfold_.lt.0.or.
     *        pfold_.gt.99999) go to 900
        endif
        if (iwant.eq.13) then
          if(ll.gt.10) go to 900
          temp3=bfill(1:11-ll)//temp(1:ll)
          read(temp3,'(1x,i10)') line_
            if (line_.lt.80.or.
     *        line_.gt.MAXBUF) go to 900
        endif
        if (iwant.eq.14) then
          if(temp(1:ll).eq.'F' .or. temp(1:ll).eq.'f' .or.
     *      temp(1:ll).eq.'0' .or. temp(1:ll).eq.'N' .or.
     *      temp(1:ll).eq.'n') then
            guesst = .false.
          else
            if(temp(1:ll).eq.'T.' .or. temp(1:ll).eq.'t' .or.
     *      temp(1:ll).eq.'1' .or. temp(1:ll).eq.'Y' .or.
     *      temp(1:ll).eq.'y') then
              guesst = .true.
            else
              go to 900
            endif
          endif
        endif
        if (iwant.eq.15) then
          if (temp(1:ll).eq.'read' .or. 
     *      temp(1:ll).eq.'READ') then
              rdbkt_ = .true.
              rdprn_ = .true.
              rdbrc_ = .true.
          else
            if (temp(1:ll).eq.'cif2read' .or. 
     *      temp(1:ll).eq.'CIF2READ' .or.
     *      temp(1:ll).eq.'CIF2read') then
               rdbrc_ = .true.
               rdbkt_ = .true.
               rdcolon_ = .true.
               rdrcqt_=.true.
               rdtq_=.true.
            else
              if (temp(1:ll).eq.'noread' .or.
     *        temp(1:ll).eq.'NOREAD') then
                rdbkt_ = .false.
                rdprn_ = .false.
                rdbrc_ = .false.  
              else
                if (temp(1:ll).eq.'nocif2read' .or. 
     *          temp(1:ll).eq.'NOCIF2READ' .or.
     *          temp(1:ll).eq.'NOCIF2read') then
                  rdbrc_ = .false.
                  rdbkt_ = .false.
                  rdcolon_ = .false.
                  rdrcqt_=.false.
                  rdtq_=.false.
                else
                  go to 900
                end if
              end if
            end if
          end if
        end if
        iwant = 0
      else
        if (temp(1:1).eq.'-') then
          temp2=temp(3:256)
          if (ll.gt.2) then
            backarg = .true.
            karg = karg-1
          else
            if (karg.lt.numarg) then
              call getarg(karg+1,temp4)
            else
              temp4 = '--'
            endif
            if (temp4(1:1).eq.'-'
     *        .and.temp4(2:2).ne.' ') then
              temp2='y'
              if (temp(2:2).eq.'i' .or. temp(2:2).eq.'o') temp2='-'
              if (temp(2:2).eq.'d') go to 900
              if (temp(2:2).eq.'e') temp2='19'
              if (temp(2:2).eq.'f') go to 900
              if (temp(2:2).eq.'m') temp2 = '80'
              if (temp(2:2).eq.'p') go to 900
              if (temp(2:2).eq.'q') go to 900
              if (temp(2:2).eq.'w') temp2 = '80'
              backarg = .true.
              karg = karg-1
            endif
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
          if (temp(2:2).eq.'u') iwant = 11
          if (temp(2:2).eq.'w') iwant = 12
          if (temp(2:2).eq.'m') iwant = 13
          if (temp(2:2).eq.'g') iwant = 14
          if (temp(2:2).eq.'B') iwant = 15
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
              outcif = temp(1:ll)
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
          if (ifound.gt.4) go to 900
        endif
      endif
      go to 100
 500  linpcif = max(1,nblen(inpcif))
      loutcif = max(1,nblen(outcif))
      lckdict = max(1,nblen(ckdict))
      lreqlst = max(1,nblen(reqlst))
      return
 900  write(iunerr,'(a)')
     *  ' cif2cif [-i input_cif] [-o output_cif] '//
     *           '[-d dictionary] [-c catck]',
     *  '         [-a aliaso_] [-e esdlim_] [-f command_file] '//
     *           '[-g guess_type] [-m maxline] [-p prefix]',
     *  '         [-q request_list] [-t tabl_]  [-u unfold] '//
     *           '[-w wrap] [-B read|cif2read] ',
     *  '         [input_cif [output_cif '//
     *           '[dictionary [request_list]]]]'
       write(iunerr,'(a)')
     * ' input_cif defaults to $CIF2CIF_INPUT_CIF or stdin',
     * ' output_cif defaults to $CIF2CIF_OUTPUT_CIF or stdout',
     * ' dictionary defaults to $CIF2CIF_CHECK_DICTIONARY',
     * ' multiple dictionaries may be specified ',
     * ' request_list defaults to $CIF2CIF_REQUEST_LIST',
     * ' input_cif of "-" is stdin, output_cif of "-" is stdout',
     * ' request_list of "-" is stdin'
       write(iunerr,'(a)')
     * ' -a has values of t or 1 or y vs. f or 0 or n,',
     * ' -c has values of t or 1 or y vs. f or 0 or n,',
     * ' -e has integer values (e.g. 9, 19(default) or 29),',
     * ' -g has values of t, 1 or y vs. f, 0 or n (default t),',
     * ' -m has values from 80 to 2048 for a maximum line width',
     * ' -p has a string value in which "_" is replaced by blank,',
     * ' -t has values of t, 1 or y vs. f, 0 or n (default f),',
     * ' -u has values of t, 1 or y vs. f, 0 or n (default f),',
     * ' -w has values of 0 for no wrap or a column to wrap on',
     * ' -B has values to read brackets for DDLm 2008 or CIF2' 
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
       subroutine c2cerr(mess)
C
C      variant of ciftbx err routine for cif2cif
C
       character*(*) mess
       call c2cmsg('error',mess)
       stop
       end
C
C
       subroutine c2cwarn(mess)
C
C      variant of ciftbx warn for cif2cif
C
       character*(*) mess
       call c2cmsg('warning',mess)
       return
       end
C
C
       subroutine c2cmsg(flag,mess)
C
       integer    nblen
       include   'cif2cif.cmn'
       include   'ciftbx.cmv'
       character*(*) flag
       character*(*) mess
       character*(MAXBUF)  tline
       integer       ll,ls,ltry,ii,i
C
       tline = ' cif2cif'
       if (nblen(flag).gt.0) 
     *   tline = ' cif2cif'//' '//flag(1:nblen(flag))
       tline= tline(1:nblen(tline))//': '
     *   //outcif(1:loutcif)//' '
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


