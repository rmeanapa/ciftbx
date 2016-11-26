       program cif2pdb
C
C      program for conversion of an mmCIF data set to a pseudo-PDB
C      entry or a WPDB entry
C
C      Version 2.0.3 30 November 2009
C
C      Herbert J. Bernstein
C      Bernstein+Sons
C      5 Brewster Lane
C      Bellport, NY 11713-2803
C      phone:  1-631-286-1339   email:  yaya@bernstein-plus-sons.com
C
C      Frances C. Bernstein
C      Bernstein+Sons
C      5 Brewster Lane
C      Bellport, NY 11713-2803
C      phone:  1-631-286-1339   email:  fcb@bernstein-plus-sons.com
C
C      This program is a version 2.0.1 of cif2pdb, capable of doing
C      a conversion of an mmCIF data set to a partial pseudo-PDB entry
C      with HEADER, TITLE, COMPND, SOURCE, KEYWRD, AUTHOR, JRNL,
C      REMARK 1, REMARK 960, SEQRES, CRYST1, ORIGX, SCALE, ATOM, 
C      SIGATM, ANISOU, SIGUIJ (starting with U's or B's), HETATM, 
C      MASTER and END records.  Fractional or orthogonal coordinates
C      may be provided in the mmCIF data set.  If an mmCIF dictionary
C      or other dictionary with the necessary aliases between mmCIF
C      and the core is provided, the program can convert core CIF
C      data sets.  This is sufficient to drive RASMOL.
C
C      With version 2, the alternative of using WPDB format
C      for output is provided, flagged by use of the -w
C      option on the command line and in the output by
C      the LEADER record in place of a HEADER record.
C
C      For more on WPDB format, see 
C         http://biomol.dowling.edu/WPDB
C
C      Comments to yaya@bernstein-plus-sons.com appreciated.
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external pdbstr
       external typeset
       logical init_
       logical ocif_
       logical dict_
       logical data_
       logical findtag_
C

C
C      Initialization of variables
C
       iunpdb = 6
       iuninp = 11
       iunout = 6
       iundac = 13
       iunerr = 0
       pdbline = 0
       pdbrec = ' '
       type_code = 'l'
       call splitstr(num_aa,
     *   "ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,"//
     *   "MET,PHE,PRO,SER,THR,TRP,TYR,VAL,ASX,GLX,UNK",
     *   aa_list,23,',')
       call splitstr(num_na,
     *   "A  ,C  ,G  ,T  ,U  ,I  ,+A ,+C ,+G ,+T ,+U ,+I ",
     *   na_list,12,',')
       call hash_init(mapstr,mapchain,
     *   NUMSTR,nmap,mhash,NUMHSH)
       result = init_(iuninp,iunout,iundac,iunerr)
       line_ = 2048
       call getfiles
       if (ckdict(1:lckdict).ne.' ')
     *   result = dict_(ckdict(1:lckdict),'valid dtype catck')
       if (inpcif(1:linpcif).eq.'-' .or. 
     *   inpcif(1:linpcif) .eq. ' ') then
         iuninp = 5
         inpcif = ' '
         result = init_(iuninp,iunout,iundac,iunerr)
       endif
       if (outent(1:loutent).eq.'-') outent = ' '
       result = ocif_(inpcif(1:linpcif))
       if(.not.result) then
         call c2perr(
     *   ' failed to open '//inpcif(1:linpcif))
       endif
       if (outent(1:1) .ne. ' ')
     * open(unit=iunpdb,file=outent(1:loutent),status='unknown')
       result = data_(' ')
       if (.not.result) then
         call c2perr(' failed to open block')
       endif
       pdbent = bloc_(1:6)
         call c2pmsg(' ',' Data Block Name '//pdbent)
       numRemark = 0
       numFtnote = 0
       numHet = 0
       numHelix = 0
       numSheet = 0
       numTurn = 0
       numSite = 0
       numXform = 0
       numCoord = 0
       numTer = 0
       numConect = 0
       numSeq = 0
       numTitle = 0
       numCompound = 0
       numSource = 0
       numCaveat = 0 
C      set up formats for simple string-based record       
       sfmt_orig = '(a6,4x,a,a)'
       sfmt_cont = '(a6,1x,i3,a)'
       if (wide_code .eq. "yes") then
         sfmt_orig = '(a,t23,a,a)'
         sfmt_cont = '(a,t19,i3,1x,a)'
       endif

       call get_hets
       call get_chains
       call proc_header
       call proc_title
       call proc_keywds
       call proc_expdta
       call proc_author
       call proc_remark
       result = findtag_('_entity_poly_seq.entity_id')
       call proc_seqres
       result = findtag_('_atom_sites_footnote.id')
       if (result) call proc_ftnote
       if (.not.result)
     *    result = findtag_('_cell.length_a')
       call proc_cryst1
       call proc_origx
       call proc_scale
       call proc_atom
       call proc_master
       stop
       end
      subroutine getfiles
C
      include 'cif2pdb.cmn'
      external pdbint
      external pdbreal
      external pdbstr
      external typeset
      logical dict_
      character*256 temp,temp2,cline
      integer iargc
      logical backarg,ffile
      myentry = '*???'
      numarg = iargc()
      call getenv("CIF2PDB_INPUT_CIF",inpcif)
      call getenv("CIF2PDB_OUTPUT_ENTRY",outent)
      call getenv("CIF2PDB_CHECK_DICTIONARY",ckdict)
      karg = 0
      kfarg = 0
      iwant = 0
      ifound = 0
      isi = 0
      iso = 0       
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
            temp = cstr(kfarg)(1:256)
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
        if (iwant.eq.2) outent = temp(1:ll)
        if (iwant.eq.3) then
          ckdict = temp(1:ll)
          lckdict = max(1,nblen(ckdict))
          result = dict_(ckdict(1:lckdict),"valid dtype catck")
          ckdict = " "
          lckdict = 1
        endif
        if (iwant.eq.4) myentry = tbxxupcs(temp(1:ll))
        if (iwant.eq.5) then
          open(unit=iuninp,file=temp(1:ll),status='OLD',err=900)
          kfarg = 1
          nfarg = 0
          ffile = .true.
        endif
        if (iwant.eq.6) then
          kpmap = nmap
          temp = nounder(temp(1:ll))
          ll = nblen(temp)
          call hash_store(temp(1:ll)//char(0),mapstr,mapchain,
     *      NUMSTR,nmap,mhash,NUMHSH,ifrom)
          if (ifrom.eq.0)
     *      call c2perr(' More than NUMSTR strings mapped ')
          if (ifrom.ne.kpmap+1)
     *      call c2pwarn(' Duplicate mapping of '//temp(1:ll))
          iwant = 7
          go to 100
        endif
        if (iwant.eq.7) then
          kpmap = nmap
          temp = nounder(temp(1:ll))
          ll = nblen(temp)
          call hash_store(temp(1:ll)//char(0),mapstr,mapchain,
     *      NUMSTR,nmap,mhash,NUMHSH,ito)
          if (ito.eq.0) 
     *      call c2perr(' More than NUMSTR strings mapped ')
          if (ito.eq.kpmap+1) mapto(ito) = 0
          mapto(ifrom) = ito
        endif
        if (iwant.eq.8) then
          if(temp(1:ll).eq."u" .or.
     *      temp(1:ll).eq."l" .or.
     *      temp(1:ll).eq."p" ) then
            type_code = temp(1:1)
          else
            go to 900
          endif
        endif
        if (iwant.eq.9) then
          if(temp(1:ll).eq."yes" .or.
     *      temp(1:ll).eq."no" ) then
            wide_code = temp(1:ll)
          else
            go to 900
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
          if (temp(2:2).eq."p") iwant = 4
          if (temp(2:2).eq."f") iwant = 5
          if (temp(2:2).eq."m") iwant = 6
          if (temp(2:2).eq."t") iwant = 8
          if (temp(2:2).eq."w") iwant = 9
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
              outent = temp(1:ll)
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
 500  linpcif = max(1,nblen(inpcif))
      loutent = max(1,nblen(outent))
      lckdict = max(1,nblen(ckdict))
      return
 900  write(iunerr,'(a)')
     *  ' cif2pdb [-i input_cif] [-o output_entry] [-d dictionary]',
     *  '         [-p pdb_entry_id] [-f command_file] [-t u|l|p]',
     *  '         [-w yes|no] [-m string_in_cif string_in_pdb]',
     *  '         [[[input_cif] [[output_entry] [dictionary]]]',
     *  ' input_cif defaults to $CIF2PDB_INPUT_CIF or stdin',
     *  ' output_cif defaults to $CIF2PDB_OUTPUT_ENTRY or stdout',
     *  ' dictionary defaults to $CIF2PDB_CHECK_DICTIONARY',
     *  ' multiple dictionaries may be specified ',
     *  ' input_cif of "-" is stdin, output_entry of "-" is stdout',
     *  ' -t has values of u for upper case, l for upper/lower,',
     *  '    p for PDB typesetting codes, (default -t l)',
     *  ' -w has the values yes or no (default -w no)'
       stop
       end

       subroutine get_hets
       include 'cif2pdb.cmn'
       external pdbreal
       external pdbint
       external typeset
       character*9 curcat
       character*10 yesno
       logical findtag_
       logical result_blk
       logical result_id
       logical result_mon_nstd_flag
       curcat = 'chem_comp'
       numhetl = 0
       result_blk = findtag_('_chem_comp.id')
       if (.not.result_blk) return
       loop_ = .false.
 100   result_id = pdbstr(curcat,'id',3,hetl(numhetl+1))
       if(result_id) then
         result_mon_nstd_flag = pdbstr(curcat,'mon_nstd_flag',10,yesno)
         if(result_mon_nstd_flag.and.
     *     (yesno(1:1).ne.'n' .and. 
     *      yesno(1:1).ne.'N' .and.
     *      yesno(1:1).ne.' ')) then
           numhetl = numhetl+1
         else
           if(.not.result_mon_nstd_flag .or.
     *        (result_mon_nstd_flag.and.yesno(1:1).eq.' ')) then
             do ii = 1,num_aa
               if (hetl(numhetl+1) .eq. aa_list(ii) ) go to 120
             enddo
             do ii = 1,num_na
               if (hetl(numhetl+1) .eq. na_list(ii) ) go to 120
             enddo
             numhetl = numhetl+1
           endif
         endif
       endif
 120   if (loop_) go to 100
       return
       end
       subroutine get_chains
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       character*33 entities(127)
       character*12 entity_type(127)
       logical result_id
       logical result_entity_id
       logical result_type
       logical result_entity
       logical ispoly
       num_ent = 0
       num_chains = 0
       loop_ = .false.
 100   result_entity = pdbstr('entity','id',33,
     *   entities(num_ent+1))
       result_type = pdbstr('entity','type',12,
     *   entity_type(num_ent+1))
       if(result_entity .and. result_type) then
         num_ent = num_ent+1
       endif
       if (loop_) go to 100
 200   result_entity_id = pdbstr('struct_asym','entity_id',33,
     *   entity_list(num_chains+1))
       if (result_entity_id) then
         ispoly = .true.
         do ii = 1,num_ent
         if (entity_list(num_chains+1).eq.
     *     entities(ii)) then
           if (entity_type(ii)(1:1).ne.'P' .and.
     *       entity_type(ii)(1:1).ne.'p') ispoly = .false.
           go to 300
         endif
         enddo
 300     if (ispoly) then
           result_id = pdbstr('struct_asym','id',10,
     *       chain_list(num_chains+1))
           if(result_id) then
             num_chains = num_chains+1
           endif
         endif
       endif
       if (loop_) go to 200
       return
       end
       subroutine proc_atom
C
C      process and atom list
C
       include 'cif2pdb.cmn'
       external pdbint
       external typeset
       character*9 curcat
       character*19 anisocat
       logical result_group_PDB,result_id,
     *    result_label_atom_id,
     *    result_label_alt_id,result_label_comp_id,
     *    result_label_asym_id,result_label_seq_id,
     *    result_auth_asym_id,
     *    result_cartn_x,result_cartn_y,result_cartn_z,
     *    result_fract_x,result_fract_y,result_fract_z,
     *    result_cartn_x_esd,result_cartn_y_esd,result_cartn_z_esd,
     *    result_fract_x_esd,result_fract_y_esd,result_fract_z_esd,
     *    result_occupancy,result_B_iso_or_equiv,
     *    result_occupancy_esd,result_B_iso_or_equiv_esd,
     *    result_footnote_id,result_type_symbol,
     *    result_label_model_id,
     *    result_aniso_U(6),
     *    result_aniso_Usd(6),
     *    doanisou, dosiguij, doanisotrop, domodel
       common/proc_atom_common/result_group_PDB,result_id,
     *    result_label_atom_id,
     *    result_label_alt_id,result_label_comp_id,
     *    result_label_asym_id,result_label_seq_id,
     *    result_auth_asym_id,
     *    result_cartn_x,result_cartn_y,result_cartn_z,
     *    result_fract_x,result_fract_y,result_fract_z,
     *    result_cartn_x_esd,result_cartn_y_esd,result_cartn_z_esd,
     *    result_fract_x_esd,result_fract_y_esd,result_fract_z_esd,
     *    result_occupancy,result_B_iso_or_equiv,
     *    result_occupancy_esd,result_B_iso_or_equiv_esd,
     *    result_footnote_id,result_type_symbol,
     *    result_label_model_id,result_bmk
       character*6 group_PDB
       character*9 idstr,oidstr,anisoidstr
       character*10 label_atom_id
       character*3 label_alt_id
       character*10 label_comp_id
       character*10 label_asym_id
       character*4 auth_asym_id
       character*10 label_seq_id
       character*13 xstr,ystr,zstr
       character*13 xesdstr,yesdstr,zesdstr
       character*9 modstr
       character*6 occstr,Bstr
       character*6 occesdstr,Besdstr
       character*6 fidstr
       character*13 anistr(6)
       character*13 anisdstr(6)
       integer anisobk, iapass, lad, alistbk,ad, oad, nad
       integer kapass
       logical result_anisotrop
       logical result_bmk
       logical findtag_
       logical bkmrk_
       logical anisosorted
       
       integer iwcoord, iwasn, iwaname, iwaltid, iwrname
       integer iwcid, iwrnum, iwfoot, iwmodel
       
       integer maxmodel, minmodel, modelmode
       
       integer iasn
       
       real*8 cartn_x,cartn_y,cartn_z,occupancy,B_iso_or_equiv
       real*8 cartn_x_esd,cartn_y_esd,cartn_z_esd
       real*8 occupancy_esd,B_iso_or_equiv_esd
       real*8 fract_xyz(3),fract_xyz_esd(3)
       real*8 aniso(6),anisd(6)
       logical fcoord, fcoord_esd
       integer ianiso(6),ianisd(6)
       integer footnote_id,modnum,pmodnum
       character*4 type_symbol
       character*10 blanks
       integer ic
       data blanks/'          '/
       data pi/3.1415926536/
       curcat = 'atom_site'
       anisocat = 'atom_site_anisotrop'
       pdbrec = 'ATOM'
       pmodnum = -999999
       text_= .false.
       fcoord = .false.
       fcoord_esd = .false.
       
       anisosorted=.false.
       iapass = 0
       iwcoord = 8
       iwasn = 5
       iwaname = 4
       iwaltid = 1
       iwrname = 3
       iwcid = 3
       iwrnum = 5
       iwfoot = 3
       iwmodel = 4
       if (wide_code.eq."yes") then
         iwcoord = 13
         iwasn = 9
         iwaname = 10
         iwaltid = 3
         iwrname = 10
         iwcid = 10
         iwrnum = 10
         iwfoot = 6
         iwmodel = 9
       endif
C
C      Scan to determine the range of models used
C
       maxmodel = -1000000000
       minmodel = 1000000000
       modelmode = 0
       domodel = .false.
       
       result_label_model_id = 
     *   findtag_('_atom_site.label_model_id')
       if (result_label_model_id) then
         modelmode = 1
         domodel = .true.
       else
         result_label_model_id = 
     *     findtag_('_atom_site.pdb2cif_label_model_id')
         if (result_label_model_id) then
           modelmode = 2
           domodel = .true.
         else
           result_label_model_id = 
     *       findtag_('_atom_site.pdbx_PDB_model_num')
           if (result_label_model_id) then
             modelmode = 3
             domodel = .true.
           endif
         endif
       endif
       if (domodel) then
       loop_ = .false.
 60    if (modelmode.eq.1) then
         result_label_model_id = pdbint(curcat,
     *     'label_model_id',
     *     iwmodel,modstr,modnum)
       else
         if (modelmode.eq.2) then
            result_label_model_id = pdbint(curcat,
     *       'pdb2cif_label_model_id',
     *       iwmodel,modstr,modnum)
         else
            result_label_model_id = pdbint(curcat,
     *        'pdbx_PDB_model_num',
     *        iwmodel,modstr,modnum)
         endif
       endif
       if (result_label_model_id .and. modstr .ne. ' ') then
         if (modnum.lt.minmodel) minmodel = modnum
         if (modnum.gt.maxmodel) maxmodel = modnum
       endif
       if (loop_) go to 60
       if (maxmodel.lt.minmodel.or.
     *  (maxmodel.eq.1 .and. minmodel.eq.1)) domodel = .false.
       endif

C
C      Prepare to deal with separate atom_site_anisotrop loop
C
       result_anisotrop = findtag_('_atom_site_anisotrop.id')
       anisobk = 0
       result_bmk = .false.
       doanisotrop = .false.
       if (result_anisotrop) then
         result_bmk = bkmrk_(anisobk)
         iapass = 0
         doanisotrop = .true.
         anisosorted = .true.
         result_anisotrop = pdbint(anisocat,'id',
     *     iwasn,anisoidstr,nad)
 80      oad = nad
         if (.not.loop_) go to 90
           result_anisotrop = pdbint(anisocat,'id',
     *       iwasn,anisoidstr,nad)
         if (nad.gt.oad) go to 80
         anisosorted = .false.
        endif
 90    iasn = 0
       ad = -1
       result_id = findtag_('_atom_site.id')
 100   pdbline = pdbline+1
       iasn = iasn+1
       result_group_PDB = pdbstr(curcat,'group_PDB',6,group_PDB)
       if (.not.result_group_PDB) then
         group_PDB = 'ATOM  '
       endif
       result_id = pdbint(curcat,'id',iwasn,idstr,id)
       if (doanisotrop) then
         alistbk = 0
         result_bmk = bkmrk_(alistbk)
       endif
       oidstr = idstr
       if (wide_code.eq."yes") then
         call rjust(idstr)
         if (.not.result_id.or.id.eq.0) then
           id = mod(iasn,1000000000)
             write (idstr,'(i9)') id
       endif
       else
         call rjust(idstr(1:5))
         if (.not.result_id.or.id.eq.0) then
           id = mod(iasn,100000)
           write (idstr,'(i5)') id
         endif
       endif
       result_label_atom_id = pdbstr(curcat,'label_atom_id',
     *   iwaname,label_atom_id)
       result_label_alt_id = pdbstr(curcat,'label_alt_id',
     *   iwaltid,label_alt_id)
       result_label_comp_id = pdbstr(curcat,'label_comp_id',
     *   iwrname,label_comp_id)
       result_label_asym_id = pdbstr(curcat,'label_asym_id',
     *   iwcid,label_asym_id)
       result_auth_asym_id = pdbstr(curcat,'auth_asym_id',4,
     *   auth_asym_id)
       ll=nblen(auth_asym_id)
       tmparg(1:max(1,ll))=auth_asym_id
       if(ll.eq.0) then
         auth_asym_id = ' '
       else
         auth_asym_id = tmparg(1:ll)
       endif
C
C      Get sequence number and insertion code.  If the result
C      is numeric, then we treat insertion code as blank
C      As of version 0.8.0 of mmCIF, we need to check for
C      and author variant, first
C
       result_label_seq_id  = pdbstr(curcat,'auth_seq_id',iwrnum,
     *   label_seq_id)
       if(.not.result_label_seq_id.or.label_seq_id.eq.'     ')
     * result_label_seq_id  = pdbstr(curcat,'label_seq_id',iwrnum,
     *   label_seq_id)
       call ctonum
       ll=nblen(label_seq_id)
       tmparg(1:max(1,ll))=label_seq_id(1:max(1,ll))
       if(type_.eq.'numb') then
         if(ll.lt.iwrnum-1) then
           label_seq_id = blanks(1:iwrnum-1-ll)//tmparg(1:ll)//' '
         else
           label_seq_id = tmparg(1:ll)//' '
         endif
       else
         if(ll.lt.iwrnum) then
           label_seq_id = blanks(1:iwrnum-ll)//tmparg(1:ll)
         endif
       endif
       if (.not.fcoord) then
       result_cartn_x = 
     *   pdbreal(curcat,'cartn_x',iwcoord,xstr,cartn_x)
       result_cartn_y = 
     *   pdbreal(curcat,'cartn_y',iwcoord,ystr,cartn_y)
       result_cartn_z =
     *   pdbreal(curcat,'cartn_z',iwcoord,zstr,cartn_z)
       result_cartn_x_esd = 
     *   pdbreal(curcat,'cartn_x_esd',iwcoord,xesdstr,
     *   cartn_x_esd)
       result_cartn_y_esd = 
     *   pdbreal(curcat,'cartn_y_esd',iwcoord,yesdstr,
     *   cartn_y_esd)
       result_cartn_z_esd = 
     *   pdbreal(curcat,'cartn_z_esd',iwcoord,zesdstr,
     *   cartn_z_esd)
       fcoord = .not.(result_cartn_x.and.result_cartn_y
     *   .and.result_cartn_z)
       fcoord_esd = (result_cartn_x_esd.and.result_cartn_y_esd
     *   .and.result_cartn_z_esd)
       endif
       if (fcoord) then
       result_fract_x = pdbreal(curcat,'fract_x',13,xstr,fract_xyz(1))
       result_fract_y = pdbreal(curcat,'fract_y',13,ystr,fract_xyz(2))
       result_fract_z = pdbreal(curcat,'fract_z',13,zstr,fract_xyz(3))
       result_fract_x_esd = 
     *   pdbreal(curcat,'fract_x_esd',13,xesdstr,fract_xyz_esd(1))
       result_fract_y_esd = 
     *   pdbreal(curcat,'fract_y_esd',13,yesdstr,fract_xyz_esd(2))
       result_fract_z_esd = 
     *   pdbreal(curcat,'fract_z_esd',13,zesdstr,fract_xyz_esd(3))
       fcoord = (result_fract_x.and.result_fract_y
     *   .and.result_fract_z)
       fcoord_esd =(result_fract_x_esd.and.result_fract_y_esd
     *   .and.result_fract_z_esd)
       endif
       if (fcoord) then
         cartn_x = 0.
         cartn_y = 0.
         cartn_z = 0.
         cartn_x_esd = 0.
         cartn_y_esd = 0.
         cartn_z_esd = 0.
         do jj = 1,3
           if (.not.fcoord_esd) fract_xyz_esd(jj)=0.
           cartn_x = cartn_x + matf2o(1,jj)*fract_xyz(jj)
           cartn_y = cartn_y + matf2o(2,jj)*fract_xyz(jj)
           cartn_z = cartn_z + matf2o(3,jj)*fract_xyz(jj)
           cartn_x_esd = cartn_x_esd +
     *       abs(matf2o(1,jj)*fract_xyz_esd(jj))
           cartn_y_esd = cartn_y_esd +
     *       abs(matf2o(2,jj)*fract_xyz_esd(jj))
           cartn_z_esd = cartn_z_esd +
     *       abs(matf2o(3,jj)*fract_xyz_esd(jj))
         enddo
         cartn_x = cartn_x + vecf2o(1)
         cartn_y = cartn_y + vecf2o(2)
         cartn_z = cartn_z + vecf2o(3)
         if (wide_code .eq. "yes") then
           write(xstr,'(f13.5)')cartn_x
           write(ystr,'(f13.5)')cartn_y
           write(zstr,'(f13.5)')cartn_z
           write(xesdstr,'(f13.5)')cartn_x_esd
           write(yesdstr,'(f13.5)')cartn_y_esd
           write(zesdstr,'(f13.5)')cartn_z_esd
         else
           write(xstr,'(f8.3)')cartn_x
           write(ystr,'(f8.3)')cartn_y
           write(zstr,'(f8.3)')cartn_z
           write(xesdstr,'(f8.3)')cartn_x_esd
           write(yesdstr,'(f8.3)')cartn_y_esd
           write(zesdstr,'(f8.3)')cartn_z_esd
         endif
       endif
       if (wide_code .eq. "yes") then
         call fixdec(xstr,5)
         call fixdec(ystr,5)
         call fixdec(zstr,5)
         if (fcoord_esd) then
           call fixdec(xesdstr,5)
           call fixdec(yesdstr,5)
           call fixdec(zesdstr,5)
         else
           xesdstr = ' '
           yesdstr = ' '
           zesdstr = ' '
         endif
       else
         call fixdec(xstr,3)
         call fixdec(ystr,3)
         call fixdec(zstr,3)
         if (fcoord_esd) then
           call fixdec(xesdstr,3)
             call fixdec(yesdstr,3)
           call fixdec(zesdstr,3)
         else
           xesdstr = ' '
           yesdstr = ' '
           zesdstr = ' '
         endif
       endif
       result_occupancy = pdbreal(curcat,'occupancy',6,occstr,
     *   occupancy)
       result_B_iso_or_equiv = pdbreal(curcat,'B_iso_or_equiv',6,
     *   Bstr,B_iso_or_equiv)
       result_occupancy_esd = pdbreal(curcat,'occupancy_esd',6,
     *   occesdstr,occupancy_esd)
       result_B_iso_or_equiv_esd = 
     *   pdbreal(curcat,'B_iso_or_equiv_esd',6,
     *   Besdstr,B_iso_or_equiv_esd)
       if (.not.result_B_iso_or_equiv) then
         result_B_iso_or_equiv = pdbreal(curcat,'U_iso_or_equiv',6,
     *     Bstr,B_iso_or_equiv)
         B_iso_or_equiv = B_iso_or_equiv*8.*pi**2
         result_B_iso_or_equiv_esd = pdbreal(curcat,
     *    'U_iso_or_equiv_esd',6,
     *     Besdstr,B_iso_or_equiv_esd)
         if (result_B_iso_or_equiv_esd) then
           B_iso_or_equiv_esd = B_iso_or_equiv_esd*8.*pi**2
           write(Besdstr,'(f6.2)')B_iso_or_equiv_esd
         else
           B_iso_or_equiv_esd = 0.
           Besdstr = ' '
         endif
       endif
       if (result_occupancy_esd) then
         write(occesdstr,'(f6.2)') occupancy_esd
       else
         occesdstr = ' '
       endif
       result_footnote_id = pdbint(curcat,'footnote_id',iwfoot,
     *   fidstr,footnote_id)
       if(fidstr.ne.' ' .and. (.not.result_footnote_id))
     *   call c2pwarn(' Non-numeric _atom_site.footnote_id '//
     *   fidstr)
       do ii = 1,6
         aniso(ii)=0.
         anisd(ii)=0.
       enddo
       doanisou = .false.
       dosiguij = .false.
       if (doanisotrop) then
         result_bmk = bkmrk_(anisobk)
         kapass = iapass
         lad = ad
 110     result_anisotrop = pdbint(anisocat,'id',
     *     iwasn,anisoidstr,ad)
         if (id.eq.ad) go to 130
         if (anisosorted .and. ad.gt.id) go to 120
         if (iapass.eq.kapass+1.and.ad.eq.lad) go to 120
              if (loop_) go to 110
         iapass = iapass+1
         if (iapass.le.kapass+1) go to 110
 120     continue
         anisobk = 0
         result_bmk = bkmrk_(anisobk)
         result_bmk = bkmrk_(alistbk)
         result_group_PDB = pdbstr(curcat,'group_PDB',6,group_PDB)
         if (.not.result_group_PDB) then
           group_PDB = 'ATOM  '
         endif
         result_id = pdbint(curcat,'id',iwasn,idstr,id)
         go to 140
 130     continue
         result_aniso_U(1) = pdbreal(anisocat,'U[1][1]',13,
     *     anistr(1),aniso(1))
         result_aniso_U(2) = pdbreal(anisocat,'U[2][2]',13,
     *     anistr(2),aniso(2))
         result_aniso_U(3) = pdbreal(anisocat,'U[3][3]',13,
     *     anistr(3),aniso(3))
         result_aniso_U(4) = pdbreal(anisocat,'U[1][2]',13,
     *     anistr(4),aniso(4))
         result_aniso_U(5) = pdbreal(anisocat,'U[1][3]',13,
     *     anistr(5),aniso(5))
         result_aniso_U(6) = pdbreal(anisocat,'U[2][3]',13,
     *     anistr(6),aniso(6))
         result_aniso_Usd(1) = pdbreal(anisocat,'U[1][1]_esd',13,
     *     anisdstr(1),anisd(1))
         result_aniso_Usd(2) = pdbreal(anisocat,'U[2][2]_esd',13,
     *     anisdstr(2),anisd(2))
         result_aniso_Usd(3) = pdbreal(anisocat,'U[3][3]_esd',13,
     *     anisdstr(3),anisd(3))
         result_aniso_Usd(4) = pdbreal(anisocat,'U[1][2]_esd',13,
     *     anisdstr(4),anisd(4))
         result_aniso_Usd(5) = pdbreal(anisocat,'U[1][3]_esd',13,
     *     anisdstr(5),anisd(5))
         result_aniso_Usd(6) = pdbreal(anisocat,'U[2][3]_esd',13,
     *     anisdstr(6),anisd(6))
         do ii = 1, 6
           if (result_aniso_U(ii) .and. aniso(ii).ne.0.) 
     *       doanisou = .true.
           if (result_aniso_Usd(ii) .and. anisd(ii).ne.0.) 
     *       dosiguij = .true.
         enddo
         if (.not.doanisou .and. .not.dosiguij) then
           result_aniso_U(1) = pdbreal(anisocat,'B[1][1]',13,
     *       anistr(1),aniso(1))
           result_aniso_U(2) = pdbreal(anisocat,'B[2][2]',13,
     *       anistr(2),aniso(2))
           result_aniso_U(3) = pdbreal(anisocat,'B[3][3]',13,
     *       anistr(3),aniso(3))
           result_aniso_U(4) = pdbreal(anisocat,'B[1][2]',13,
     *       anistr(4),aniso(4))
           result_aniso_U(5) = pdbreal(anisocat,'B[1][3]',13,
     *       anistr(5),aniso(5))
           result_aniso_U(6) = pdbreal(anisocat,'B[2][3]',13,
     *       anistr(6),aniso(6))
           result_aniso_Usd(1) = pdbreal(anisocat,'B[1][1]_esd',13,
     *       anisdstr(1),anisd(1))
           result_aniso_Usd(2) = pdbreal(anisocat,'B[2][2]_esd',13,
     *       anisdstr(2),anisd(2))
           result_aniso_Usd(3) = pdbreal(anisocat,'B[3][3]_esd',13,
     *       anisdstr(3),anisd(3))
           result_aniso_Usd(4) = pdbreal(anisocat,'B[1][2]_esd',13,
     *       anisdstr(4),anisd(4))
           result_aniso_Usd(5) = pdbreal(anisocat,'B[1][3]_esd',13,
     *       anisdstr(5),anisd(5))
           result_aniso_Usd(6) = pdbreal(anisocat,'B[2][3]_esd',13,
     *       anisdstr(6),anisd(6))
           do ii = 1, 6
             aniso(ii) = aniso(ii)/(8.*pi**2)
             anisd(ii) = anisd(ii)/(8.*pi**2)
             if (result_aniso_U(ii) .and. aniso(ii).ne.0.) 
     *         doanisou = .true.
             if (result_aniso_Usd(ii) .and. anisd(ii).ne.0.) 
     *         dosiguij = .true.
           enddo
         endif
         go to 120
       else
         result_aniso_U(1) = pdbreal(curcat,'aniso_U[1][1]',13,
     *     anistr(1),aniso(1))
         result_aniso_U(2) = pdbreal(curcat,'aniso_U[2][2]',13,
     *     anistr(2),aniso(2))
         result_aniso_U(3) = pdbreal(curcat,'aniso_U[3][3]',13,
     *     anistr(3),aniso(3))
         result_aniso_U(4) = pdbreal(curcat,'aniso_U[1][2]',13,
     *     anistr(4),aniso(4))
         result_aniso_U(5) = pdbreal(curcat,'aniso_U[1][3]',13,
     *     anistr(5),aniso(5))
         result_aniso_U(6) = pdbreal(curcat,'aniso_U[2][3]',13,
     *     anistr(6),aniso(6))
         result_aniso_Usd(1) = pdbreal(curcat,'aniso_U[1][1]_esd',13,
     *     anisdstr(1),anisd(1))
         result_aniso_Usd(2) = pdbreal(curcat,'aniso_U[2][2]_esd',13,
     *     anisdstr(2),anisd(2))
         result_aniso_Usd(3) = pdbreal(curcat,'aniso_U[3][3]_esd',13,
     *     anisdstr(3),anisd(3))
         result_aniso_Usd(4) = pdbreal(curcat,'aniso_U[1][2]_esd',13,
     *     anisdstr(4),anisd(4))
         result_aniso_Usd(5) = pdbreal(curcat,'aniso_U[1][3]_esd',13,
     *     anisdstr(5),anisd(5))
         result_aniso_Usd(6) = pdbreal(curcat,'aniso_U[2][3]_esd',13,
     *     anisdstr(6),anisd(6))
         do ii = 1, 6
           if (result_aniso_U(ii) .and. aniso(ii).ne.0.) 
     *       doanisou = .true.
           if (result_aniso_Usd(ii) .and. anisd(ii).ne.0.) 
     *       dosiguij = .true.
         enddo
         if (.not.doanisou .and. .not.dosiguij) then
           result_aniso_U(1) = pdbreal(curcat,'aniso_B[1][1]',13,
     *       anistr(1),aniso(1))
           result_aniso_U(2) = pdbreal(curcat,'aniso_B[2][2]',13,
     *       anistr(2),aniso(2))
           result_aniso_U(3) = pdbreal(curcat,'aniso_B[3][3]',13,
     *       anistr(3),aniso(3))
           result_aniso_U(4) = pdbreal(curcat,'aniso_B[1][2]',13,
     *       anistr(4),aniso(4))
           result_aniso_U(5) = pdbreal(curcat,'aniso_B[1][3]',13,
     *       anistr(5),aniso(5))
           result_aniso_U(6) = pdbreal(curcat,'aniso_B[2][3]',13,
     *       anistr(6),aniso(6))
           result_aniso_Usd(1) = pdbreal(curcat,'aniso_B[1][1]_esd',13,
     *       anisdstr(1),anisd(1))
           result_aniso_Usd(2) = pdbreal(curcat,'aniso_B[2][2]_esd',13,
     *       anisdstr(2),anisd(2))
           result_aniso_Usd(3) = pdbreal(curcat,'aniso_B[3][3]_esd',13,
     *       anisdstr(3),anisd(3))
           result_aniso_Usd(4) = pdbreal(curcat,'aniso_B[1][2]_esd',13,
     *       anisdstr(4),anisd(4))
           result_aniso_Usd(5) = pdbreal(curcat,'aniso_B[1][3]_esd',13,
     *       anisdstr(5),anisd(5))
           result_aniso_Usd(6) = pdbreal(curcat,'aniso_B[2][3]_esd',13,
     *       anisdstr(6),anisd(6))
           do ii = 1, 6
             aniso(ii) = aniso(ii)/(8.*pi**2)
             anisd(ii) = anisd(ii)/(8.*pi**2)
             if (result_aniso_U(ii) .and. aniso(ii).ne.0.) 
     *         doanisou = .true.
             if (result_aniso_Usd(ii) .and. anisd(ii).ne.0.) 
     *         dosiguij = .true.
           enddo

         endif
       endif
 140   result_type_symbol = pdbstr(curcat,'type_symbol',4,
     *  type_symbol)
C
C     Set up MODEL/ENDMDL brackets for NMR entries
C
       if (domodel) then
       if (modelmode.eq.1) then
         result_label_model_id = pdbint(curcat,
     *     'label_model_id',
     *     iwmodel,modstr,modnum)
       else
         if (modelmode.eq.2) then
            result_label_model_id = pdbint(curcat,
     *       'pdb2cif_label_model_id',
     *       iwmodel,modstr,modnum)
         else
            result_label_model_id = pdbint(curcat,
     *        'pdbx_PDB_model_num',
     *        iwmodel,modstr,modnum)
         endif
       endif

       if (.not. result_label_model_id) then
         modstr = '    '
         modnum = -999999
       endif
       if (modstr .eq. '    ') modnum = -999999
       if (pmodnum .ne. modnum) then
         if (pmodnum.ne.-999999) write(iunpdb,'(6HENDMDL)')
         pmodnum = modnum
         if (modnum.ne.-999999) then 
           if (wide_code .eq."yes") then
             write(iunpdb,'(6HMODEL ,2x,i9)') modnum
           else
             write(iunpdb,'(6HMODEL ,4x,i4)') modnum
           endif
         endif
       endif
       endif
         
C
C     Fix up the atom name, if it begins with the atom type
C
      lat = nblen(type_symbol)
C
C     The type_symbol may have more structure than simply
C     an element symbol.  There may be very general types,
C     by the most likely are element symbol followed by
C     a charge
C
      if (lat.gt.1) then
        if(type_symbol(lat:lat).eq."+" .or.
     *    type_symbol(lat:lat).eq."-") then
          if (index("0123456789",type_symbol(lat-1:lat-1))
     *      .gt.0) then
            lat = lat-2
          else
            lat = lat-1
          endif
        endif
      endif
C
      if (label_atom_id.eq.' ') label_atom_id = oidstr 
      laid = nblen(label_atom_id)
      tmparg = label_atom_id(1:max(1,laid))
      ic = ichar(label_atom_id(1:1)) - ichar('0')
      if ( ic.lt.0 .or. ic.gt.9) then
        if(laid.gt.0.and.laid.lt.iwaname.and.lat.gt.0) then
          if(label_atom_id(1:lat).eq.type_symbol(1:lat)
     *      .and. lat.eq.1) then
            label_atom_id = ' '//tmparg(1:max(1,laid))
            if (label_atom_id.eq." O''") then
              label_atom_id = ' OXT'
            endif
          endif
        else
          if(laid.gt.0.and.laid.lt.iwaname) then
            label_atom_id = ' '//tmparg(1:max(1,laid))
          endif
          if(laid.eq.iwaname.and.lat.eq.1.and.iwaname.eq.4) then
            ic = ichar(label_atom_id(4:4)) - ichar('0')
            if (ic.ge.0 .and. ic.le.9) then
              label_atom_id = tmparg(4:4)//tmparg(1:3)
            endif
          endif
        endif
      endif
      tmparg = type_symbol(1:max(1,lat))
      if (lat.eq.1) type_symbol = " "//tmparg(1:1)
C
C     Fix up group_PDB
C
      if(.not.result_group_PDB.or.group_PDB .eq. " ") then
        group_PDB = 'ATOM'
        do ii = 1,numhetl
        if (label_comp_id.eq.hetl(ii)) group_PDB = 'HETATM'
        enddo
      else
      if (group_PDB(1:1).eq.'A' .or. group_PDB(1:1).eq.'a')
     *  group_PDB = 'ATOM'
      if (group_PDB(1:1).eq.'H' .or. group_PDB(1:1).eq.'h')
     *  group_PDB = 'HETATM'
      if (group_PDB(1:1).eq.'T' .or. group_PDB(1:1).eq.'t')
     *  group_PDB = 'TER'
      endif
      if (group_PDB(1:1).eq.'A' .or.group_PDB(1:1).eq.'H')
     *  numCoord = numCoord+1
      if (group_PDB(1:1).eq.'T') numTer = numTer+1
      if (wide_code.eq."yes") then
        call rjust(label_asym_id)
        call rjust(label_comp_id)
        call rjust(fidstr)
        write(iunpdb,
     * '(a6,2x,i9,1x,3x,1x,a10,a3,3a10,3a13,2f6.2,a6,2(1x,a4))')
     * group_PDB,id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,xstr,ystr,zstr,
     * occupancy,B_iso_or_equiv,fidstr,auth_asym_id,type_symbol
      else
        call rjust(label_comp_id(1:3))
        call rjust(fidstr(1:3))
        write(iunpdb,
     * '(a6,i5,1x,a4,a1,a3,1x,a1,a5,3x,3f8.3,2f6.2,1x,a3,2x,2a4)')
     * group_PDB,id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,cartn_x,cartn_y,cartn_z,
     * occupancy,B_iso_or_equiv,fidstr,auth_asym_id,
     * type_symbol
      endif
C
C     Process SIGATM
C
      if (fcoord_esd.or.result_B_iso_or_equiv_esd.or.
     *  result_occupancy_esd) then
        if (wide_code.eq."yes") then
          call rjust(label_asym_id)
          call rjust(label_comp_id)
          call rjust(fidstr)
          write(iunpdb,
     *     '(a6,2x,i9,1x,3x,1x,a10,a3,3a10,3a13,2a6,a6,2(1x,a4))')
     *     'SIGATM',id,label_atom_id,label_alt_id,label_comp_id,
     *     label_asym_id,label_seq_id,
     *     xesdstr,yesdstr,zesdstr,
     *     occesdstr,Besdstr,fidstr,auth_asym_id,type_symbol
        else
           call rjust(label_comp_id(1:3))
           call rjust(fidstr(1:3))
           write(iunpdb,
     *     '(a6,i5,1x,a4,a1,a3,1x,a1,a5,3x,3a8,2a6,1x,a3,2x,2a4)')
     *     'SIGATM',id,label_atom_id,label_alt_id,label_comp_id,
     *     label_asym_id,label_seq_id,
     *     xesdstr,yesdstr,zesdstr,
     *     occesdstr,Besdstr,fidstr,auth_asym_id,
     *     type_symbol
        endif
      endif
C
C     Process ANISOU
C
      if (doanisou) then
         do ii = 1,6
           if (aniso(ii).gt.0.) then
             ianiso(ii) = int(aniso(ii)*1.e4+.5)
           else
             ianiso(ii) = int(aniso(ii)*1.e4-.5)
           endif
         enddo
         if (wide_code.eq."yes") then
         write(iunpdb,
     * '(a6,2x,i9,1x,3x,1x,a10,a3,3a10,6(2x,i6),10x,a4,1x,a4)')
     * 'ANISOU',id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,ianiso,auth_asym_id,
     * type_symbol
         else
         write(iunpdb,
     * '(a6,i5,1x,a4,a1,a3,1x,a1,a5,1x,6I7,2x,2a4)')
     * 'ANISOU',id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,ianiso,auth_asym_id,
     * type_symbol
       endif
       endif
C
C     PROCESS SIGUIJ
C
       if (dosiguij) then
         do ii = 1,6
           if (anisd(ii).gt.0.) then
             ianisd(ii) = int(anisd(ii)*1.e4+.5)
           else
             ianisd(ii) = int(anisd(ii)*1.e4-.5)
           endif
         enddo
         if (wide_code.eq."yes") then
         write(iunpdb,
     * '(a6,2x,i9,1x,3x,1x,a10,a3,3a10,6(2x,i6),10x,a4,1x,a4)')
     * 'SIGUIJ',id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,ianisd,auth_asym_id,
     * type_symbol
         else
         write(iunpdb,
     * '(a6,i5,1x,a4,a1,a3,1x,a1,a5,1x,6I7,2x,2a4)')
     * 'SIGUIJ',id,label_atom_id,label_alt_id,label_comp_id,
     * label_asym_id,label_seq_id,ianisd,auth_asym_id,
     * type_symbol
         endif
       endif      
       if (loop_) go to 100
       if (domodel .and. pmodnum.ne.-999999) 
     *   write(iunpdb,'(6HENDMDL)')
       return
       end
       subroutine proc_author
C
C      process AUTHOR
C
       include 'cif2pdb.cmn'
       external pdbreal
       external pdbint
       logical result_author
       character*2048 autstr,nautstr,xtemp,tail
       character*110 author,blanks
       integer kcpos
       data blanks/
     *   '                                                            '/
       pdbrec = 'AUTHOR'
       icont = 0
       ipos = 1
       iepos = 60
       kcpos = 2
       ltail = 0
       autstr = ' '
       text_ = .false.
       if (wide_code .eq. "yes") then
         iepos = 110
         kcpos = 1
       endif
       result_author= pdbstr('audit_author',
     * 'name',2048,nautstr)
       if(.not.result_author) return
 100   continue
       autstr=nautstr
 110   if(loop_.or.text_) then
         result_author= pdbstr('audit_author',
     *   'name',2048,nautstr)
       else
         nautstr = ' '
         result_author=.false.
       endif
       if (result_author .and. nautstr.eq. ' ') go to 110
       lw = nblen(autstr)
       ii = index(autstr(1:lw),',')
       if (ii.ne.0) then
         autstr(ii:ii) = ' '
         if(ii.gt.1.and.ii.lt.lw) then
           if(autstr(lw:lw).eq.'.') then
             xtemp=autstr(ii+1:lw)//autstr(1:ii-1)
             lw=lw-1
           else
             xtemp=autstr(ii+1:lw)//' '//autstr(1:ii-1)
           endif
           autstr=xtemp(1:lw)//char(0)
         endif
       endif
       if (type_code.eq.'u') autstr = tbxxupcs(autstr)
       call splitstr(nword,autstr,cstr,1024,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:max(1,nblen(cstr(ii))))//" ")
       enddo
       endif
       if(nword.gt.0.and.nautstr.ne.' ' .and. .not.text_) then
         lw = nblen(cstr(nword))
         cstr(nword)=cstr(nword)(1:lw)//','//char(0)
       endif
       do ii = 2,nword
         xtemp=cstr(1)(1:nblen(cstr(1)))//' '//
     *     cstr(ii)(1:nblen(cstr(ii)))
         cstr(1)=xtemp
       enddo
       nword = 1
       do ii = 1,nword
 200     lw = nblen(cstr(ii))
         if (lw+ipos.le.iepos+1) then
           author(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if((ipos.eq.1.and.icont.eq.0)
     *       .or.(ipos.eq.kcpos.and.icont.gt.0)) then
             author(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           author(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if (ipos.gt.iepos.or.ltail.gt.0) then               
           pdbline = pdbline+1
           if (icont.eq.0)
     *       write(iunpdb,
     *       sfmt_orig)
     *       pdbrec,author(1:iepos)
           if (icont.gt.0)
     *       write(iunpdb,
     *       sfmt_cont)
     *       pdbrec,icont+1,author(1:iepos)
           author = ' '
           ipos = kcpos
           icont = icont+1
           if (ltail.gt.0) then
             cstr(ii) = tail
             ltail = 0
             go to 200
           endif
         endif
       continue
       enddo
       if (nautstr.ne.' ') go to 100
       if ((ipos.gt.1.and.icont.eq.0)
     *   .or.(ipos.gt.kcpos.and.icont.gt.0)) then               
         pdbline = pdbline+1
         if (icont.eq.0)
     *     write(iunpdb,
     *     sfmt_orig)
     *     pdbrec,author(1:ipos-1)
         if (icont.gt.0)
     *     write(iunpdb,
     *     sfmt_cont)
     *     pdbrec,icont+1,author(1:ipos-1)
       endif
       return
       end
       subroutine proc_authrm(recstr,reccat,thisid,wcidstr)
C
C      process AUTHOR's or EDITORS for JRNL or REMARK 1
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       character*12 recstr
       character*7 thisid
       character*4 reccat
       character*15 mycat
       logical result_author,result_cid
       character*7 mycid
       character*5 wcidstr
       character*2048 autstr,nautstr,xtemp,tail
       character*100 author,blanks
       data blanks/
     *   '                                                   '/
       icont = 0
       ipos = 1
       iepos = 51
       ltail = 0
       autstr = ' '
       if (wide_code .eq. "yes") then
         iepos = 100
       endif
       text_ = .false.
       mycat = 'citation_author'
       if (reccat.eq.'EDIT') 
     *   mycat = 'citation_editor'
       result_cid = pdbstr(mycat,
     * 'citation_id',7,mycid)
       if(.not.result_cid) return
       call rjust(mycid)
       if(mycid.ne.thisid) return
       result_author= pdbstr(mycat,
     * 'name',2048,nautstr)
       if(.not.result_author) return
 100   continue
       autstr=nautstr
 110   if(loop_.or.text_) then
         result_author= pdbstr(mycat,
     *   'name',2048,nautstr)
         result_cid = pdbstr(mycat,
     *   'citation_id',7,mycid)
         if(.not.result_cid) result_author = .false.
         call rjust(mycid)
         if(mycid.ne.thisid) result_author = .false.
         if(.not.result_author) nautstr = ' '
       else
         nautstr = ' '
         result_author=.false.
       endif
       if (result_author .and. nautstr.eq. ' ') go to 110
       lw = nblen(autstr)
       ii = index(autstr(1:lw),',')
       if (ii.ne.0) then
         autstr(ii:ii) = ' '
         if(ii.gt.1.and.ii.lt.lw) then
           if(autstr(lw:lw).eq.'.') then
             xtemp=autstr(ii+1:lw)//autstr(1:ii-1)
             lw=lw-1
           else
             xtemp=autstr(ii+1:lw)//' '//autstr(1:ii-1)
           endif
           autstr=xtemp(1:lw)//char(0)
         endif
       endif
       if (type_code.eq.'u') autstr = tbxxupcs(autstr)
       call csplitstr(nword,autstr,cstr,lcstr,1024,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:lcstr(ii))//" ")
         lcstr(ii) = nblen(cstr(ii))
       enddo
       endif
       if(nword.gt.0.and.nautstr.ne.' ' .and. .not.text_) then
         lw = lcstr(nword)
         cstr(nword)=cstr(nword)(1:lw)//','//char(0)
         lcstr(nword) = lw+1
       endif
       do ii = 2,nword
         xtemp=cstr(1)(1:lcstr(1))//' '//
     *     cstr(ii)(1:lcstr(ii))
         cstr(1)=xtemp
         lcstr(1) = nblen(cstr(1))
       enddo
       nword = 1
       do ii = 1,nword
 200     lw = lcstr(ii)
         if (lw+ipos.le.iepos+1) then
           author(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if(ipos.eq.1) then
             author(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           author(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if (ipos.gt.iepos.or.ltail.gt.0) then 
           call write_remark(recstr(7:10),author(1:iepos),reccat,
     *       wcidstr,icont)              
           author = ' '
           ipos = 1
           if (ltail.gt.0) then
             cstr(ii) = tail
             ltail = 0
             go to 200
           endif
         endif
       enddo
       if (nautstr.ne.' ') go to 100
       if (ipos.gt.1) then
         call write_remark(recstr(7:10),author(1:ipos-1),reccat,
     *    wcidstr, icont)       
       endif
       return
       end
       subroutine proc_cryst1
C
C      process cryst1
C
       include 'cif2pdb.cmn'
       external typeset
       character*4 curcat
       logical result_length_a,
     *    result_length_b,
     *    result_length_c,
     *    result_angle_alpha,
     *    result_angle_beta,
     *    result_angle_gamma,
     *    result_Z_PDB,
     *    result_space_group_name
       character*13 axstr,bxstr,cxstr
       character*9 alphaxstr,betaxstr,gammaxstr
       character*5 Zstr
       character*13 SGstr
       integer kedge,kangle,ksg,kz
       real*8 a,b,c,alpha,beta,gamma
       real*8 cell(6)
       real*8 det
       integer Z
       text_ = .false.
       kedge = 9
       kangle = 7
       ksg = 11
       kz = 4
       if (wide_code .eq. 'yes' ) then
         kedge = 13
         kangle = 9
         ksg = 13
         kz = 5
       endif

C
C      initialize default cell
C
       cell_a = 1.d0
       cell_b = 1.d0
       cell_c = 1.d0
       cell_alpha = 90.d0
       cell_beta = 90.d0
       cell_gamma = 90.d0
       do ii = 1,3
       vecf2o(ii) = 0.d0
       veco2f(ii) = 0.d0
       do jj = 1,3
         matf2o(ii,jj) = 0.d0
         mato2f(ii,jj) = 0.d0
         if (ii.eq.jj) then
           matf2o(ii,jj) = 1.d0
           mato2f(jj,jj) = 1.d0
         endif
       enddo
       enddo
       a=cell_a
       b=cell_b
       c=cell_b
       alpha=cell_alpha
       beta=cell_beta
       gamma=cell_gamma
C
       result_space_group_name = pdbstr('symmetry',
     * 'space_group_name_H-M',ksg,SGstr)
       curcat = 'cell'
       pdbrec = 'CRYST1'
       pdbline = pdbline+1
       result_length_a = pdbreal(curcat,'length_a',kedge,axstr,a)       
       result_length_b = pdbreal(curcat,'length_b',kedge,bxstr,b)
       result_length_c = pdbreal(curcat,'length_c',kedge,cxstr,c)
       result_angle_alpha = pdbreal(curcat,'angle_alpha',kangle,
     *  alphaxstr,alpha)
       result_angle_beta = pdbreal(curcat,'angle_beta',kangle,
     *  betaxstr,beta)
       result_angle_gamma = pdbreal(curcat,'angle_gamma',kangle,
     *  gammaxstr,gamma)
       if (abs(alpha) .lt. 1.) alpha = 90.
       if (abs(beta) .lt. 1.) beta = 90.
       if (abs(gamma) .lt. 1.) gamma = 90.
       result_Z_PDB = pdbint(curcat,'Z_PDB',kz,Zstr,Z)
       if (result_length_a) cell_a = a
       if (result_length_b) cell_b = b
       if (result_length_c) cell_c = c
       if (result_angle_alpha) cell_alpha = alpha
       if (result_angle_beta) cell_beta = beta
       if (result_angle_gamma) cell_gamma = gamma
       if (.not.(result_length_a .and. result_length_b .and.
     *   result_length_c .and. result_angle_alpha .and.
     *   result_angle_beta .and. result_angle_gamma)) then
          call c2pwarn(' One or more default cell parameters used ')
       endif
       if (.not.result_space_group_name) then
          SGstr = 'P 1'
          call c2pwarn(' Default space group [P 1] used ')
       endif
       if (.not.result_Z_PDB) then
          Z = 1
          call c2pwarn(' Default Z [1] used ')
       endif
       if (wide_code .eq. 'yes') then
         call decjust(axstr,7)
         call decjust(bxstr,7)
         call decjust(cxstr,7)
         call decjust(alphaxstr,4)
         call decjust(betaxstr,4)
         call decjust(gammaxstr,4)
         write(iunpdb,
     *     '(6hCRYST1,12x,3x,3(1x,a13),3(1x,a9),1x,a13,1x,i5)')
     *     axstr,bxstr,cxstr,alphaxstr,betaxstr,gammaxstr,SGstr,Z
       else          
       write(iunpdb,
     * '(6hCRYST1,3f9.3,3f7.2,1x,a11,i4)')
     *  a,b,c,alpha,beta,gamma,SGstr,Z
       endif
       cell(1)=cell_a
       cell(2)=cell_b
       cell(3)=cell_c
       cell(4)=cell_alpha
       cell(5)=cell_beta
       cell(6)=cell_gamma
       call cell2mat(cell,matf2o,mato2f)
       cell_vol = det(matf2o)
       return
       end
       subroutine proc_expdta
C
C      process EXPDTA
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       logical result_method,result_details,tryd
       character*2048 expstr,nexpstr,tail
       character*110 expwds,blanks
       integer kcpos
       data blanks/
     *   '                                                            '/
       pdbrec = 'EXPDTA'
       icont = 0
       ipos = 1
       iepos = 60
       ltail = 0
       expstr = ' '
       text_ = .false.
       tryd = .false.
       kcpos = 2
       if (wide_code .eq. "yes") then
         kcpos = 1
         iepos = 110
       endif
       result_method= pdbstr('exptl',
     * 'method',2048,nexpstr)
       result_details=.false.
       if(.not.result_method) return
 100   continue
       expstr=nexpstr
 110   if(text_) then
         if(.not.tryd) then
           result_method = pdbstr('exptl',
     *     'method',2048,nexpstr)
         else
           result_details = pdbstr('exptl',
     *     'details',2048,nexpstr)
         endif
       else
         if(loop_.or.(.not.tryd)) then
           tryd = .not.tryd
           go to 110
         else
           nexpstr = ' '
           result_method=.false.
           result_details=.false.
         endif
       endif
       if ((result_method .or.result_details)
     *      .and. nexpstr.eq. ' ') go to 110
       if (type_code.eq.'u') expstr = tbxxupcs(expstr)
       call splitstr(nword,expstr,cstr,1024,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:nblen(cstr(ii)))//" ")
       enddo
       endif
       if(nword.gt.0.and.nexpstr.ne.' ' .and. .not.text_) then
         lw = nblen(cstr(nword))
         if(result_method) then
           cstr(nword)=cstr(nword)(1:lw)//','//char(0)
         else
           cstr(nword)=cstr(nword)(1:lw)//';'//char(0)
         endif
       endif
       do ii = 1,nword
 200     lw = nblen(cstr(ii))
         if (lw+ipos.le.iepos+1) then
           expwds(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if((ipos.eq.1.and.icont.eq.0)
     *       .or.(ipos.eq.kcpos.and.icont.gt.0)) then
             expwds(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           expwds(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if(ipos.le.iepos) then
           expwds(ipos:ipos) = " "
           ipos = ipos+1
         endif
         if (ipos.gt.iepos.or.ltail.gt.0) then               
           pdbline = pdbline+1
           if (icont.eq.0)
     *       write(iunpdb,
     *       sfmt_orig)
     *       pdbrec,expwds(1:iepos)
           if (icont.gt.0)
     *       write(iunpdb,
     *       sfmt_cont)
     *       pdbrec,icont+1,expwds(1:iepos)
           expwds = ' '
           ipos = kcpos
           icont = icont+1
           if (ltail.gt.0) then
             cstr(ii) = tail
             ltail = 0
             go to 200
           endif
         endif
       enddo
       if (nexpstr.ne.' ') go to 100
       if ((ipos.gt.1.and.icont.eq.0)
     *   .or.(ipos.gt.kcpos.and.icont.gt.0)) then               
         pdbline = pdbline+1
         if (icont.eq.0)
     *     write(iunpdb,
     *     sfmt_orig)
     *     pdbrec,expwds(1:ipos-1)
         if (icont.gt.0)
     *     write(iunpdb,
     *     sfmt_cont)
     *     pdbrec,icont+1,expwds(1:ipos-1)
       endif
       return
       end

       subroutine proc_ftnote
C
C      process FTNOTE
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       logical result_foot, result_footid
       common/proc_ftnote_common/result_foot, result_footid
       character*2048 footstr,tail
       character*6 footids
       character*3 contstr
       character*110 foot,blanks
       integer kcpos, kl, ikk,ik
       data blanks/
     *   '                                                            '/
       pdbrec = 'FTNOTE'
       ipos = 1
       iepos = 59
       ltail = 0
       indent = 0
       ik = 0
       
       text_ = .false.
       kcpos = 1
       if (wide_code .eq. "yes") then
         iepos = 110
       endif
 90    result_footid = pdbstr('atom_sites_footnote',
     *   'id',6,footids)
       if (.not. result_footid) return
       icont = 0
       contstr = " "
       ipos = 1
       ltail = 0
       foot = " "
       call rjust(footids)
 100   result_foot = pdbstr('atom_sites_footnote',
     *   'text',2048,footstr)
       kl = nblen(footstr)
       if (footstr(1:3).eq.'   ' .or. kl .eq. 0) indent = 1
       if (type_code.eq.'u') footstr = tbxxupcs(footstr)
       if (indent.eq.1) then
         if (kl.lt.2) then
           cstr(1) = ' '//char(0)
         else
             cstr(1) = footstr(2:kl)//char(0)
         endif
         nword = 1
       else
         call splitstr(nword,footstr,cstr,1024,' ')
       endif
       if (type_code.eq.'p') then
         do ii = 1,nword
           if (nblen(cstr(ii)).gt.0)
     *       cstr(ii) = typeset(cstr(ii)(1:nblen(cstr(ii)))//" ")
         enddo
       endif       
C
       do ii = 1,nword
         lw = max(1,nblen(cstr(ii)))
         if ((ii.eq.1 .and. cstr(1)(lw:lw).eq.":").or.indent.eq.1) then
           if ((ipos.gt.1.and.icont.eq.0)
     *     .or.(ipos.gt.kcpos.and.icont.gt.0)) then               
             pdbline = pdbline+1
             if (icont.gt.0) write(contstr,'(i3)')icont+1
             if (wide_code.eq."yes") then
               write(iunpdb,'(a)')pdbrec//'     '//
     *           footids//' '//contstr//' '//foot(1:ipos-1)
             else
               write(iunpdb,'(a)')pdbrec//footids(3:6)//
     *           ' '//foot(1:ipos-1)
             endif
             foot = ' '
             ipos = kcpos
             icont = icont+1
             numFtnote = numFtnote+1
           endif
         endif
 200   lw = nblen(cstr(ii))
       if (lw.eq.0) go to 300
       if (lw+ipos.le.iepos+1) then
         foot(ipos:lw+ipos-1) = cstr(ii)(1:lw)
       else
         if((ipos.eq.1.and.icont.eq.0)
     *     .or.(ipos.le.kcpos.and.icont.gt.0)) then
           do ikk = 1,iepos-ipos+1
             ik = iepos-ipos+2-ikk
             if (cstr(ii)(ik:ik).eq.' ') go to 210
           enddo
 210       if (ik.lt.5)ik=min(5,iepos-ipos+1)
           foot(ipos:iepos) = cstr(ii)(1:ik)
           ltail = lw - (ik)
           tail = cstr(ii)(ik+1:lw)//char(0)
           ipos = iepos
           lw = lw-ltail
         else
           tail = cstr(ii)(1:lw)//char(0)
           ltail = lw
           lw = 0
         endif
       endif
       if ((ltail.gt.0.and.ipos+lw-1.lt.iepos) .or.
     *   (indent.eq.1.and.ii.eq.nword)) then
         foot(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
       endif
       ipos = ipos+lw
       if(ipos.le.iepos) then
         foot(ipos:ipos) = " "
         ipos = ipos+1
       endif
       if ((ipos.gt.iepos.or.ltail.gt.0) .or.
     *     (ii.eq.nword.and.cstr(ii)(lw:lw).eq.";")) then
           kl = nblen(foot(1:iepos))          
           pdbline = pdbline+1
           if (icont.gt.0) write(contstr,'(i3)')icont+1
           if (kl.eq.0) then
             if (wide_code.eq."yes") then
               write(iunpdb,'(a)')pdbrec//'     '//
     *         footids//' '//contstr
             else
               write(iunpdb,'(a)')pdbrec//footids(3:5)
           endif
           else
             if (wide_code.eq."yes") then
               write(iunpdb,'(a)')pdbrec//'     '//
     *           footids//' '//contstr//' '//foot(1:kl)
             else
               write(iunpdb,'(a)')pdbrec//footids(3:6)//
     *           ' '//foot(1:kl)
             endif
           endif
           icont = icont+1
           foot = ' '
           ipos = kcpos
           numFtnote = numFtnote+1
           if (ltail.gt.0) then
             cstr(ii) = tail
             ltail = 0
             go to 200
           endif
         endif
 300   continue
       enddo
       if (text_) go to 100
       if ( (ipos.gt.1.and.icont.eq.0)
     *   .or.(ipos.gt.kcpos.and.icont.gt.0)) then               
         pdbline = pdbline+1
         if (icont.gt.0) write(contstr,'(i3)')icont+1
         if (wide_code.eq."yes") then
           write(iunpdb,'(a)')pdbrec//'     '//
     *       footids//' '//contstr//' '//foot(1:ipos-1)
         else
           write(iunpdb,'(a)')pdbrec//footids(3:6)//
     *       ' '//foot(1:ipos-1)
         endif
         numFtnote = numFtnote+1
       endif
       if (loop_) go to 90
       return
       end
       subroutine proc_header
C
C      process header 
C      (or wpdb leader using proc_leader)
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       logical result_details,
     *    result_date,
     *    result_idcode
       character*40 classstr
       character*10 cdate
       character*9 pdate,pdbdate
       character*80 idcode
       if (wide_code.eq.'yes') then
         call proc_leader
         return
       endif
       pdbrec = 'HEADER'
       pdbline = pdbline+1
       text_ = .false.
       result_details = pdbstr('struct_keywords','pdbx_keywords',
     *    40, classstr)
       if (.not.result_details) then
         result_details = pdbstr('struct_biol','details',40,classstr)
       endif
       if (type_code.eq.'u') classstr = tbxxupcs(classstr)
       if (type_code.eq.'p') classstr = tbxxupcs(classstr)
       result_date = pdbstr('database_PDB_rev','date_original',
     *  10,cdate)
       if (.not. result_date)
     * result_date = pdbstr('audit','creation_date',
     *  10,cdate)
       pdate = pdbdate(cdate)
       lme = nblen(myentry)
       ito = 0
       if (lme.ne.0) then
         kpmap = nmap
         call hash_store(myentry(1:lme)//char(0),mapstr,mapchain,
     *     NUMSTR,nmap,mhash,NUMHSH,ito)
         if (ito.eq.0)
     *     call c2perr(' More than NUMSTR strings mapped ')
         if (nmap.eq.kpmap+1) mapto(ito) = 0
       endif
       result_idcode = pdbstr('database_2','database_code',80,
     *  idcode)
       if (.not. result_idcode)
     * result_idcode = pdbstr('struct_biol','id',80,
     *  idcode)
       lidc = nblen(idcode)
       if (lidc.gt.0 .and. lme.gt.0) then
         if (idcode .ne. myentry .or. lme.gt.4) then
         if (idcode .ne. myentry )
     *     call c2pwarn(' Entry id code '//idcode//' does not match '
     *    //'command line -p argument')
         kpmap = nmap
         call hash_store
     *     ('Entry:'//idcode(1:lidc)//char(0),mapstr,mapchain,
     *     NUMSTR,nmap,mhash,NUMHSH,ifrom)
         if (ito.eq.0)
     *     call c2perr(' More than NUMSTR strings mapped ')
         if (ifrom.ne.kpmap+1.and.mapto(ifrom).ne.0) then
           call c2perr(' Previous mapping of '//idcode(1:lidc))
         endif
         mapto(ifrom) = ito
         idcode = myentry(1:4)
         lidc = lme
       endif
       endif
       if (lidc.eq.0) then
         lidc = 1
         idcode(1:1) = ' '
       endif
       write(iunpdb,
     * '(6hHEADER,4x,a40,a9,3x,a4)') classstr,pdate,idcode(1:lidc)
       return
       end
C
       subroutine proc_leader
C
C      process wpdb leader
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       logical result_details,
     *    result_date,
     *    result_idcode
       character*60 classstr
       character*10 cdate
       character*11 wpdate,wpdbdate
       character*80 idcode
       pdbrec = 'LEADER'
       pdbline = pdbline+1
       text_ = .false.
       result_details = pdbstr('struct_keywords','pdbx_keywords',
     *    60, classstr)
       if (.not.result_details) then
         result_details = pdbstr('struct_biol','details',60,classstr)
       endif
       if (type_code.eq.'u') classstr = tbxxupcs(classstr)
       if (type_code.eq.'p') classstr = tbxxupcs(classstr)
       result_date = pdbstr('database_PDB_rev','date_original',
     *  10,cdate)
       if (.not. result_date)
     * result_date = pdbstr('audit','creation_date',
     *  10,cdate)
       wpdate = wpdbdate(cdate)
       lme = nblen(myentry)
       ito = 0
       if (lme.ne.0) then
         kpmap = nmap
         call hash_store(myentry(1:lme)//char(0),mapstr,mapchain,
     *     NUMSTR,nmap,mhash,NUMHSH,ito)
         if (ito.eq.0)
     *     call c2perr(' More than NUMSTR strings mapped ')
         if (nmap.eq.kpmap+1) mapto(ito) = 0
       endif
       result_idcode = pdbstr('database_2','database_code',80,
     *  idcode)
       if (.not. result_idcode)
     * result_idcode = pdbstr('struct_biol','id',80,
     *  idcode)
       lidc = nblen(idcode)
       if (lidc.gt.0) then
         if (idcode .ne. myentry .and.lme.gt.0) then
         call c2pwarn(' Entry id code '//idcode//' does not match '
     *  //'command line -p argument')
         kpmap = nmap
         call hash_store
     *     ('Entry:'//idcode(1:lidc)//char(0),mapstr,mapchain,
     *     NUMSTR,nmap,mhash,NUMHSH,ifrom)
         if (ito.eq.0)
     *     call c2perr(' More than NUMSTR strings mapped ')
         if (ifrom.ne.kpmap+1.and.mapto(ifrom).ne.0) then
           call c2perr(' Previous mapping of '//idcode(1:lidc))
         endif
         mapto(ifrom) = ito
         idcode = myentry
         lidc = lme
       endif
       endif
       if (lidc.eq.0) then
         lidc = 1
         idcode(1:1) = ' '
       endif
       write(iunpdb,
     * '(6HLEADER,7x,4H   1,1x,3x,1x,a60,1x,a11,1x,a15)') classstr,
     *  wpdate, idcode(1:lidc)
       return
       end
C
       subroutine proc_keywds
C
C      process KEYWDS
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       logical result_keywds
       character*2048 keystr,nkeystr,tail
       character*110 keywds,blanks
       integer kcpos
       data blanks/
     *   '                                                            '/
       pdbrec = 'KEYWDS'
       icont = 0
       ipos = 1
       iepos = 60
       kcpos = 2
       ltail = 0
       keystr = ' '
       text_ = .false.
       if (wide_code .eq. "yes") then
         iepos = 110
         kcpos = 1
       endif
       result_keywds= pdbstr('struct_keywords',
     * 'text',2048,nkeystr)
       if(.not.result_keywds) return
 100   continue
       keystr=nkeystr
 110   if(loop_.or.text_) then
         result_keywds= pdbstr('struct_keywords',
     *   'text',2048,nkeystr)
       else
         nkeystr = ' '
         result_keywds=.false.
       endif
       if (result_keywds .and. nkeystr.eq. ' ') go to 110
       if (type_code.eq.'u') keystr = tbxxupcs(keystr)
       call splitstr(nword,keystr,cstr,1024,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:nblen(cstr(ii)))//" ")
       enddo
       endif
       if(nword.gt.0.and.nkeystr.ne.' ' .and. .not.text_) then
         lw = nblen(cstr(nword))
         cstr(nword)=cstr(nword)(1:lw)//','//char(0)
       endif
       do ii = 1,nword
 200     lw = nblen(cstr(ii))
         if (lw+ipos.le.iepos+1) then
           keywds(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if((ipos.eq.1.and.icont.eq.0)
     *       .or.(ipos.le.kcpos.and.icont.gt.0)) then
             keywds(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           keywds(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if(ipos.le.iepos) then
           keywds(ipos:ipos) = " "
           ipos = ipos+1
         endif
         if (ipos.gt.iepos.or.ltail.gt.0) then               
           pdbline = pdbline+1
           if (icont.eq.0)
     *       write(iunpdb,
     *       sfmt_orig)
     *       pdbrec,keywds(1:iepos)
           if (icont.gt.0)
     *       write(iunpdb,
     *       sfmt_cont)
     *       pdbrec,icont+1,keywds(1:iepos)
           keywds = ' '
           ipos = kcpos
           icont = icont+1
           if (ltail.gt.0) then
             cstr(ii) = tail
             ltail = 0
             go to 200
           endif
         endif
       continue
       enddo
       if (nkeystr.ne.' ') go to 100
       if ((ipos.gt.1.and.icont.eq.0)
     *   .or.(ipos.gt.kcpos.and.icont.gt.0)) then               
         pdbline = pdbline+1
         if (icont.eq.0)
     *     write(iunpdb,
     *     sfmt_orig)
     *     pdbrec,keywds(1:ipos-1)
         if (icont.gt.0)
     *     write(iunpdb,
     *     sfmt_cont)
     *     pdbrec,icont+1,keywds(1:ipos-1)
       endif
       return
       end

       subroutine proc_master
C
C      process master and end
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external pdbstr
       external typeset
       pdbrec = 'MASTER'
       pdbline = pdbline+1
       if (wide_code.eq."yes") then
         write(iunpdb,
     *   '(6hMASTER,16x,12I9)')
     *   numRemark,numFtnote,numHet,numHelix,numSheet,numTurn,
     *   numSite,numXform,numCoord,numTer,numConect,numSeq
       else
         write(iunpdb,
     *   '(6hMASTER,4x,12I5)')
     *   numRemark,numFtnote,numHet,numHelix,numSheet,numTurn,
     *   numSite,numXform,numCoord,numTer,numConect,numSeq
       endif
       pdbrec = 'END'
       pdbline = pdbline+1
       write(iunpdb,'(3hEND)')
       call c2pmsg(' ',' Processing Finished')
       return
       end
       subroutine proc_origx
C
C      process origx
C
       include 'cif2pdb.cmn'
       external typeset
       external pdbint
       external pdbstr
       character*19 curcat
       logical result_origx11,
     *    result_origx12,
     *    result_origx13,
     *    result_origx21,
     *    result_origx22,
     *    result_origx23,
     *    result_origx31,
     *    result_origx32,
     *    result_origx33,
     *    result_origxv1,
     *    result_origxv2,
     *    result_origxv3
       common/proc_origx_common/result_origx11,
     *    result_origx12,
     *    result_origx13,
     *    result_origx21,
     *    result_origx22,
     *    result_origx23,
     *    result_origx31,
     *    result_origx32,
     *    result_origx33,
     *    result_origxv1,
     *    result_origxv2,
     *    result_origxv3
       character*13 x11str,x12str,x13str,x21str,x22str,x23str,
     *    x31str,x32str,x33str,v1str,v2str,v3str
       real*8 x11,x12,x13,x21,x22,x23,x31,x32,x33,v1,v2,v3
       integer knum
       curcat = 'database_pdb_matrix'
       pdbrec = 'ORIGX'
       pdbline = pdbline+1
       knum = 10
       if (wide_code.eq.'yes') knum = 13
       result_origx11 = pdbreal(curcat,'origx[1][1]',knum,x11str,x11)
       result_origx12 = pdbreal(curcat,'origx[1][2]',knum,x12str,x12)
       result_origx13 = pdbreal(curcat,'origx[1][3]',knum,x13str,x13)
       result_origxv1 = pdbreal(curcat,'origx_vector[1]',knum,
     *   v1str,v1)
       pdbline = pdbline+1
       result_origx21 = pdbreal(curcat,'origx[2][1]',knum,x21str,x21)
       result_origx22 = pdbreal(curcat,'origx[2][2]',knum,x22str,x22)
       result_origx23 = pdbreal(curcat,'origx[2][3]',knum,x23str,x23)
       result_origxv2 = pdbreal(curcat,'origx_vector[2]',knum,
     *   v2str,v2)
       pdbline = pdbline+1
       result_origx31 = pdbreal(curcat,'origx[3][1]',knum,x31str,x31)
       result_origx32 = pdbreal(curcat,'origx[3][2]',knum,x32str,x32)
       result_origx33 = pdbreal(curcat,'origx[3][3]',knum,x33str,x33)
       result_origxv3 = pdbreal(curcat,'origx_vector[3]',knum,
     *   v3str,v3)
       if(result_origx11) then
         if (wide_code .eq. 'yes' ) then
           call decjust(x11str,7)
           call decjust(x12str,7)
           call decjust(x13str,7)
           call decjust(v1str,8)
           call decjust(x21str,7)
           call decjust(x22str,7)
           call decjust(x23str,7)
           call decjust(v2str,8)
           call decjust(x31str,7)
           call decjust(x32str,7)
           call decjust(x33str,7)
           call decjust(v3str,8)
           write(iunpdb,
     *       '(5hORIGX,i1,12x,3x,3(1x,a13),1x,a13)') 
     *       1,x11str,x12str,x13str,v1str
           write(iunpdb,
     *       '(5hORIGX,i1,12x,3x,3(1x,a13),1x,a13)') 
     *       2,x21str,x22str,x23str,v2str
           write(iunpdb,
     *       '(5hORIGX,i1,12x,3x,3(1x,a13),1x,a13)') 
     *       3,x31str,x32str,x33str,v3str
         else
           write(iunpdb,
     *       '(5hORIGX,i1,4x,3f10.6,5x,f10.5)') 
     *       1,x11,x12,x13,v1,
     *       2,x21,x22,x23,v2,
     *       3,x31,x32,x33,v3
         endif
         numXform = numXform+3
       else
         pdbline = pdbline-3
       endif
       return
       end
       subroutine proc_scale
C
C      process scale
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbstr
       external typeset
       character*10 curcat
       logical result_scale11,
     *    result_scale12,
     *    result_scale13,
     *    result_scale21,
     *    result_scale22,
     *    result_scale23,
     *    result_scale31,
     *    result_scale32,
     *    result_scale33,
     *    result_scalev1,
     *    result_scalev2,
     *    result_scalev3
       character*13 x11str,x12str,x13str,x21str,x22str,x23str,
     *    x31str,x32str,x33str,v1str,v2str,v3str
       real*8 x11,x12,x13,x21,x22,x23,x31,x32,x33,v1,v2,v3
       integer knum
       character*12 fts
       fts = 'fract_transf'
       curcat = 'atom_sites'
       pdbrec = 'SCALE'
       pdbline = pdbline+1
       knum = 10
       if (wide_code.eq.'yes') knum = 13
       result_scale11 = pdbreal(curcat,fts//'_matrix[1][1]',
     *   knum,x11str,x11)
       result_scale12 = pdbreal(curcat,fts//'_matrix[1][2]',
     *   knum,x12str,x12)
       result_scale13 = pdbreal(curcat,fts//'_matrix[1][3]',
     *   knum,x13str,x13)
       result_scalev1 = pdbreal(curcat,fts//'_vector[1]',knum,
     *   v1str,v1)
       pdbline = pdbline+1
       result_scale21 = pdbreal(curcat,fts//'_matrix[2][1]',
     *   knum,x21str,x21)
       result_scale22 = pdbreal(curcat,fts//'_matrix[2][2]',
     *   knum,x22str,x22)
       result_scale23 = pdbreal(curcat,fts//'_matrix[2][3]',
     *   knum,x23str,x23)
       result_scalev2 = pdbreal(curcat,fts//'_vector[2]',knum,
     *   v2str,v2)
       pdbline = pdbline+1
       result_scale31 = pdbreal(curcat,fts//'_matrix[3][1]',
     *   knum,x31str,x31)
       result_scale32 = pdbreal(curcat,fts//'_matrix[3][2]',
     *   knum,x32str,x32)
       result_scale33 = pdbreal(curcat,fts//'_matrix[3][3]',
     *   knum,x33str,x33)
       result_scalev3 = pdbreal(curcat,fts//'_vector[3]',knum,
     *   v3str,v3)
       if (.not.(result_scale11.and.result_scale12.and.result_scale13
     *    .and.result_scale21.and.result_scale22.and.result_scale23
     *    .and.result_scale31.and.result_scale32.and.result_scale33
     *    .and.result_scalev1.and.result_scalev2.and.result_scalev3))
     *    then
       x11 = mato2f(1,1)
       x12 = mato2f(1,2)
       x13 = mato2f(1,3)
       x21 = mato2f(2,1)
       x22 = mato2f(2,2)
       x23 = mato2f(2,3)
       x31 = mato2f(3,1)
       x32 = mato2f(3,2)
       x33 = mato2f(3,3)
       v1 = veco2f(1)
       v2 = veco2f(2)
       v3 = veco2f(3)
       call c2pwarn(' Default SCALE used ')
       else
       mato2f(1,1) = x11
       mato2f(1,2) = x12
       mato2f(1,3) = x13
       mato2f(2,1) = x12
       mato2f(2,2) = x22
       mato2f(2,3) = x23
       mato2f(3,1) = x31
       mato2f(3,2) = x32
       mato2f(3,3) = x33
       veco2f(1) = v1
       veco2f(2) = v2
       veco2f(3) = v3
       endif
       call invxfrm(mato2f,veco2f,matf2o,vecf2o)
       if (wide_code .eq.'yes' ) then
         write(iunpdb,
     *   '(5hSCALE,i1,12x,3x,3(1x,f13.6),1x,f13.5)') 
     *    1,x11,x12,x13,v1
         write(iunpdb,
     *   '(5hSCALE,i1,12x,3x,3(1x,f13.6),1x,f13.5)') 
     *    2,x21,x22,x23,v2
         write(iunpdb,
     *   '(5hSCALE,i1,12x,3x,3(1x,f13.6),1x,f13.5)') 
     *    3,x31,x32,x33,v3
       else
         write(iunpdb,
     *   '(5hSCALE,i1,4x,3f10.6,5x,f10.5)') 1,x11,x12,x13,v1,
     *    2,x21,x22,x23,v2,3,x31,x32,x33,v3
       endif
       numXform = numXform+3
       return
       end
       subroutine write_remark(remids,remark,subtype,wcidstr,icont)
C
C      write the next REMARK line with id remids
C      and text remark with continuation flag icont
C
C      The general wide format structure is:
C         Columns  1 --  6:  pdbrec
C         Columns  8 -- 11:  remids
C         Columns 12 -- 17:  wcidstr
C         Columns 19 -- 21:  continuation flag
C         Columns 23 -- 26:  subtype, if non-blank
C         Columns 28 --132:  remark
C
       include 'cif2pdb.cmn'
       external pdbreal
       external pdbint
       external pdbstr
       external typeset
       character *(*) remids,remark,subtype,wcidstr
       character *132 remstr
       integer icont
       integer lc
       character*3 contstr
       
       contstr = ' '
       if (wide_code.eq."yes" .or. remids .eq. "   1") then
         remstr = subtype//' '//remark
       else
         remstr = remark
       endif
       if (icont .eq. 0 .and. remids.ne."   1"
     *    .and. pdbrec.eq."REMARK") then
         if (wide_code.eq."yes") then
           write(iunpdb,
     *     '(a6,1x,a4)')
     *     pdbrec,remids
         else
           write(iunpdb,
     *       '(a6,a4)')
     *        pdbrec,remids
         endif
         icont = icont+1
         if (pdbrec.eq."REMARK") numRemark = numRemark+1
         pdbline = pdbline+1
       endif
       if (icont.gt.0) write(contstr,'(i3)') icont+1
       lc = nblen(remstr)
       pdbline = pdbline+1
       if (wide_code.eq."yes") then
         if (lc.gt.0) then
         write(iunpdb,
     *   '(a6,1x,a4,1x,a5,1x,a3,1x,a)')
     *   pdbrec,remids,wcidstr,contstr,remstr(1:lc)
         else
         write(iunpdb,
     *   '(a6,1x,a4,1x,a5,1x,a3)')
     *   pdbrec,remids,wcidstr,contstr
         endif
       else
         if (subtype.eq.'    ') then
           if (lc.gt.0) then
           write(iunpdb,
     *       '(a6,a4,1x,a)')
     *       pdbrec,remids,remstr(1:lc)
           else
           write(iunpdb,
     *       '(a6,a4)')
     *       pdbrec,remids
           endif
         else
           lc = nblen(remark)
           if (lc.gt.0) then
           write(iunpdb,
     *       '(a6,a4,2x,a4,a2,1x,a)')
     *       pdbrec,remids,subtype,contstr(2:3),remark(1:lc)
           else
           write(iunpdb,
     *       '(a6,a4,2x,a4,a2)')
     *       pdbrec,remids,subtype,contstr(2:3)
           endif
         endif
       endif
       icont = icont+1
       if (pdbrec.eq."REMARK") numRemark = numRemark+1
       return
       end
       subroutine proc_remark
C
C      process REMARK and JRNL
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       logical bkmrk_
       logical result_remark,result_remid
       logical result_cid, result_bmk, have_primary
       dimension iremno(999)
       parameter (MARKS=200)
       character*7 citjl(MARKS),cital(MARKS),citel(MARKS)
       character*5 resstr
       integer markjl(MARKS),markal(MARKS),markel(MARKS)
       integer markjl2(MARKS)
       integer nordjl(MARKS),nordal(MARKS),nordel(MARKS)
       integer rembks(999)
       integer kremno,ii,iii,kk,jj,kremnot,ifirst
       integer kcpos, icont, lc
       character*4 remids
       character*7 cidstr
       character*5 wcidstr
       character*2048 remstr,tail
       character*111 remark,blanks
       common/proc_remark_common/result_remark,result_remid,result_cid,
     *   result_bmk
       data iremno /999*0/
       data blanks/
     *   '                                                           '/
       pdbrec = 'JRNL'
       ncid = 0
       ncidal = 0
       ncidel = 0
       ifirst = 0
       have_primary = .false.
       iepos = 59
       kcpos = 1
       if (wide_code .eq. "yes") then
         iepos = 111
         kcpos = 2
       endif
  10   result_cid = pdbstr('citation','id',7,cidstr)
       if (.not.result_cid) go to 20
       if (cidstr.eq.'primary') have_primary = .true.
       call rjust(cidstr)
       if (ncid.eq.0) then
         citjl(1) = cidstr
         markjl(1) = 0
         result_bmk = bkmrk_(markjl(1))
         markjl2(1) = 0
         result_bmk = bkmrk_(markjl2(1))
         ncid = 1
       else
         if (citjl(ncid).ne.cidstr) then
           ncid = ncid+1
           if (ncid.gt.MARKS) 
     *       call c2perr(' Overflow of citjl array, increase MARKS')
           citjl(ncid) = cidstr
           markjl(ncid) = 0
           markjl2(ncid) = 0
           result_bmk = bkmrk_(markjl(ncid))
           result_bmk = bkmrk_(markjl2(ncid))
         endif
       endif
       if (loop_) go to 10
       call ssort(citjl,ncid,7,nordjl)
  20   result_cid = pdbstr('citation_author','citation_id',7,cidstr)
       if(.not.result_cid) go to 30
       call rjust(cidstr)
       if (cidstr.eq.'primary') have_primary = .true.
       if (ncidal.eq.0) then
         cital(1) = cidstr
         markal(1) = 0
         result_bmk = bkmrk_(markal(1))
         ncidal = 1
       else
         if (cital(ncidal).ne.cidstr) then
           ncidal = ncidal+1
           if (ncidal.gt.MARKS) 
     *       call c2perr(' Overflow of cital array, increase MARKS')
           cital(ncidal) = cidstr
           markal(ncidal) = 0
           result_bmk = bkmrk_(markal(ncidal))
         endif
       endif
       if (loop_) go to 20
       call ssort(cital,ncidal,7,nordal)
  30   result_cid = pdbstr('citation_editor','citation_id',7,cidstr)
       if(.not.result_cid) go to 40
       if (cidstr.eq.'primary') have_primary = .true.
       call rjust(cidstr)
       if (ncidel.eq.0) then
         citel(1) = cidstr
         markel(1) = 0
         result_bmk = bkmrk_(markel(1))
         ncidel = 1
       else
         if (citel(ncidel).ne.cidstr) then
           ncidel = ncidel+1
           if (ncidel.gt.MARKS) 
     *       call c2perr(' Overflow of citel array, increase MARKS')
           citel(ncidel) = cidstr
           markel(ncidel) = 0
           result_bmk = bkmrk_(markel(ncidel))
         endif
       endif
       if (loop_) go to 30
       call ssort(citel,ncidel,7,nordel)
  40   continue
       if (have_primary) then
         do ii = 1,ncidal
         if (cital(ii).eq.'primary') then
           result_bmk = bkmrk_(markal(ii))
           call proc_authrm('JRNL        ','AUTH',cital(ii),'     ')
         endif
         enddo
         do ii = 1,ncid
         if (citjl(ii).eq.'primary') then
           result_bmk = bkmrk_(markjl(ii))
           call proc_titlrm('JRNL        ','TITL',citjl(ii),'     ')
         endif
         enddo
         do ii = 1,ncidel
         if (citel(ii).eq.'primary') then
           result_bmk = bkmrk_(markel(ii))
           call proc_authrm('JRNL        ','EDIT',citel(ii),'     ')
         endif
         enddo
         do ii = 1,ncid
         if (citjl(ii).eq.'primary') then
           result_bmk = bkmrk_(markjl2(ii))
           call proc_refrm('JRNL        ','REF ',citjl(ii),'     ')
         endif
         enddo
       endif
       pdbrec = 'REMARK'
       do iii = 1,ncid
         kk = nordjl(iii)
         if (citjl(kk).ne.'primary') then
         cidstr = citjl(kk)
         call ljust(cidstr)
         wcidstr = cidstr(1:5)
         call rjust(wcidstr)
         if (wide_code .eq. "yes") then
           write(iunpdb,'(a)') 'REMARK    1 '//wcidstr
           write(iunpdb,'(a)') 'REMARK    1 '//wcidstr
     *      //'          REFERENCE '//wcidstr
           pdbline = pdbline+2
           numRemark = numRemark+2
         else
           if (ifirst.eq.0) then
             write(iunpdb,'(a)') 'REMARK   1'
             pdbline = pdbline+1
             numRemark = numRemark+1
           endif
           write(iunpdb,'(a)') 'REMARK   1 REFERENCE '//cidstr
           pdbline = pdbline+1
           numRemark = numRemark+1
         endif
         ifirst = 1
         do ii = 1,ncidal
         jj = nordal(ii)
         if (cital(jj).eq.citjl(kk)) then
           result_bmk = bkmrk_(markal(jj))
           call proc_authrm('REMARK   1  ','AUTH',cital(jj),wcidstr)
         endif
         enddo
         result_bmk = bkmrk_(markjl(kk))
         call proc_titlrm('REMARK   1  ','TITL',citjl(kk),wcidstr)
         do ii = 1,ncidel
         jj = nordel(ii)
         if (citel(jj).eq.citjl(kk)) then
           result_bmk = bkmrk_(markel(jj))
           call proc_authrm('REMARK   1  ','EDIT',citel(jj),wcidstr)
         endif
         enddo
         result_bmk = bkmrk_(markjl2(kk))
         call proc_refrm('REMARK   1  ','REF ',citjl(kk),wcidstr)
         endif
       enddo
C
C      Prescan for remarks in the range of 1-999
C
       do kremno = 1,999
         iremno(kremno) = 0
         rembks(kremno) = 0
       enddo
  80   result_remid = pdbstr('database_PDB_remark','id',4,remids)
       if(.not.result_remid) go to 83
       call rjust(remids)
       read(remids,'(i4)',err=81) kremno
       go to 82
  81   kremno = 999
  82   if(kremno.ge.1.and.kremno.le.999) then
         iremno(kremno)=iremno(kremno)+1
         if (iremno(kremno).eq.1) result_remid = bkmrk_(rembks(kremno))
       endif
       if (loop_) go to 80
C
C      Prepare the resultion for REMARK 2
C
  83   result_cid = pdbstr('reflns','d_resolution_high',5,resstr)

C
       ipos = 1
       ltail = 0
       indent = 0
       text_ = .false.
       remark = ' '
       kremnot = 0
  85   kremnot = kremnot+1
       if (kremnot.eq.2.and.iremno(2).eq.0) then
         icont = 0
         if (result_cid) then
           call rjust(resstr)
           if (index(resstr,'.').eq.0) then
             call write_remark("   2",
     *         "RESOLUTION."//resstr(2:5)//". ANGSTROMS.",
     *         '    ','     ',icont)
           else
             call write_remark("   2",
     *   "RESOLUTION."//resstr//" ANGSTROMS.",
     *   '    ','     ',icont)
           endif
         else
           call write_remark("   2",
     *     "RESOLUTION. NOT APPLICABLE.",'    ',
     *     '     ',icont)
         endif
       endif
       if (iremno(kremnot).gt.0) then
       result_remid = bkmrk_(rembks(kremnot))
  90   result_remid = pdbstr('database_PDB_remark','id',4,remids)
       if(.not.result_remid) goto 85
       call rjust(remids)
       read(remids,'(i4)',err=92) kremno
       go to 95
  92   kremno = 999
  95   if(kremno .ne. kremnot) go to 320
       icont = 0
 100   continue
       result_remark= pdbstr('database_PDB_remark',
     * 'text',2048,remstr)
       indent = 0
       if (remstr(1:2).eq.'  ') indent = 1
       lwr = nblen(remstr)
       if (lwr .eq. 0 .and. icont.eq.0) go to 310
       if (type_code.eq.'u') remstr = tbxxupcs(remstr)
       lc = index(remstr,": ")
       if (remstr(1:1).eq.' ' .and.
     *   lwr.gt.1 .and.
     *   lwr-1+ipos.le.iepos+1 .and.
     *   ipos .eq. 1 .and.
     *   (lc.eq.0.or.lc.gt.lwr-3)) then
         nword = 1
         cstr(1) = remstr(2:lwr)//char(0)
         lcstr(1) = lwr-1
       else
         if (remstr(1:6).eq.'REMARK') then
           indent = 1
           nword = 1
           cstr(1) = remstr(12:70)//char(0)
           lcstr(1) = 59
         else
           if ((lc .ne. 0 .and. lc .le. lwr-3) .or. 
     *       kremno .lt. 6 .or. kremno .gt. 99) then
             indent = 1
             nword = 1
             if (remstr(1:1).ne.' ') then
               cstr(1) = remstr(1:lwr)//char(0)
               lcstr(1) = lwr
             else
               cstr(1) = remstr(2:lwr)//char(0)
               lcstr(1) = lwr-1
             endif
           else
             call csplitstr(nword,remstr,cstr,lcstr,1024,' ')
           endif
         endif
       endif
       if (indent.ne.0 .and. ipos .gt. 1) then
         lc = max(1,nblen(remark(1:ipos-1)))
         call write_remark(remids,remark(1:lc),'    ','     ',icont)
         remark = ' '
         ipos = 1
       endif
       if ((indent .ne. 0 
     *       .and. remstr(1:1) .eq. ' ') .or. lwr .eq. 0 ) then
         if (lwr.gt.1) remark(1:lwr) = remstr(2:lwr)//char(0)
         lc = max(1,nblen(remark))
         call write_remark(remids,remark(1:lc),'    ','     ',icont)
         remark = ' '
         ipos = 1
         go to 310
       endif
       do ii = 1,nword
 200     lw = lcstr(ii)
         if (lw+ipos.le.iepos+1) then
           remark(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if((ipos.eq.1.and.icont.eq.0)
     *       .or.(ipos.le.kcpos.and.icont.gt.0)) then
             remark(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             tail = cstr(ii)(iepos-ipos:lw)
             ltail = lw - (iepos-ipos+1)
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         ipos = ipos+lw
         if(ipos.le.iepos) then
           remark(ipos:ipos) = " "
           ipos = ipos+1
         endif
         if ((ipos.gt.iepos.or.ltail.gt.0) .or.
     *     (indent.eq.1.and.ii.eq.nword)) then
           call write_remark(remids,
     *       remark(1:min(iepos,max(1,ipos-1))),
     *       '    ','     ',icont)
           remark = ' '
           ipos = 1
           if (ltail.gt.0) then
             cstr(ii) = tail
             lcstr(ii) = nblen(cstr(ii))
             ltail = 0
             go to 200
           endif
         endif
       enddo
 310   continue
       if (text_) go to 100
       if (ipos.gt.1) then
         call write_remark(remids,remark(1:ipos-1),
     *    '    ','     ',icont)
         ipos = 1
         remark = ' '
       endif
 320   if (loop_.and.iremno(kremnot).gt.1) go to 90
       endif
       if (kremnot.lt.999) go to 85
       if (wide_code.eq.'yes') then
         ii = 9600
         go to 410
       endif
       do ii = 960,950,-1
       if (iremno(ii).eq.0) go to 410
       enddo
       ii=96
       call c2pwarn(' REMARK 950-960 already used ')
 410   call tagchk(ii,iunpdb,wide_code,pdbline,numRemark)
       return
       end
C
       subroutine write_csdet(pdbtag,cscat,cstag)
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       character *(*) pdbtag,cscat,cstag
       character*2049 csstr
       character*2149 xcsstr
       character*80 alttag
       integer icont, iepos, kltcs, ibpos, i, ii
       logical result_csdet
       
       ibpos = 11
       iepos = 60
       if (wide_code .eq. "yes") then
         ibpos = 23
         iepos = 110
       endif
       alttag = cstag
       result_csdet = pdbstr(cscat,cstag,2048,csstr)
       if (.not.result_csdet) then
         alttag = 'pdbx_'//cstag
         result_csdet = 
     *         pdbstr(cscat,alttag(1:nblen(alttag)),2048,csstr)
       endif
       if (.not.result_csdet) return
       if (type_code .eq. 'u') 
     *       csstr(1:2048) = tbxxupcs(csstr(1:2048))
       if (type_code .eq. 'p') then
         csstr(1:2048) = typeset(csstr(1:max(1,nblen(csstr)))//' ')
           endif
       kl = nblen(csstr(1:2048))
           csstr(kl+1:kl+1)=';'       
       if (pdbrec.eq.'SOURCE') icont = numSource
       if (pdbrec.eq.'COMPND') icont = numCompound
       xcsstr =pdbtag//': '//csstr(1:kl+1) 
       call write_title(xcsstr(1:nblen(xcsstr)),icont,iepos)
 100   if (.not. text_) return
       result_csdet = pdbstr(cscat,alttag(1:nblen(alttag)),2048,csstr)
       if (.not.result_csdet) return
       if (type_code .eq. 'u') 
     *      csstr(1:2048) = tbxxupcs(csstr(1:2048))
       if (type_code .eq. 'p') then
         csstr(1:2048) = typeset(csstr(1:max(1,nblen(csstr)))//' ')
       endif
       kltcs = 0
       if (tcsbuf(1:1).ne.char(0)) then
         kltcs = max(1,nblen(tcsbuf))
         if (tcsbuf(kltcs:kltcs).eq.';') then
           tcsbuf(kltcs:kltcs) = ' '
           kltcs = max(1,nblen(tcsbuf))
         endif
       endif
       kl = nblen(csstr(1:2048))
       csstr(kl+1:kl+1) = ';'
       if (pdbrec.eq.'SOURCE') icont = numSource
       if (pdbrec.eq.'COMPND') icont = numCompound
       if (csstr(1:2).eq.'  ' 
     *   .or. kltcs.ge.iepos+ibpos-5 .or. kltcs.eq.0) then
         call write_title(csstr(1:kl+1),icont,iepos)
         go to 100
       endif
       do ii = 1,min(kl+1,iepos+ibpos-kltcs-1)
         i = 1+min(kl+1,iepos+ibpos-kltcs-1)-ii
         if (csstr(i:i).eq.' ') then
           tcsbuf(kltcs+2:ibpos+iepos-1)=csstr(1:i)
           if (i.ge.kl) go to 100
           call write_title(csstr(i+1:kl+1),icont,iepos)
           go to 100
         endif
       enddo
       call write_title(csstr(1:kl+1),icont,iepos)
       go to 100
       end

       subroutine proc_compnd
C
C      Process COMPND
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       external pdbstr
C
C      MOL_ID:             _entity.id
C      MOLECULE:           _entity.pdbx_description or
C                             value of _pdbx_entity_name.name for which
C                            _pdbx_entity_name.entity_id matches
C                            _entity.id and for which
C                            _pdbx_entity_name.name_type is
C                            'RCSB_NAME'
C      CHAIN:              values of _struct_asym.id for which
C                            _struct_asym.entity_id matches
C                            _entity.id
C      FRAGMENT:           _entity.pdbx_fragment
C      SYNONYM:            values of _pdbx_entity_name.name for which
C                            _pdbx_entity_name.entity_id matches
C                            _entity.id and for which
C                            _pdbx_entity_name.name_type is
C                            'RCSB_SYNONYM'
C                          or values of _entity_name_com.name for which
C                            _entity_name_com.entity_id matches
C                            _entity.id
C
C      ENGINEERED: YES;    if _entity.src_method is not "nat" 
C      BIOLOGICAL_UNIT:     
C      MUTATION:           _entity.pdbx_mutation            
C      EC:                 _entity.pdbx_ec
       return
       end
       
       subroutine proc_source
C
C      Process SOURCE
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       character*2049 eidstr
       logical result_source
       integer iepos, kl, icont
       
       iepos = 60
       if (wide_code.eq."yes") iepos = 110
       
       pdbrec = 'SOURCE'

 100   result_source= pdbstr('entity_src_nat',
     *  'entity_id',2048,eidstr)
       kl = nblen(eidstr(1:2048))
       eidstr(kl+1:kl+1)=';'
       if (.not.result_source) go to 300
       icont = numSource
       call write_title('MOL_ID: '//eidstr(1:kl+1),icont,iepos)
       call write_csdet('FRAGMENT','entity_src_nat',
     *   'fragment')
       call write_csdet('ORGANISM_SCIENTIFIC','entity_src_nat',
     *   'organism_scientific')
       call write_csdet('ORGANISM_COMMON','entity_src_nat',
     *   'organism_common')
       call write_csdet('STRAIN','entity_src_nat',
     *   'strain')
       call write_csdet('VARIANT','entity_src_nat',
     *   'variant')
       call write_csdet('CELL_LINE','entity_src_nat',
     *   'cell_line')
       call write_csdet('ATCC','entity_src_nat',
     *   'ATCC')
       call write_csdet('ORGAN','entity_src_nat',
     *   'organ')
       call write_csdet('TISSUE','entity_src_nat',
     *   'tissue')
       call write_csdet('TISSUE_FRACTION','entity_src_nat',
     *   'tissue_fraction')
       call write_csdet('CELL','entity_src_nat',
     *   'cell')
       call write_csdet('ORGANELLE','entity_src_nat',
     *   'organelle')
       call write_csdet('SECRETION','entity_src_nat',
     *   'secretion')
       call write_csdet('CELLULAR_LOCATION','entity_src_nat',
     *   'cellular_location')
       call write_csdet('PLASMID','entity_src_nat',
     *   'plasmid_name')
       call write_csdet('GENE','entity_src_nat',
     *   'gene')
       if (loop_) go to 100
       

300       result_source= pdbstr('entity_src_gen',
     * 'entity_id',2048,eidstr)
       if (.not.result_source) go to 600
       kl = nblen(eidstr(1:2048))
       eidstr(kl+1:kl+1)=';'
       icont = numSource
       call write_title('MOL_ID: '//eidstr(1:kl+1),icont,iepos)
       call write_csdet('FRAGMENT','entity_src_gen',
     *   'gene_src_fragment')
       call write_csdet('ORGANISM_SCIENTIFIC','entity_src_gen',
     *   'gene_src_scientific_name')
       call write_csdet('ORGANISM_COMMON','entity_src_gen',
     *   'gene_src_common_name')
       call write_csdet('STRAIN','entity_src_gen',
     *   'gene_src_strain')
       call write_csdet('VARIANT','entity_src_gen',
     *   'gene_src_variant')
       call write_csdet('CELL_LINE','entity_src_gen',
     *   'gene_src_cell_line')
       call write_csdet('ATCC','entity_src_gen',
     *   'gene_src_ATCC')
       call write_csdet('ORGAN','entity_src_gen',
     *   'gene_src_organ')
       call write_csdet('TISSUE','entity_src_gen',
     *   'gene_src_tissue')
       call write_csdet('TISSUE_FRACTION','entity_src_gen',
     *   'gene_src_tissue_fraction')
       call write_csdet('CELL','entity_src_gen',
     *   'gene_src_cell')
       call write_csdet('ORGANELLE','entity_src_gen',
     *   'gene_src_organelle')
       call write_csdet('SECRETION','entity_src_gen',
     *   'gene_src_secretion')
       call write_csdet('CELLULAR_LOCATION','entity_src_gen',
     *   'gene_src_cellular_location')
       call write_csdet('GENE','entity_src_gen',
     *   'gene_src_gene')
       call write_csdet('EXPRESSION_SYSTEM','entity_src_gen',
     *   'host_org_scientific_name')
       call write_csdet('EXPRESSION_SYSTEM_COMMON','entity_src_gen',
     *   'host_org_common_name')
       call write_csdet('EXPRESSION_SYSTEM_FRAGMENT',
     *'entity_src_gen',
     *   'host_org_fragment')
       call write_csdet('EXPRESSION_SYSTEM_STRAIN','entity_src_gen',
     *   'host_org_strain')
       call write_csdet('EXPRESSION_SYSTEM_VARIANT','entity_src_gen',
     *   'host_org_variant')
       call write_csdet('EXPRESSION_SYSTEM_CELL_LINE','entity_src_gen',
     *   'host_org_cell_line')
       call write_csdet('EXPRESSION_SYSTEM_ATCC','entity_src_gen',
     *   'host_org_ATCC')
       call write_csdet('EXPRESSION_SYSTEM_ORGAN','entity_src_gen',
     *   'host_org_organ')
       call write_csdet('EXPRESSION_SYSTEM_TISSUE','entity_src_gen',
     *   'host_org_tissue')
       call write_csdet('EXPRESSION_SYSTEM_TISSUE_FRACTION',
     *   'entity_src_gen',
     *   'host_org_tissue_fraction')
       call write_csdet('EXPRESSION_SYSTEM_CELL','entity_src_gen',
     *   'host_org_cell')
       call write_csdet('EXPRESSION_SYSTEM_ORGANELLE','entity_src_gen',
     *   'host_org_organelle')
       call write_csdet('EXPRESSION_SYSTEM_SECRETION','entity_src_gen',
     *   'host_org_secretion')
       call write_csdet('EXPRESSION_SYSTEM_CELLULAR_LOCATION',
     *   'entity_src_gen',
     *   'host_org_cellular_location')
       call write_csdet('EXPRESSION_SYSTEM_VECTOR_TYPE',
     *     'entity_src_gen',
     *   'host_org_vector_type')
       call write_csdet('EXPRESSION_SYSTEM_VECTOR','entity_src_gen',
     *   'host_org_vector')
       call write_csdet('EXPRESSION_SYSTEM_PLASMID','entity_src_gen',
     *   'plasmid_name')
       call write_csdet('EXPRESSION_SYSTEM_GENE','entity_src_gen',
     *   'host_org_gene')
       call write_csdet('OTHER_DETAILS','entity_src_gen',
     *   'description')
       if (loop_) go to 300

       
600    if (tcsbuf(1:1).eq.char(0)) return
       kl = max(1,nblen(tcsbuf))
       if (tcsbuf(kl:kl).eq.';') tcsbuf(kl:kl) = ' '


       return
       end
       subroutine proc_titlrm(recstr,reccat,thisid,wcidstr)
C
C      process TITLE for JRNL or REMARK 1
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       character*12 recstr
       character*5 wcidstr
       character*7 thisid
       character*4 reccat
       character*8 mycat
       logical result_title,result_cid
       common/proc_titlrm_comment/result_title,result_cid
       character*2048 titstr,ntitstr,tail
       character*100 title,blanks
       data blanks/
     *   '                                                   '/
       icont = 0
       ipos = 1
       iepos = 51
       ltail = 0
       titstr = ' '
       title = ' '
       if (wide_code .eq. "yes") then
         iepos = 100
       endif
       text_ = .false.
       mycat = 'citation'
       result_title= pdbstr(mycat,
     * 'title',2048,ntitstr)
       if(.not.result_title) return
 100   continue
       titstr=ntitstr
 110   if(text_) then
         result_title= pdbstr(mycat,
     *   'title',2048,ntitstr)
         if(.not.result_title) ntitstr = ' '
       else
         ntitstr = ' '
         result_title=.false.
       endif
       if (result_title .and. ntitstr.eq. ' ') go to 110
       lw = nblen(titstr)
       if (type_code.eq.'u') titstr = tbxxupcs(titstr)
       call csplitstr(nword,titstr,cstr,lcstr,1024,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:lcstr(ii))//" ")
         lcstr(ii) = nblen(cstr(ii))
       enddo
       endif
       do ii = 1,nword
 200     lw = lcstr(ii)
         if (lw+ipos.le.iepos+1) then
           title(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if(ipos.eq.1) then
             title(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           title(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if(ipos.le.iepos) then
           title(ipos:ipos) = " "
           ipos = ipos+1
         endif
         if (ipos.gt.iepos.or.ltail.gt.0) then
           call write_remark(recstr(7:10),title(1:iepos),reccat,
     *       wcidstr,icont)
           title = ' '
           ipos = 1
           if (ltail.gt.0) then
             cstr(ii) = tail
             lcstr(ii) = ltail
             ltail = 0
             go to 200
           endif
         endif
       enddo
       if (ntitstr.ne.' ') go to 100
       if (ipos.gt.1) then               
         call write_remark(recstr(7:10),title(1:iepos),reccat,
     *    wcidstr, icont)
       endif
       return
       end
       subroutine proc_refrm(recstr,reccat,thisid,wcidstr)
C
C      process REF, PUBL and REFN for JRNL or REMARK 1
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       character*12 recstr
       character*7 thisid
       character*4 reccat
       character*8 mycat
       character*5 wcidstr
       logical result_title,result_cid
       logical result_vol, result_page, result_abbrev
       logical result_year, result_publ, result_ASTM
       logical result_country, result_ISSN, result_ISBN
       logical result_CSD, book_title
       character volstr*4,pagestr*5,
     *   yearstr*4,ASTMstr*6,countrystr*2,
     *   ISSNstr*25,ISBNstr*25,CSDstr*4
       character*2048 titstr,ntitstr,tail
       character*100 xtemp
       character*2 xvstr
       character*4 xISxNstr
       character*4 xASTMstr
       character*74 title,blanks
       common/proc_refrm_common/result_title,result_cid,
     *  result_vol, result_page, result_abbrev,
     *  result_year, result_publ, result_ASTM,
     *  result_country, result_ISSN, result_ISBN,
     *  result_CSD
       data blanks/
     *   '                                                   '/
       mycat='citation'
       result_vol=pdbstr(mycat,'journal_volume',4,volstr)
       result_page=pdbstr(mycat,'page_first',5,pagestr)
       result_year=pdbstr(mycat,'year',4,yearstr)
       result_country=pdbstr(mycat,'country',2,countrystr)
       result_ASTM=pdbstr(mycat,'journal_id_ASTM',
     *   6,ASTMstr)
       if (.not.result_ASTM)
     * result_ASTM=pdbstr(mycat,'journal_coden_ASTM',
     *   6,ASTMstr)
       result_ISSN=pdbstr(mycat,'journal_id_ISSN',
     *   25,ISSNstr)
       if (.not.result_ISSN)
     * result_ISSN=pdbstr(mycat,'journal_coden_ISSN',
     *   25,ISSNstr)
       result_ISBN=pdbstr(mycat,'book_id_ISBN',
     *   25,ISBNstr)
       if (.not.result_ISBN)
     *  result_ISBN=pdbstr(mycat,'book_coden_ISBN',
     *   25,ISBNstr)
       result_CSD=pdbstr(mycat,'journal_id_CSD',
     *   4,CSDstr)
       if (.not.result_CSD)
     *  result_CSD=pdbstr(mycat,'journal_coden_CSD',
     *   4,CSDstr)
       xASTMstr='ASTM'
       if(ASTMstr.eq.' ' .or.(.not.result_ASTM))
     *   xASTMstr=' '
       call rjust(volstr)
       xvstr = 'V.'
       if(volstr.eq.' '.or.(.not.result_vol)) xvstr='  '
       xISxNstr = '    '
       if(result_ISSN.and.ISSNstr.ne.' ') then
         xISxNstr = 'ISSN'
       else
         if(result_ISBN.and.ISBNstr.ne.' ') then
           xISxNstr = 'ISBN'
           ISSNstr=ISBNstr
         endif
       endif
       call rjust(pagestr)
       call rjust(yearstr)
       call rjust(CSDstr)
       do ipass = 1,2
       icont = 0
       ipos = 1
       iepos = 51
       if(ipass.eq.1) iepos = 28
       if (wide_code .eq. "yes") then
         iepos = 74
       endif
       ltail = 0
       titstr = ' '
       text_ = .false.
       mycat = 'citation'
       book_title = .false.
       if(ipass.eq.1) then
         result_title = pdbstr(mycat,
     *   'journal_abbrev',2048,ntitstr)
         if(.not.result_title) then
           book_title = pdbstr(mycat,
     *     'book_title',2048,ntitstr)
           result_title=book_title
         endif
       else
         result_title = pdbstr(mycat,
     *   'book_publisher',2048,ntitstr)
       endif
       if(.not.result_title) goto 900
       title=' '
 100   continue
       titstr=ntitstr
 110   if(text_) then
         if(ipass.eq.1) then
         if (.not.book_title) then
         result_title= pdbstr(mycat,
     *   'journal_abbrev',2048,ntitstr)
         else
         result_title= pdbstr(mycat,
     *   'book_title',2048,ntitstr)
         endif
         else
         result_title= pdbstr(mycat,
     *   'book_publisher',2048,ntitstr)
         endif
         if(.not.result_title) ntitstr = ' '
       else
         ntitstr = ' '
         result_title=.false.
       endif
       if (result_title .and. ntitstr.eq. ' ') go to 110
       lw = nblen(titstr)
       tmparg = titstr(1:max(1,lw))
       if (type_code.eq.'u') titstr = tbxxupcs(tmparg(1:max(1,lw)))
       call csplitstr(nword,titstr,cstr,lcstr,128,' ')
       if (type_code.eq.'p') then
       do ii = 1,nword
         cstr(ii) = typeset(cstr(ii)(1:lcstr(ii))//" ")
         lcstr(ii) = nblen(cstr(ii))
       enddo
       endif
       do ii = 1,nword
 200     lw = lcstr(ii)
         if (lw+ipos.le.iepos+1) then
           title(ipos:lw+ipos-1) = cstr(ii)(1:lw)
         else
           if(ipos.eq.1) then
             title(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
             ltail = lw - (iepos-ipos+1)
             ipos = iepos
             lw = lw-ltail
           else
             tail = cstr(ii)(1:lw)
             ltail = lw
             lw = 0
           endif
         endif
         if (ltail.gt.0.and.ipos+lw-1.lt.iepos) then
           title(ipos+lw:iepos) = blanks(1:iepos-(ipos+lw)+1)         
         endif
         ipos = ipos+lw
         if(ipos.le.iepos) then
           title(ipos:ipos) = " "
           ipos = ipos+1
         endif
         if (ipos.gt.iepos.or.ltail.gt.0) then               
           if (wide_code.eq."yes") then
             if (icont.eq.0 .and. ipass.eq.1) then
               write (xtemp,'(a74,1x,a2,a6,1x,a9,1x,a4)') 
     *           title(1:iepos),
     *           xvstr,volstr,pagestr,yearstr
             else
               xtemp = title(1:iepos)
             endif
             if (ipass.eq.1) then
               call write_remark(recstr(7:10),xtemp,'REF ',
     *           wcidstr,icont)
             else
               call write_remark(recstr(7:10),xtemp,'PUBL',
     *           wcidstr,icont)
             endif
           else 
           pdbline = pdbline+1
           if (recstr(1:4).eq.'REMA') numRemark = numRemark+1
           if (icont.eq.0) then
             if(ipass.eq.1) then
               write(iunpdb,
     *           '(a12,a4,3x,a28,2x,a2,a4,1x,a5,1x,a4)')
     *           recstr,'REF ',title(1:iepos),
     *           xvstr,volstr,pagestr,yearstr
             else
               write(iunpdb,'(a12,a4,3x,a)')
     *           recstr,'PUBL',title(1:iepos)
             endif
           else
             if(ipass.eq.1) then
               write(iunpdb,
     *         '(a12,a4,i2,1x,a)')
     *         recstr,'REF ',icont+1,title(1:iepos)
             else
               write(iunpdb,
     *         '(a12,a4,i2,1x,a)')
     *         recstr,'PUBL',icont+1,title(1:iepos)
             endif
           endif
           endif
           title = ' '
           ipos = 1
           icont = icont+1
           if (ltail.gt.0) then
             cstr(ii) = tail
             lcstr(ii) = ltail
             ltail = 0
             go to 200
           endif
         endif
       enddo
       if (ntitstr.ne.' ') go to 100
       if (ipos.gt.1) then               
         if (wide_code.eq."yes") then
             if (icont.eq.0 .and. ipass.eq.1) then
               write (xtemp,'(a74,1x,a2,1x,a6,1x,a9,2x,a4)') 
     *           title(1:iepos),
     *           xvstr,volstr,pagestr,yearstr
             else
               xtemp = title(1:iepos)
             endif
             if (ipass.eq.1) then
               call write_remark(recstr(7:10),xtemp,'REF ',
     *           wcidstr,icont)
             else
               call write_remark(recstr(7:10),xtemp,'PUBL',
     *           wcidstr,icont)
             endif
         else 
         pdbline = pdbline+1
         if (recstr(1:4).eq.'REMA') numRemark = numRemark+1
         if (icont.eq.0) then
           if(ipass.eq.1) then
             write(iunpdb,
     *         '(a12,a4,3x,a28,2x,a2,a4,1x,a5,1x,a4)')
     *         recstr,'REF ',title(1:iepos),
     *         xvstr,volstr,pagestr,yearstr
           else
             write(iunpdb,'(a12,a4,3x,a)')
     *         recstr,'PUBL',title(1:iepos)
           endif
         else
           if(ipass.eq.1) then
             write(iunpdb,
     *       '(a12,a4,i2,1x,a)')
     *       recstr,'REF ',icont+1,title(1:iepos)
           else
             write(iunpdb,
     *       '(a12,a4,i2,1x,a)')
     *       recstr,'PUBL',icont+1,title(1:iepos)
           endif
         endif
         endif
       endif
 900   if(ipass.eq.2) then
         pdbline = pdbline+1
         if (recstr(1:4).eq.'REMA') numRemark = numRemark+1
         if (wide_code.eq."yes") then
         call rjust(CSDstr)
         write(iunpdb,
     * '(a6,1x,a4,1x,a5,4x,1x,a4,1x,a4,1x,a6,1x,a2,1x,a4,1x,a25,1x,a6)')
     *   recstr(1:6),recstr(7:10),wcidstr,'REFN',xASTMstr,
     *   ASTMstr,countrystr,
     *   xISxNstr,ISSNstr,CSDstr
         else
         write(iunpdb,
     *   '(a12,a4,3x,a4,1x,a6,2x,a2,1x,a4,1x,a25,1x,a4)')
     *   recstr,'REFN',xASTMstr,ASTMstr,countrystr,
     *   xISxNstr,ISSNstr,CSDstr
         endif
       endif
       enddo
       return
       end
       subroutine proc_seqres
C
C      process SEQRES
C
       include 'cif2pdb.cmn'
       external pdbreal
       external pdbint
       external typeset
       logical result_entity_id
       logical result_num
       logical result_mon_id
       character*10 monidstr
       character*10 mychain
       character*33 entstr,pentstr
       character*9 snumstr,psnumstr
       character*90 resstr(500)
       integer iwresstr,iwmonid,iwsnum
       
       pdbrec = 'SEQRES'
       ipos = 1
       iser = 1
       numres = 0
       loop_ = .false.
       psnumstr = ' '
       pentstr = ' '
       resstr(1) = ' '
       text_=.false.
       iwresstr = 51
       iwmonid = 3
       iwsnum = 7
       if (wide_code .eq. "yes") then
         iwresstr = 90
         iwmonid = 10
         iwsnum = 9
       endif
 100   result_mon_id = pdbstr('entity_poly_seq','mon_id',
     *  iwmonid,monidstr)
       result_entity_id = pdbstr('entity_poly_seq','entity_id',
     *  33,entstr)
       result_num = pdbstr('entity_poly_seq','num',iwsnum,snumstr)
       if (result_mon_id .and. result_entity_id .and.
     *   result_num .and. snumstr .ne.psnumstr) then
         if (entstr.ne.pentstr )then
           if (ipos.eq.1) iser = iser-1
           if (numres.gt.0) then
           numhit = 0
           do ii = 1,num_chains
             if (entity_list(ii).eq.pentstr) then
               numhit = numhit+1
               mychain = chain_list(ii)
               call rjust(mychain)
               if (wide_code .eq. "yes" ) then
                 write(iunpdb,'(6HSEQRES,2x,i9,1x,3x,1x,a10,i9,a90)')
     *           (kser,
     *           mychain,numres,resstr(kser),kser=1,iser)
               else
                 write(iunpdb,'(6HSEQRES,2x,i2,1x,a1,1x,i4,2x,a51)')
     *           (kser,
     *           chain_list(ii),numres,resstr(kser),kser=1,iser)
               endif
               numSeq = numSeq+1
             endif
           enddo
           if (numhit.eq.0)
     *       call c2pwarn(' Found no chains for SEQRES entity '//
     *       entstr)
           endif
           pentstr = entstr
           numres = 0
           iser = 1
           ipos = 1
           resstr(1) = ' '
         endif
         psnumstr = snumstr
         if (entstr.eq.pentstr) then
           call rjust(monidstr(1:iwmonid))
           resstr(iser)(ipos:ipos+iwmonid-1) = monidstr(1:iwmonid)
           numres = numres+1
           ipos = ipos+iwmonid
           if (wide_code.ne."yes") then
             resstr(iser)(ipos+iwmonid:ipos+iwmonid) = ' '
             ipos = ipos+1
           endif
           if (ipos.gt.iwresstr) then
             iser = iser+1
             ipos = 1
             if (iser.gt.500) call c2perr(' Too many SEQRES cards')
             resstr(iser) = ' '
           endif
         endif
       endif
       if (loop_) go to 100
       if (ipos.eq.1) iser = iser-1
       if (numres.gt.0) then
         numhit = 0
         do ii = 1,num_chains
           if (entity_list(ii).eq.pentstr) then
             numhit = numhit+1
             if (wide_code .eq. "yes" ) then
               mychain = chain_list(ii)
               call rjust(mychain)
               write(iunpdb,'(6HSEQRES,2x,i9,1x,3x,1X,a10,i9,a90)')
     *           (kser,
     *           mychain,numres,resstr(kser),kser=1,iser)
             else
               write(iunpdb,'(6HSEQRES,2x,i2,1x,a1,1x,i4,2x,a51)')
     *         (kser,
     *         chain_list(ii),numres,resstr(kser),kser=1,iser)
             endif
             numSeq = numSeq+1
           endif
         enddo
         if (numhit.eq.0)
     *     call c2pwarn(' Found no chains for SEQRES entity '//
     *     entstr)
       endif
       return
       end
       subroutine write_title(title,icont,iepos)
C
C      Write TITLE, COMPND, SOURCE or CAVEAT record
C      always write a TITLE record if none has been written
C
C      This code does not write immediately.  It stores
C      the next line to be written in tcsbuf, first writing
C      out tcsbuf if tcsbuf(1:1) is not char(0).  To force
C      out the pending tcsbuf, call with icont = -1
C
       include 'cif2pdb.cmn'
       external pdbstr
       external pdbreal
       external pdbint
       external typeset
       character *(*) title
       integer icont, ktcs, kl, kll, ipos, iopos
       integer i, ik
       character *110 titout
       
       if (tcsbuf(1:1) .ne. char(0)) then
         ktcs = max(1,nblen(tcsbuf))
         write(iunpdb,'(a)') tcsbuf(1:ktcs)
         tcsbuf = char(0)
       endif

       if (icont.eq.-1) return
       kl = max(1,nblen(title))
       kll = min(iepos-4,kl)
       if (numTitle.eq.0.and.pdbrec.ne.'TITLE') then
         pdbline = pdbline+1
         write(tcsbuf,sfmt_orig)
     *     "TITLE ","TITLE not provided, see COMPND and SOURCE"
         numTitle = numTitle+1
       endif
       
       iopos = 1

 100   continue
       ipos = 1
       titout(1:iepos) = ' '
       if (wide_code.ne."yes" .and. icont.gt.0 .and. 
     *   title(iopos:iopos).ne.' ') ipos=2
       if (kl - iopos +1 .le. iepos-ipos+1) then
         titout(ipos:ipos+kl-iopos) = title(iopos:kl)
         iopos = kl+1
       else
         do i = iopos,min(kl,iopos+iepos-ipos)
         ik = min(kl,iopos+iepos-ipos)+iopos-i+1
         if (title(ik:ik).eq.' ') go to 110
         enddo
         titout(ipos:ipos+min(kl,iopos+iepos-ipos)-iopos) = 
     *     title(iopos:min(kl,iopos+iepos-ipos))
         iopos = min(kl,iopos+iepos-ipos)+1
         go to 120
 110     titout(ipos:ipos+ik-iopos) = title(iopos:ik)
         iopos = ik+1
       endif
 120   continue       
       if (tcsbuf(1:1) .ne. char(0)) then
         ktcs = max(1,nblen(tcsbuf))
         write(iunpdb,'(a)') tcsbuf(1:ktcs)
         tcsbuf = char(0)
       endif
       pdbline = pdbline+1
       kll = max(1,nblen(titout(1:iepos)))
       if (icont.eq.0)
     *   write(tcsbuf,
     *     sfmt_orig)
     *     pdbrec,titout(1:kll)
       if (icont.gt.0)
     *   write(tcsbuf,
     *     sfmt_cont)
     *     pdbrec,icont+1,titout(1:kll)
       if (pdbrec .eq. 'TITLE') numTitle = numTitle+1
       if (pdbrec .eq. 'COMPND') numCompound = numCompound+1
       if (pdbrec .eq. 'SOURCE') numSource = numSource+1
       if (pdbrec .eq. 'CAVEAT') numCaveat = numCaveat+1
       icont = icont+1
       if (iopos.le.kl) go to 100
       return
       end
       
       subroutine proc_title
C
C      process TITLE
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       logical result_title,result_cavid
       logical bkmrk_
       character*2048 titstr,tail
       character*110 title,blanks
       character*7 cidstr
       logical result_bk, result_cid
       integer marktl
       integer kcpos, kl
       common/proc_title_common/result_title,result_cavid
     *  result_bk,result_cid
       data blanks/
     *   '                                                            '/
       pdbrec = 'TITLE'
       icont = numTitle
       title = ' '
       ipos = 1
       iepos = 60
       ltail = 0
       mode = 1
       text_ = .false.
       kcpos = 2
       marktl = 0
       if (wide_code .eq. "yes") then
         kcpos = 1
         iepos = 110
       endif
       result_bk = bkmrk_(marktl)
 100   continue
C      mode           CIF tag
C      -1             _citation.title
C      0,1            _struct.title
C      2              _database_pdb_caveat.text
C
C      We start in mode 1.  If we leave the TITLE record
C      without find any TITLE text, we divert to mode -1
C      trying for a citation title.  If we fail at that
C      try for the COMPOUND text and use that for the TITLE
C      as mode 0. If we fail at that we put in an explanatory note
C      When we are done we use the bookmark in marktl
C      to return to the original mode 1 processing
C
       if (mode.eq.-1) then
         result_title = pdbstr('citation',
     *     'title',2048,titstr)
       endif
       if (mode.eq.1.or.mode.eq.0) then
         result_title= pdbstr('struct',
     *     'title',2048,titstr)
       endif
       if (mode.eq.2) then
         result_title = pdbstr('database_pdb_caveat',
     *     'text',2048,titstr)
         if (titstr(1:6).eq.'CAVEAT' .and.
     *     titstr(12:15).eq.myentry) then
           titstr(1:15) = '               '
         endif
       endif
       if (type_code.eq.'u') titstr = tbxxupcs(titstr)
       call csplitstr(nword,titstr,cstr,lcstr,1024,' ')
       if (type_code.eq.'p') then
         do ii = 1,nword
           cstr(ii) = typeset(cstr(ii)(1:lcstr(ii))//" ")
           lcstr(ii) = nblen(cstr(ii))
         enddo
       endif       
C
       do ii = 1,nword
         lw = max(1,lcstr(ii))
         if (ii.eq.1 .and. cstr(1)(lw:lw).eq.":") then
           if ((ipos.gt.1.and.icont.eq.0)
     *     .or.(ipos.gt.kcpos.and.icont.gt.0)) then 
             call write_title(title(1:ipos-1),icont,iepos)
             title = ' '
             ipos = kcpos          
           endif
         endif
         if (cstr(ii)(1:lcstr(ii)).eq.'Compound::' .or.
     *   cstr(ii)(1:lcstr(ii)).eq.'Source::' .or.
     *   cstr(ii)(1:lcstr(ii)).eq.'COMPOUND::' .or.
     *   cstr(ii)(1:lcstr(ii)).eq.'SOURCE::' .or.
     *   cstr(ii)(1:lcstr(ii)).eq.'Error::' .or.
     *   cstr(ii)(1:lcstr(ii)).eq.'ERROR::' ) then
           if (cstr(ii)(1:lcstr(ii)).eq.'Compound::' .or.
     *       cstr(ii)(1:lcstr(ii)).eq.'COMPOUND::') then
             pdbrec = 'COMPND'
             icont = numCompound
           endif
           if (cstr(ii)(1:lcstr(ii)).eq.'Source::' .or.
     *       cstr(ii)(1:lcstr(ii)).eq.'SOURCE::') then
             pdbrec = 'SOURCE'
             icont = numSource
           endif
           if (cstr(ii)(1:lcstr(ii)).eq.'Error::' .or.
     *       cstr(ii)(1:lcstr(ii)).eq.'ERROR::') then
             pdbrec = 'CAVEAT'
             icont = numCaveat
           endif
           if (pdbrec.ne.'TITLE' .and.
     *       mode.eq.0 .and.result_bk .and. numTitle.eq.0) then
             icont = 0
             ipos = 1
             pdbrec = 'TITLE '
             go to 100
           endif
           if (pdbrec.ne.'COMPND' .and. pdbrec.ne.'TITLE' .and.
     *       mode.eq.0 .and. result_bk) then
             result_bk = bkmrk_(marktl)
             result_bk = .false.
             mode = 1
             pdbrec = 'TITLE '
             loop_ = .false.
             text_ = .false.
             icont = numTitle
             ltail = 0
             ipos = 1
             if (icont.gt.0) ipos = kcpos
             go to 100
           endif
           if (pdbrec.ne.'TITLE ' .and. 
     *       numTitle.eq.0 .and.
     *       mode .eq. 1 .and. result_bk) then
             result_bk = bkmrk_(marktl)
             result_bk = .false.
             mode = -1
             pdbrec = 'TITLE '
             loop_ = .false.
             icont = 0
             ltail = 0
             ipos = 1
             if (icont.gt.0) ipos = kcpos
 110         result_cid = pdbstr('citation','id',7,cidstr)
             if (.not.result_cid) go to 150
             if (cidstr.ne.'primary') then
               if (loop_) go to 110
               go to 150
             endif
             text_ = .false.
             go to 100
C
C            The attempt at getting the title from
C            the citation failed
C
 150         mark_tl = 0
             result_bk = bkmrk_(mark_tl)
             mode = 0
             loop_ = .false.
             text_ = .false.
             go to 100
           endif
           ipos = 1
           if (icont.gt.0) ipos = kcpos
           if (pdbrec .eq. 'CAVEAT') then
             if (wide_code.eq."yes") then
               ipos = 17
               title(1:15) = myentry
             else
               ipos = 10
               title(1:10) = ' '//myentry
             endif
             go to 100
           endif
           go to 300
         endif
 200   lw = lcstr(ii)
       if (lw+ipos.le.iepos+1) then
         title(ipos:lw+ipos-1) = cstr(ii)(1:lw)
       else
         if((ipos.eq.1.and.icont.eq.0)
     *     .or.(ipos.le.kcpos.and.icont.gt.0)) then
           title(ipos:iepos) = cstr(ii)(1:iepos-ipos+1)
           ltail = lw - (iepos-ipos+1)
           tail = cstr(ii)(iepos-ipos:lw)
           lw = lw-ltail
         else
           tail = cstr(ii)(1:lw)
           ltail = lw
           lw = 0
         endif
       endif
       ipos = ipos+lw
       if(ipos.le.iepos) then
         title(ipos:ipos) = " "
         ipos = ipos+1
       endif
       if ((ipos.gt.iepos.or.ltail.gt.0) .or.
     *     (ii.eq.nword.and.lw.ge.1.and.
     *     cstr(ii)(max(1,lw):max(1,lw)).eq.";")) then
           kl = max(1,nblen(title(1:ipos-1)))
           call write_title(title(1:kl),icont,iepos)        
           title = ' '
           ipos = kcpos
           if (ltail.gt.0) then
             cstr(ii) = tail
             lcstr(ii) = ltail
             ltail = 0
             go to 200
           endif
       endif
 300   continue
       enddo
       if (text_) go to 100
       if ((pdbrec.ne.'CAVEAT' .and.
     *   ((ipos.gt.1.and.icont.eq.0)
     *   .or.(ipos.gt.kcpos.and.icont.gt.0))) .or.
     *   (pdbrec.eq.'CAVEAT' .and. ipos.gt.10)) then 
         call write_title(title(1:ipos-1),icont,iepos)  
         title = ' '
         ipos = kcpos
       endif
       if (mode.eq.2) then
         if (result_bk) result_bk = bkmrk_(marktl)
         if (numSource.eq.0) call proc_source
         call write_title(" ",-1,-1)
         return
       endif
       titstr = ' '
       if (mode.eq. -1) then
         if (numTitle.eq.0) then
           mode = 0
           ipos = 1
           go to 100
         else
           mode = 1
           ipos = 1
           go to 100
         endif
       endif
       if (mode.eq.0) then
         mode = 1
         ipos = 1
         go to 100
       endif
       result_cavid = pdbstr('database_pdb_caveat',
     * 'id',2048,titstr)
       if (result_cavid .and. titstr(1:4).ne.'    ') then
         mode = 2
         pdbrec = 'CAVEAT'
         icont = numCaveat
         if (wide_code .eq. "yes") then
           ipos = 17
           title(1:15) = myentry
         else
           ipos = 10
           title(1:10) = ' '//myentry
         endif
         go to 100
       endif
       if (numSource.eq.0) call proc_source
       call write_title(" ",-1,-1)
       return
       end
       function pdbstr(curcat,curnam,curlim,gotnam)
C
C      extracts a value for token "_"//curcat//"."//curnam
C      enforcing a limit of curlim characters
C
C      The strange code on result is to avoid rejection of a
C      numeric field if we asked for it as a character field
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbreal
       external typeset
       logical char_
       character*(*) curcat
       character*(*) curnam
       integer curlim,ll
       character*(*) gotnam
       character*2048 temp
       gotnam = ' '
       ll = 1
       if (.not.text_) strg_ = " "
       tmparg = "_"//curcat//"."//curnam
       result = char_(tmparg(1:nblen(tmparg)),temp)
       if (result.or.(type_.eq.'numb')) then
         ll = nblen(temp)
         if (ll.ge.curlim) then
           kpmap = nmap
           tmparg = curnam//":"//temp(1:ll)//char(0)
           call hash_store(tmparg,mapstr,mapchain,
     *       NUMSTR,nmap,mhash,NUMHSH,ifrom)
           if (ifrom.eq.0)
     *       call c2perr(' More than NUMSTR strings mapped ')
           if (ifrom.eq.kpmap+1) mapto(ifrom) = 0
           kfrom = mapto(ifrom)
           if (kfrom.ne.0) then
             temp = mapstr(kfrom)
             ll = nblen(temp)
           endif
           if (ll .gt. curlim) then
             call c2pwarn(' String '//temp(1:ll)//' too long'//
     *       ' truncating to '//temp(1:curlim))
             ll = curlim
             if (kfrom.eq.0) then
               call hash_store(temp(1:ll)//char(0),mapstr,mapchain,
     *           NUMSTR,nmap,mhash,NUMHSH,ito)
               if (ito.eq.0)
     *           call c2perr(' More than NUMSTR strings mapped ')
               mapto(ifrom)=ito
             endif
           endif
         else
           if(ll.eq.0) then
             ll=1
             temp(1:1) = ' '
           endif
         endif
         gotnam = temp(1:ll)
       endif
       pdbstr = result.or.(type_ .eq. 'numb')
       return
       end
       function pdbint(curcat,curnam,curlim,gotint,intval)
C
C      extracts a value for token "_"//curcat//"."//curnam
C      required to be a right justified integer in a field
C      curlim characters wide
C
C
C      The strange code on result is to avoid rejection of a
C      numeric field if we asked for it as a character field
C
       include 'cif2pdb.cmn'
       external pdbreal
       external pdbstr
       external typeset
       logical char_
       character*(*) curcat
       character*(*) curnam
       integer curlim,ll
       character*(*) gotint
       character*20 stars,blanks
       character*2048 temp
       character*80 temp2
       data stars/'********************'/
       data blanks/'                    '/
       gotint = ' '
       intval = 0
       strg_ = " "
       tmparg = "_"//curcat//"."//curnam
       result = char_(tmparg,temp)
       if ((result.or.(type_.ne.'null')).and.
     *   strg_(1:long_).ne." ") then
         ll = nblen(temp)
         if (ll .gt. curlim) then
           call c2pwarn(' found '//temp(1:ll)//
     *   ' too large, converted to'//stars(1:curlim))
           ll = curlim
           temp(1:curlim) = stars(1:curlim)
         endif
         temp2 = temp(1:80)
         call rjust(temp2)
         read(temp2,'(i80)',err=105) intval
         go to 110
105      tmparg=' non-numeric field '//
     *     "_"//curcat//"."//curnam//" = "//temp2(min(80,81-ll):80)
         call c2pwarn(tmparg(1:nblen(tmparg)))
         intval=0
110      gotint = temp2(min(80,81-ll):80)
         if (ll.eq.0 .or. 
     *     (ll.eq.1 .and. type_.eq.'null')) gotint = ' '
       endif
       pdbint = (result.or.(type_ .ne. 'null'))
       return
       end
       function pdbreal(curcat,curnam,curlim,gotreal,realval)
C
C      extracts a value for token "_"//curcat//"."//curnam
C      required to be a right justified integer in a field
C      curlim characters wide
C
C
C      The strange code on result is to avoid rejection of a
C      numeric field if we asked for it as a character field
C
       include 'cif2pdb.cmn'
       external pdbint
       external pdbstr
       external typeset
       logical numd_
       character*(*) curcat
       character*(*) curnam
       integer curlim,ll
       character*(*) gotreal
       real*8 xnumb,xesd
       real*8 realval
       character*80 temp
       character*20 stars,blanks
       data stars/'********************'/
       data blanks/'                    '/
       gotreal = ' '
       realval = 0.
       strg_ = " "
       tmparg = "_"//curcat//"."//curnam
       xesd = 0.
       result = numd_(tmparg,xnumb,xesd)
       temp = strg_(1:long_)
       if ((result.or.(type_.ne.'null')) .and.
     *   strg_ .ne. " ") then
       ll = nblen(temp)
       if (xesd .ne.0.) then
          ll = index(temp,'(') -1
       endif
       if (ll .gt. curlim) then
         call c2pwarn(' found '//temp(1:ll)//
     *     ' too large, converted to'//stars(1:curlim))
         ll = curlim
         temp(1:curlim) = stars(1:curlim)
       endif
       if (type_.ne.'numb') then
         tmparg=' non-numeric field '//
     *     "_"//curcat//"."//curnam//" = "//temp(1:ll)
         call c2pwarn(tmparg)
       endif
       if (ll.lt.curlim.and.ll.gt.0) then
         temp = blanks(1:curlim-ll)//temp(1:ll)
         ll=curlim
       endif
       gotreal = temp(1:ll)
       if (ll.eq.1 .and. type_.eq.'null') gotreal = ' '
       if(type_.eq.'numb')
     *   realval = xnumb
       else
C         print *,' failed to find '//"_"//curcat//"."//curnam
       endif
       pdbreal = (result.or.(type_ .ne. 'null')).and.
     *   strg_ .ne." "
       return
       end
       function nblen(str)
C
C      variant of lastnb which also detects a null character to
C      terminate a string
C
       character*(*) str
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
       subroutine decjust(str,decpos)
C
C      justify string str to place the decimal point at decpos
C
       character*(*) str
       integer decpos
       integer i, ii, kd
       
       kl = len(str)
       kd = index(str,".")
       if (kd.eq.decpos) return
       if (kd.eq.0) then
         call rjust(str)
         return
       endif
       if (kd.lt.decpos) then
         do ii = 1,kl-(decpos-kd)
           i = kl-ii+1
           str(i:i) = str(i-(decpos-kd):i-(decpos-kd))
         enddo
         do ii = 1,decpos-kd
           str(ii:ii) = ' '
         enddo
       else
         do i = 1,kl-(kd-decpos)
           str(i:i) = str(i+(kd-decpos):i+(kd-decpos))
         enddo
         do i = kl-(kd-decpos),kl
           str(i:i) = ' '
         enddo       
       endif
       return
       end
       subroutine rjust(str)
C
C      right justify string str, in place
C
       character*(*) str
       character*2048 blanks
       character*2048 xxstr
       blanks = ' '
       ll = len(str)
       ln = nblen(str)
       if (ln.lt.ll .and. ln.gt.0) then
         xxstr(1:ln) = str(1:ln)
         str = blanks(1:ll-ln)//xxstr(1:ln)
       endif
       return
       end
       subroutine ljust(str)
C
C      left justify str in place
C
       character*(*) str
       ll = len(str)
       jj = 1
       do ii = 1,ll
       if (str(ii:ii).ne.' ') go to 100
       enddo
       return
 100   if (ii.eq.1) return
       do kk = ii,ll
         str(jj:jj)=str(kk:kk)
         jj=jj+1
       enddo
       if(jj.le.ll) str(jj:ll) = ' '     
       return
       end

         subroutine c2perr(mess)
C
C        variant of ciftbx err routine for cif2pdb
C
         character*(*) mess
         call c2pmsg('error',mess)
         stop
         end
         subroutine c2pwarn(mess)
C
C        variant of ciftbx warn for cif2pdb
C
         character*(*) mess
         call c2pmsg('warning',mess)
         return
         end
         subroutine c2pmsg(flag,mess)
C
         include   'cif2pdb.cmn'
         external pdbint
         external pdbreal
         external pdbstr
         external typeset
         character*(*) flag
         character*(*) mess
         character*(MAXBUF)  tline
         integer       ll,ls,ltry,ii,i
C
         tline = ' cif2pdb'
         if (nblen(flag).gt.0) 
     *     tline = tline//' '//flag(1:nblen(flag))
         tline= tline//': '
     *   //outent(1:loutent)//' '//pdbrec
     *   //' line:'
         ll = max(1,lastnb(tline))
         write(iunerr,'(a,i7)')tline(1:ll),pdbline
         ll=len(mess)
         ls=1
100      if(ll-ls.le.79) then
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
         
       function pdbdate(cifdate)
C
C      return a PDB format date (dd-mmm-yy) from a cif format date
C      yyyy-mm-dd
C
       character*9 pdbdate
       character*(*) cifdate
       character*3 months(12)
       character*256 tmparg
       data months/'JAN','FEB','MAR','APR','MAY','JUN','JUL',
     *  'AUG','SEP','OCT','NOV','DEC'/
       read(cifdate,'(i4,1x,i2,1x,i2)',err=100) iyyyy,imm,idd
       if (imm.lt.1.or.imm.gt.12.or.idd.lt.1.or.idd.gt.31)
     *  go to 100
       if (iyyyy.lt.1970.or.iyyyy.gt.2069) go to 100
       write(pdbdate,'(i2.2,1h-,a3,1h-,i2.2)')idd,months(imm),
     * mod(iyyyy,100)
       return
 100   tmparg = ' Unable to translate date: '//cifdate
       call c2pwarn(tmparg)
       pdbdate = '??-???-??'
       return
       end
C         
       function wpdbdate(cifdate)
C
C      return a WPDB format date (dd-mmm-yyyy) from a cif format date
C      yyyy-mm-dd
C
       character*11 wpdbdate
       character*(*) cifdate
       character*3 months(12)
       character*256 tmparg
       data months/'JAN','FEB','MAR','APR','MAY','JUN','JUL',
     *  'AUG','SEP','OCT','NOV','DEC'/
       read(cifdate,'(i4,1x,i2,1x,i2)',err=100) iyyyy,imm,idd
       if (imm.lt.1.or.imm.gt.12.or.idd.lt.1.or.idd.gt.31)
     *  go to 100
       write(wpdbdate,'(i2.2,1h-,a3,1h-,i4.4)')idd,months(imm),
     * iyyyy
       return
 100   tmparg = ' Unable to translate date: '//cifdate
       call c2pwarn(tmparg)
       wpdbdate = '??-???-????'
       return
       end
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
C
       subroutine csplitstr(numf,str,sarry,lsarry,maxf,fs)
C
C      split the string str into a maximum of maxf fields using the
C      field separator fs.  The number of fields found is reported in
C      numf.  The fields are split into the array sarry, with each
C      field terminated by a null character.  The numbers of
C      characters in each field is reported in the array lsarry.
C      An empty field is reported with a length of 0.
C
C      This version uses counted strings
C
C      The field separator ' ' is treated as a special case, with
C      full blank and tab removal, and only the count of non-blank
C      fields returned.
C
       integer numf, maxf
       character*(*) str, sarry(maxf), fs
       integer lsarry(maxf)
       ll = len(str)
       kl = index(str,char(0))
       if (kl.gt.0) ll = kl-1
       lfs = len(fs)
       if (fs.eq.' ') lfs = 1
       is = 1
       numf = 0
       if (ll.le.0) return
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
           lsarry(numf) = min(ll-is+1,len(sarry(numf))-1)
           sarry(numf)(1:lsarry(numf)+1) = str(is:ll)//char(0)
         endif
         return
       else
         if (fs.ne.' ' .or. jj.gt.1) then
         numf = numf+1
           if (numf.le.maxf) then
             if (jj.gt.1) then
               lsarry(numf) = min(jj-1,len(sarry(numf))-1)
               sarry(numf)(1:lsarry(numf)+1) 
     *           = str(is:is-1+jj-1)//char(0)
             else
               sarry(numf)(1:1) = char(0)
               lsarry(numf) = 0
             endif
           endif
         endif
         is = is+jj+lfs-1
         if (is.le.ll) go to 100
         if (fs.ne.' ') then 
           numf = numf+1
           if (numf.le.maxf) then
             sarry(numf)(1:1) = char(0)
             lsarry(numf) = 0
           endif
         endif
         return
       endif
       end
C
       function nounder(str)
C
C      convert a string with underscores into a string
C      with underscores replaced by blanks, except when
C      an underscore follows an underscore, in which case
C      the single underscore in included.
C
       character*255 nounder
       character*(*) str
       ll = nblen(str)
       lp = ll
       kp = 1
       ks = 1
       if (ll.gt.0) then
 100     kt = index(str(ks:ll),'_')
         if (kt.eq.0) then
           nounder(kp:lp+1) = str(ks:ll)//char(0)
           return
         else
           nounder(kp:kp+kt-1) = str(ks:ks+kt-1)
           kp=kp+kt
           ks=ks+kt
           if (ks.lt.ll) then
             if (str(ks:ks).eq.'_') then
               ks = ks+1
               lp = lp-1
             else
               nounder(kp-1:kp-1) = ' '
             endif
           else
             nounder(kp-1:kp-1) = ' '
           endif
         endif
         if (ks.le.ll) go to 100
       endif
       nounder(kp:kp) = char(0)
       return
       end
       function typeset(str)
C
C      convert a string with upper and lower case characters
C      into an upper case string with the following typesetting
C      conventions (from pre 1993 PDB format descriptions)
C
C      The string is treated as a blank-delimited set of words
C      where, unless flagged otherwise, the first character is
C      intended to be upper-case and the rest lower-case.  These
C      assumptions are over-ridden by the following codes:
C
C        blank  -- immediately following character is upper-case
C        comma  -- immediately following character is upper-case
C        period -- immediately following character is upper-case
C        left parenthesis
C               -- immediately following character is upper-case
C        asterisk
C               -- immediately following character is upper-case
C               and the asterisk is dropped from the output string
C               (therefore "*" is changed to "(STAR)")
C        slash  -- acts like a caps lock until the end of string,
C               a hyphen or a dollar sign, whichever comes first
C               (therefore "/" is changed to "(SLASH)")
C        dollar sign
C               -- terminates effect of slash, blank, comma period, etc.
C               (therefore "$" is changed to "(DOLLAR)")
C
       character*2048 typeset
       character*2148 tmparg
       character*(*) str
       character*1 mode,c
       logical upper,lower
       ll = nblen(str)
       kp = 1
       ks = 1
       mode = " "
       if (ll.gt.0) then
 100     c = str(ks:ks)
         ic = ichar(c)
         upper = (ic.ge.ichar('A')).and.(ic.le.ichar('Z'))
         lower = (ic.ge.ichar('a')).and.(ic.le.ichar('z'))
         if(upper) then
           if(mode.eq.'l') then
             typeset(kp:kp) = "*"
             mode = " "
             if(kp.gt.2) then
               if(typeset(kp-2:kp-2).eq."*") then
                 typeset(kp-2:kp-2) = "/"
                 mode = "/"
                 go to 200
               endif
             endif
             kp=kp+1
             go to 200
           endif
           go to 200
         endif
         if(lower) then
           if(mode.eq." ".or.mode.eq."/") then
             typeset(kp:kp) = "$"
             kp = kp+1
             mode = "l"
           endif
           go to 200
         endif
         if(c.eq."$") then
           typeset(kp:kp+7) = "(DOLLAR)"
           kp = kp+8
           mode = "l"
           go to 300
         endif
         if(c.eq."*") then
           typeset(kp:kp+5) = "(STAR)"
           kp = kp+6
           mode = "l"
           go to 300
         endif
         if(c.eq."/") then
           typeset(kp:kp+6) = "(SLASH)"
           kp = kp+7
           mode = "l"
           go to 300
         endif
         if(c.eq."-") then
           if(mode.eq."/") mode = "l"
         endif
 200     if(lower) c = char(ic+ichar('A')-ichar('a'))
         if (mode.ne."/".and.
     *     (c.eq." ".or.c.eq.",".or.c.eq."(".or.c.eq.".")) 
     *     mode = " "
         if (mode.eq."/".and.
     *     (c.eq." ".or.c.eq.",".or.c.eq."(".or.c.eq.".")) then
           typeset(kp:kp) = "$"
           kp = kp+1
           mode = " "
         endif
         if(mode.eq." " .and.
     *     c.ne." " .and. c.ne."," .and. c.ne.".") mode = "l"
         typeset(kp:kp) = c
         kp = kp+1
 300     ks = ks+1
         tmparg=' Unable to typeset '//str
         if (kp.gt.245) call c2perr(tmparg)
       if (ks.le.ll) go to 100
       endif
       typeset(kp:kp) = char(0)
       return
       end
      function findtag_(name)

      logical findtag_
      include 'ciftbx.sys'
      logical find_,data_
      character*(MAXBUF) tbxxlocs
      logical result
      character*(*) name
      character*(MAXBUF) string, temp
      integer ii, ifirst
      ifirst = 0
      temp = tbxxlocs(name)
 200  if (nname.gt.0) then
        call hash_find(temp,
     *    dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *    ii)
        if (ii.gt.0) then
          findtag_ = .true.
          return
        endif
      endif
      if (ifirst.eq.0) then
        result=find_(' ','head',string)
        ifirst = 1
      endif
      result = data_(' ')
      if (.not.result) then
        findtag_ = .false.
        return
      endif
      go to 200
      end

      subroutine tagchk(iremno,iunpdb,wide_code,pdbline,numRemark)
C
C     routine to check data block names aginst the list of tags
C     processed by cif2pdb and prints out all unprocessed tags
C     as part of REMARK iremno
C
      PARAMETER (NUMusdtgs=189)
      include 'ciftbx.sys'
      character*(*) wide_code
      integer pdbline,numRemark
      logical procnam(NUMBLOCK)
C
C     The array usdtgs must contain the names of each of the tags
C     processed by cif2pdb, given in lower case. 
C
      character*(NUMCHAR) usdtgs(NUMusdtgs)
      character*(MAXBUF) string,temp
      character*2048 tmparg
      character*4 xdtp
      character*22 rempref
      logical result
      logical ploop_
      logical pdata_
      logical prefx_
      logical test_
      logical pfile_
      logical char_
      logical ptext_
      logical pchar_
      logical numd_
      logical pnumd_
      logical find_
      logical data_
      logical dtype_
      integer nolin
      double precision dpn, dpesd
      data usdtgs /
     * '_atom_site.b_iso_or_equiv', '_atom_site.u_iso_or_equiv',
     * '_atom_site.b_iso_or_equiv_esd', 
     * '_atom_site.u_iso_or_equiv_esd',
     * '_atom_site.aniso_b[1][1]', '_atom_site.aniso_b[1][2]',
     * '_atom_site.aniso_b[1][3]', '_atom_site.aniso_b[2][2]',
     * '_atom_site.aniso_b[2][3]', '_atom_site.aniso_b[3][3]',
     * '_atom_site.aniso_b[1][1]_esd', '_atom_site.aniso_b[1][2]_esd',
     * '_atom_site.aniso_b[1][3]_esd', '_atom_site.aniso_b[2][2]_esd',
     * '_atom_site.aniso_b[2][3]_esd', '_atom_site.aniso_b[3][3]_esd',
     * '_atom_site.aniso_u[1][1]', '_atom_site.aniso_u[1][2]',
     * '_atom_site.aniso_u[1][3]', '_atom_site.aniso_u[2][2]',
     * '_atom_site.aniso_u[2][3]', '_atom_site.aniso_u[3][3]',
     * '_atom_site.aniso_u[1][1]_esd', '_atom_site.aniso_u[1][2]_esd',
     * '_atom_site.aniso_u[1][3]_esd', '_atom_site.aniso_u[2][2]_esd',
     * '_atom_site.aniso_u[2][3]_esd', '_atom_site.aniso_u[3][3]_esd',
     * '_atom_site.auth_atom_id', '_atom_site.auth_comp_id',
     * '_atom_site.auth_asym_id', '_atom_site.auth_seq_id',
     * '_atom_site.cartn_x', '_atom_site.cartn_y', '_atom_site.cartn_z',
     * '_atom_site.fract_x', '_atom_site.fract_y', '_atom_site.fract_z',
     * '_atom_site.cartn_x_esd', '_atom_site.cartn_y_esd',
     * '_atom_site.cartn_z_esd',
     * '_atom_site.fract_x_esd', '_atom_site.fract_y_esd',
     * '_atom_site.fract_z_esd',
     * '_atom_site.footnote_id',
     * '_atom_site.group_pdb',
     * '_atom_site.id',
     * '_atom_site.label_alt_id', '_atom_site.label_asym_id',
     * '_atom_site.label_atom_id', '_atom_site.label_comp_id',
     * '_atom_site.label_entity_id', '_atom_site.label_seq_id',
     * '_atom_site.label_model_id', '_atom_site.occupancy',
     * '_atom_site.occupancy_esd',
     * '_atom_site.pdb2cif_label_model_id',
     * '_atom_site.pdbx_pdb_model_num',
     * '_atom_site.type_symbol',
     * '_atom_site_anisotrop.id','_atom_site_anisotrop.type_symbol',
     * '_atom_site_anisotrop.pdbx_auth_atom_id',
     * '_atom_site_anisotrop.pdbx_auth_alt_id',
     * '_atom_site_anisotrop.pdbx_auth_comp_id',
     * '_atom_site_anisotrop.pdbx_auth_asym_id',
     * '_atom_site_anisotrop.pdbx_auth_seq_id',
     * '_atom_site_anisotrop.pdbx_label_atom_id',
     * '_atom_site_anisotrop.pdbx_label_alt_id',
     * '_atom_site_anisotrop.pdbx_label_comp_id',
     * '_atom_site_anisotrop.pdbx_label_asym_id',
     * '_atom_site_anisotrop.pdbx_label_seq_id',
     * '_atom_site_anisotrop.b[1][1]','_atom_site_anisotrop.b[1][2]',
     * '_atom_site_anisotrop.b[1][3]',
     * '_atom_site_anisotrop.b[2][1]','_atom_site_anisotrop.b[2][2]',
     * '_atom_site_anisotrop.b[2][3]',
     * '_atom_site_anisotrop.b[3][1]','_atom_site_anisotrop.b[3][2]',
     * '_atom_site_anisotrop.b[3][3]',
     * '_atom_site_anisotrop.u[1][1]','_atom_site_anisotrop.u[1][2]',
     * '_atom_site_anisotrop.u[1][3]',
     * '_atom_site_anisotrop.u[2][1]','_atom_site_anisotrop.u[2][2]',
     * '_atom_site_anisotrop.u[2][3]',
     * '_atom_site_anisotrop.u[3][1]','_atom_site_anisotrop.u[3][2]',
     * '_atom_site_anisotrop.u[3][3]',
     * '_atom_site_anisotrop.b[1][1]_esd',
     * '_atom_site_anisotrop.b[1][2]_esd',
     * '_atom_site_anisotrop.b[1][3]_esd',
     * '_atom_site_anisotrop.b[2][1]_esd',
     * '_atom_site_anisotrop.b[2][2]_esd',
     * '_atom_site_anisotrop.b[2][3]_esd',
     * '_atom_site_anisotrop.b[3][1]_esd',
     * '_atom_site_anisotrop.b[3][2]_esd',
     * '_atom_site_anisotrop.b[3][3]_esd',
     * '_atom_site_anisotrop.u[1][1]_esd',
     * '_atom_site_anisotrop.u[1][2]_esd',
     * '_atom_site_anisotrop.u[1][3]_esd',
     * '_atom_site_anisotrop.u[2][1]_esd',
     * '_atom_site_anisotrop.u[2][2]_esd',
     * '_atom_site_anisotrop.u[2][3]_esd',
     * '_atom_site_anisotrop.u[3][1]_esd',
     * '_atom_site_anisotrop.u[3][2]_esd',
     * '_atom_site_anisotrop.u[3][3]_esd',
     * '_atom_site.pdbx_pdb_ins_code',
     * '_atom_sites.fract_transf_matrix[1][1]',
     * '_atom_sites.fract_transf_matrix[1][2]',
     * '_atom_sites.fract_transf_matrix[1][3]',
     * '_atom_sites.fract_transf_matrix[2][1]',
     * '_atom_sites.fract_transf_matrix[2][2]',
     * '_atom_sites.fract_transf_matrix[2][3]',
     * '_atom_sites.fract_transf_matrix[3][1]',
     * '_atom_sites.fract_transf_matrix[3][2]',
     * '_atom_sites.fract_transf_matrix[3][3]',
     * '_atom_sites.fract_transf_vector[1]',
     * '_atom_sites.fract_transf_vector[2]',
     * '_atom_sites.fract_transf_vector[3]',
     * '_atom_sites_footnote.id','_atom_sites_footnote.text',
     * '_audit_author.name', '_cell.z_pdb',
     * '_cell.angle_alpha', '_cell.angle_beta', '_cell.angle_gamma',
     * '_cell.length_a', '_cell.length_b', '_cell.length_c',
     * '_chem_comp.id',
     * '_chem_comp.mon_nstd_flag',
     * '_citation.id',
     * '_citation.coordinate_linkage',
     * '_citation.title', '_citation.country',
     * '_citation.journal_abbrev', '_citation.journal_volume',
     * '_citation.journal_issue', '_citation.page_first',
     * '_citation.year',
     * '_citation.journal_coden_astm',
     * '_citation.journal_coden_issn',
     * '_citation.journal_coden_csd', '_citation.journal_id_astm',
     * '_citation.journal_id_issn', '_citation.journal_id_csd',
     * '_citation.book_id_isbn',
     * '_citation.book_title',
     * '_citation.book_publisher',
     * '_citation.book_coden_isbn',
     * '_citation.details',
     * '_citation_editor.citation_id',
     * '_citation_editor.name',
     * '_citation_author.citation_id',
     * '_citation_author.name',
     * '_citation_author.ordinal',
     * '_database_2.database_code',
     * '_database_pdb_remark.id',
     * '_database_pdb_remark.text',
     * '_database_pdb_rev.date_original',
     * '_database_pdb_caveat.id',
     * '_database_pdb_matrix.origx[1][1]',
     * '_database_pdb_matrix.origx[1][2]',
     * '_database_pdb_matrix.origx[1][3]',
     * '_database_pdb_matrix.origx[2][1]',
     * '_database_pdb_matrix.origx[2][2]',
     * '_database_pdb_matrix.origx[2][3]',
     * '_database_pdb_matrix.origx[3][1]',
     * '_database_pdb_matrix.origx[3][2]',
     * '_database_pdb_matrix.origx[3][3]',
     * '_database_pdb_matrix.origx_vector[1]',
     * '_database_pdb_matrix.origx_vector[2]',
     * '_database_pdb_matrix.origx_vector[3]',
     * '_entity.id',
     * '_entity.type',
     * '_entity_poly_seq.entity_id',
     * '_entity_poly_seq.mon_id',
     * '_entity_poly_seq.num',
     * '_exptl.method',
     * '_exptl.details',
     * '_reflns.d_resolution_high',
     * '_struct.title',
     * '_struct_asym.entity_id',
     * '_struct_asym.id',
     * '_struct_biol.details',
     * '_struct_keywords.entry_id',
     * '_struct_keywords.text',
     * '_symmetry.space_group_name_h-m'/
C
C     Mark names as false in procnam
C
      do ii = 1,NUMBLOCK
      procnam(ii) = .false.
      enddo
C
C     Mark the names that cif2pdb processes as true
C
      do ii = 1, NUMusdtgs
        call hash_find(usdtgs(ii),
     *    dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *    iname)
          if (iname.gt.0) then
            procnam(iname) = .true.
          endif
      end do
      line_=2048
      nolin = precn_
      if (wide_code .eq. "yes") then
        rempref = ' '
        write(rempref(1:11),'(a6,1x,i4)')'REMARK',iremno
        write(iunpdb,'(a)') rempref(1:11)
        write(iunpdb,'(a)') rempref//'EMBEDDED CIF'
        result=prefx_(rempref,22)
        pfold_ = 119
      else
        rempref = ' '
        write(rempref(1:10),'(a6,i4)')'REMARK',iremno
        write(iunpdb,'(a)') rempref(1:10)
        write(iunpdb,'(a)') rempref(1:11)//'EMBEDDED CIF'
        result=prefx_(rempref,11)
        pfold_=68
      endif
      result=pfile_(' ')
      if(.not.result) then
        call c2perr(
     *  ' failed to open file for insertion of cif')
      endif
      result=find_(' ','head',string)
200   call cpcmnt
      result=data_(' ')
      if (result) go to 210
        call cpcmnt
        pdbline = pdbline + 2 + precn_-nolin
        numRemark = numRemark + 2 + precn_-nolin
        call close_
        return
210   string=bloc_
      saveo_=save_
      result=pdata_(string)
      if(.not.result) then
          call c2perr(
     *    ' duplicate data block '//bloc_)
      endif
      tabl_=.false.
      align_=.false.
      ixname=1
250   if(ixname.le.nname) then
        loop_ = .false.
        text_ = .false.
        if(nloop(ixname).ne.0) then
          iyname=ixname
          igood=0
          do ii = ixname,nname
            if(nloop(ixname).ne.nloop(ii)) go to 260
            if(.not.procnam(ii)) igood=igood+1
            iyname=ii
          enddo
260       continue
          if(igood.ne.0) then
C
C         We have a loop to process, force category keys back in
C
          if (ndict.gt.0) then
          do ii = ixname,iyname
            kd = ddict(ii)
            if (kd.ne.0) then
              if (aroot(kd).ne.0) kd=aroot(kd)
              if (catkey(kd)) procnam(ii) = .false.
            endif
          enddo
          endif
          do ii = ixname,iyname
            call cpcmnt
            result=find_(dname(ii),'name',string)
            if(.not.procnam(ii)) then
              pposval_ = 0
              pposnam_= 0
              result=ploop_(string(1:max(1,lastnb(string))))
            endif
          enddo
          tabl_ = .false.
          align_ = .false.
          if(igood.lt.iyname-ixname+1) then
            tabl_ = .true.
          endif
270       do ii = ixname,iyname
            if(.not.procnam(ii)) then
            call cpcmnt
            result=test_(dname(ii))
            pposnam_=posnam_
            pposval_=posval_
            pposdec_=posdec_
            pposend_=posend_
            if(posend_.gt.59) then
              pposnam_=max(0,pposnam_-21)
              pposval_=max(0,pposval_-21)
              pposdec_=max(0,pposdec_-21)
              pposend_=max(0,pposend_-21)
            endif
            if (tabl_) then
              pposnam_=0
              pposval_=0
              pposdec_=0
              pposend_=0
            endif
            if(type_.eq.'null') then
              if (long_.eq.1.and.strg_(1:1).eq.'?') then
                 result=pchar_(' ','?')
              else
                 result=pchar_(' ','.')
              endif
              type_=' '
            endif
            if(type_.eq.'char'.or.type_.eq.'numb') then
               result=dtype_(dname(ii),xdtp)
               if (result) then
                 if (xdtp.eq.'numb') then
                   dpesd=0.0D0
                   result=numd_(dname(ii),dpn,dpesd)
                   plzero_=lzero_
                   result=pnumd_(' ',dpn,dpesd)
                   type_ = ' '
                 endif
               endif
               if (type_.ne.' ') then
               result=char_(dname(ii),string)
               if(string.eq.'?'.or.string.eq.'.') then
                 tmparg=''''//
     *             string(1:1)//''''
                 result=pchar_(' ',tmparg)
               else
                 result=pchar_(' ',string(1:max(1,lastnb(string))))
               endif
               endif
               type_=' '
            endif
C           if(type_.eq.'numb') then
C             dpesd=0.0D0
C             result=numd_(dname(ii),dpn,dpesd)
C             plzero_=lzero_
C             result=pnumd_(' ',dpn,dpesd)
C           endif
            if(type_.eq.'text') then
280           result=char_(dname(ii),string)
              ll=lastnb(string)
              ls=1
              if(string(1:2).eq.'  ') ls=3
              if(ll-ls+1.le.59) then
              result=ptext_(' ',string(ls:max(ls,lastnb(string))))
              else
              result=ptext_(' ',string(ls:59))
              result=ptext_(' ',string(60:max(ls,lastnb(string))))
              endif
              if(text_) go to 280
              call tbxxeot
              type_=' '
            endif
            endif
          enddo
          if(loop_) go to 270
          endif
          ixname=iyname+1
          go to 250
        else
          call cpcmnt
          result=find_(dname(ixname),'name',temp)
          call cpcmnt
          if(.not.procnam(ixname)) then
          tabl_ = .false.
          result=test_(dname(ixname))
          pposnam_=posnam_
          pposval_=posval_
          pposdec_=posdec_
          pposend_=posend_
          if(posend_.gt.59) then
            pposnam_=max(0,pposnam_-21)
            pposval_=max(0,pposval_-21)
            pposdec_=max(0,pposdec_-21)
            pposend_=max(0,pposend_-21)
          endif
          if(type_.eq.'null') then
            if (long_.eq.1.and.strg_(1:1).eq.'?') then
               result=pchar_(temp,'?')
            else
               result=pchar_(temp,'.')
            endif
            type_=' '
          endif  
          if(type_.eq.'char'.or.type_.eq.'numb') then
            result=dtype_(temp,xdtp)
            if (result) then
              if (xdtp.eq.'numb') then
                dpesd=0.0D0
                result=numd_(temp,dpn,dpesd)
                plzero_=lzero_
                result=pnumd_(temp,dpn,dpesd)
                type_ = ' '
              endif
            endif
            if (type_.ne.' ') then
            result=char_(temp,string)
            if(string.eq.'?'.or.string.eq.'.') then
              result=pchar_(temp,''''//
     *          string(1:1)//'''')
            else
              result=pchar_(temp,string(1:max(1,lastnb(string))))
            endif
            endif
            type_=' '
          endif
C          if(type_.eq.'numb') then
C            dpesd=0.0D0
C            result=numd_(temp,dpn,dpesd)
C            plzero_=lzero_
C            result=pnumd_(temp,dpn,dpesd)
C          endif
          if(type_.eq.'text') then
290         result=char_(temp,string)
            ll=lastnb(string)
            ls=1
            if(string(1:2).eq.'  ') ls=3
            if(ll-ls+1.le.59) then
            result=ptext_(' ',string(ls:max(ls,lastnb(string))))
            else
            result=ptext_(' ',string(ls:59))
            result=ptext_(' ',string(60:max(ls,lastnb(string))))
            endif
            if(text_) go to 290
            call tbxxeot
            type_=' '
          endif
          endif
        endif
        ixname=ixname+1
        go to 250
      endif
      go to 200
C
      END
C
C
      subroutine cpcmnt
C
C     routine to copy a set of comments (if any) from the
C     input cif to the output cif
C
      include 'cif2pdb.cmn'
      external pdbreal
      external pdbint
      external pdbstr
      external typeset
      logical pcmnt_
      logical cmnt_
      character*(MAXBUF) string
      call tbxxbtab
100   result=cmnt_(string)
      if (result) then
        if(long_+posnam_.le.55) then
        pposnam_=posnam_
        result=pcmnt_(string(1:long_))
        else
        pposnam_=min(17,posnam_)
        result=pcmnt_(string(1:40))
        pposnam_=17
        result=pcmnt_(string(41:long_))
        endif
        goto 100
      endif
      call tbxxetab
      return
      end
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
       subroutine cell2mat(cell,matf2o,mato2f)
C
C      cell2mat -- convert a cell to matrices
C
C      Herbert J. Bernstein, yaya@bernstein-plus-sons.com
C      20 May 1998
C
C      Creates matrices for an orthogonal system with
C         a along x, b in the x-y plane
C
C      cell   - real*8 array of a, b, c, alpha, beta, gamma
C      matf2o - real*8 3x3 matrix to convert from fractional
C               to orthogonal
C      mato2f - real*8 3x3 matrix to convert from orthogonal
C               to fractional
C
       real*8 cell(6),matf2o(3,3),mato2f(3,3)
       real*8 cell_a, cell_b, cell_c
       real*8 cell_alpha, cell_beta, cell_gamma
       real*8 pi,torad
       data pi/3.14159 26535 9 d0/
       torad = pi/180.d0
       cell_a = cell(1)
       cell_b = cell(2)
       cell_c = cell(3)
       cell_alpha = cell(4)
       if (abs(cell(4)).lt.1.d0) 
     *   cell_alpha = acos(cell(4))*180.d0/pi
       cell_beta  = cell(5)
       if (abs(cell(5)).lt.1.d0)
     *   cell_beta  = acos(cell(5))*180.d0/pi
       cell_gamma = cell(6)
       if (abs(cell(6)).lt.1.d0)
     *   cell_gamma = acos(cell(6))*180.d0/pi
       do ii = 1,3
         do jj = 1, 3
           matf2o(ii,jj) = 0.d0
           mato2f(ii,jj) = 0.d0
         enddo
       enddo
       matf2o(1,1) = cell_a
       matf2o(1,2) = cell_b*cos(torad*cell_gamma)
       matf2o(2,2) = cell_b*sin(torad*cell_gamma)
       matf2o(1,3) = cell_c*cos(torad*cell_beta)
       matf2o(2,3) = cell_c*(cos(torad*cell_alpha)
     *                 -cos(torad*cell_beta)
     *                 *cos(torad*cell_gamma))
     *                 /sin(torad*cell_gamma)
       matf2o(3,3) = sqrt(cell_c**2-matf2o(1,3)**2
     *                -matf2o(2,3)**2)
       mato2f(1,1) = 1.d0/matf2o(1,1)
       mato2f(2,2) = 1.d0/matf2o(2,2)
       mato2f(3,3) = 1.d0/matf2o(3,3)
       mato2f(1,2) = -matf2o(1,2)
     *                 /(matf2o(1,1)*matf2o(2,2))
       mato2f(2,3) = -matf2o(2,3)
     *                 /(matf2o(2,2)*matf2o(3,3))
       mato2f(1,3) = -mato2f(1,1)*
     *                 (matf2o(1,2)*mato2f(2,3)
     *                 +matf2o(1,3)*mato2f(3,3))
       return
       end
       subroutine matmul(mat1,mat2,mat3)
C
C      multiply mat1 times mat2 and return mat3
C
C      Herbert J. Bernstein, yaya@bernstein-plus-sons.com
C      20 May 1998
C              
       real*8 mat1(3,3), mat2(3,3), mat3(3,3)
       do ii = 1,3
       do jj = 1,3
       mat3(ii,jj) = 0
       do kk = 1,3
       mat3(ii,jj) = mat3(ii,jj) + mat1(ii,kk)*mat2(kk,jj)
       enddo
       enddo
       enddo
       return
       end

       subroutine fixdec(coordstr,itarg)
       integer iper, ll, kper
       character*(*) coordstr
       character*13 blanks, ncoord
       
       blanks = '             '
       ll = len(coordstr)
       iper = index(coordstr,".")
       if (iper.eq.0) return
       kper = ll-iper
       if (kper.eq.itarg) return
       if (kper.lt.itarg) then
         if (coordstr(1:itarg-kper).ne.blanks) return
         ncoord = coordstr(itarg-kper+1:ll)//blanks(1:itarg-kper)
         coordstr = ncoord
         return
       endif
       if (coordstr(ll-(kper-itarg):ll).eq.blanks) then
         ncoord = blanks(1:kper-itarg)//coordstr
         coordstr = ncoord
       endif
       return
       end
       
       real*8 function det(mat)
C
C      compute the real*8 determinant of a real*8 matrix
C
C      Herbert J. Bernstein, yaya@bernstein-plus-sons.com
C      6 March 1998
C              
       real*8 mat(3,3)       
       det =
     *  (mat(1,1)*mat(2,2)*mat(3,3)
     *  + mat(2,1)*mat(3,2)*mat(1,3)
     *  + mat(3,1)*mat(1,2)*mat(2,3)
     *  - mat(1,3)*mat(2,2)*mat(3,1)
     *  - mat(1,2)*mat(2,1)*mat(3,3)
     *  - mat(1,1)*mat(2,3)*mat(3,2))
       return
       end
       subroutine invxfrm(mato2f,veco2f,matf2o,vecf2o)
C
C      compute the inverse of a tranform
C
C      Herbert J. Bernstein, yaya@bernstein-plus-sons.com
C      6 March 1998
C
C      mato2f - 3x3 real*8 input matrix
C      veco2f - 3 element real*8 input vector
C      matf2o - 3x3 real*8 output matrix
C      vecf2o - 3 element real*8 output vector
C
C      matf2o is the inverse of mato2f
C      vecf2o = -matf2o veco2f

       real*8 mato2f(3,3),veco2f(3),matf2o(3,3),vecf2o(3)
       real*8 newvol, det
       newvol = 1.d0/det(mato2f)
       matf2o(1,1) = newvol
     *  * (mato2f(2,2)*mato2f(3,3)-mato2f(3,2)*mato2f(2,3))
       matf2o(2,1) = newvol
     *  * (mato2f(2,3)*mato2f(3,1)-mato2f(3,3)*mato2f(2,1))
       matf2o(3,1) = newvol
     *  * (mato2f(2,1)*mato2f(3,2)-mato2f(3,1)*mato2f(2,2))

       matf2o(1,2) = newvol
     *  * (mato2f(1,3)*mato2f(3,2)-mato2f(1,2)*mato2f(3,3))
       matf2o(2,2) = newvol
     *  * (mato2f(1,1)*mato2f(3,3)-mato2f(1,3)*mato2f(3,1))
       matf2o(3,2) = newvol
     *  * (mato2f(1,2)*mato2f(3,1)-mato2f(1,1)*mato2f(3,2))

       matf2o(1,3) = newvol
     *  * (mato2f(1,2)*mato2f(2,3)-mato2f(1,3)*mato2f(2,2))
       matf2o(2,3) = newvol
     *  * (mato2f(1,3)*mato2f(2,1)-mato2f(1,1)*mato2f(2,3))
       matf2o(3,3) = newvol
     *  * (mato2f(1,1)*mato2f(2,2)-mato2f(1,2)*mato2f(2,1))
       
       do ii = 1,3
       vecf2o(ii) = 0.d0
       do jj = 1,3
       vecf2o(ii) = vecf2o(ii) - matf2o(ii,jj)*veco2f(jj)
       enddo
       enddo
       return
       end
