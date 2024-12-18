########################################################################
#                                                                      #
#             USER DEFINITIONS FOR CIF2PDB INSTALLATION                #
#                                                                      #
########################################################################
#  In the m4 document Makefile.m4 (NOT Makefile), 
#  define values for:
#
#    expand -- A directory in which to expand compressed files
#              The default is `/var/tmp'.  On some systems `/tmp'
#              may be preferable.  If expanded files are to be
#              kept long term, then holding the expanded file in
#              the current directory `.' may be preferable.
#              The directory used must allow for at least 40
#              megabytes of files for a full build and test.
#    decomp -- The httpd virtual path to the cgi-bin script 
#              `decomp'.cgi.  This is a special script which
#              decompresses files on the fly.  If this script is
#              not available on your system, define `decomp' as
#              `NODECOMP'
#    vpath  -- The httpd virtual path to the directory holding
#              this version of the directory, i.e. containing 
#              the cif2pdb_n.n.n directory.  If `decomp' is
#              defined as `NODECOMP' then `vpath' may simply
#              be set to `..'  MANIFEST.html and 
#              README.cif2pdb.html must be rebuilt if `vpath'
#              changes.
#    binpath -- The path to utility programs such as sed
#
#     
#
########################################################################
define(`expand',`/var/tmp')`#' Directory for uncompressed files, EXPAND = expand()
define(`decomp',`NODECOMP')`#' Path to `decomp'.cgi, DECOMP = decomp()
define(`vpath',`..')`#' Path to program directory, VPATH = vpath()
define(`binpath',`/usr/bin/')`#' Path to utilities such as sed, BINPATH = binpath() 
########################################################################
#                                                                      #
#                      END OF USER DEFINITIONS                         #
#                                                                      #
########################################################################
########################################################################
#                                                                      #
#            SPECIAL DEFINITIONS FOR CIF2PDB INSTALLATION              #
#                                                                      #
########################################################################
#  The following definitions may need to be modified at
#  some sites, but most often are best left alone.
#  Make changes in the m4 document Makefile.m4
#  (NOT Makefile)
#
#  There is a strong assumption that the name of the directory
#  within which this build is being made is cif2pdb_n.n.n
#  where n.n.n is the current version, and that ciftbx.src is
#  installed on the same level.
#
#  This is the path to the CIFtbx2 directory
TBXPATH         =       ../ciftbx.src
#
#
#  This is a path to cif_mm.dic
#
MMDICPATH	=	../dictionaries/cif_mm_cif2pdb.dic
#
#  provide a path to the shell sh here
SHELL           =       /bin/sh
#
#  Add any necessary compilation flags on the next line
FFLAGS          = -g
#
#  On some systems, csh is broken, and tcsh should be used.  In
#  that case, uncomment the next line
#CSHELL		=	tcsh
#
#  In can be useful to time the performance of the package
#
#  If no timer is desired, comment out the next line and uncomment
#  the one after
TIMER		= time
#TIMER		=

########################################################################
#                                                                      #
#          NO CHANGES SHOULD BE NEEDED BELOW THIS LINE                 #
#                                                                      #
########################################################################

define(`vnumb',`2.0.3')define(`vdate',`30 November 2009')#
define(concat,`$1'`$2')#
ifdef(`SITE',,`define(`SITE',`LOCAL')')`#' Makefile built for `SITE' = SITE()
#
`#'  Makefile for cif2pdb vnumb()
`#'  Version of vdate()
#
#
#
#  H. J. Bernstein, Bernstein+Sons
#  yaya@bernstein-plus-sons.com
#
#
#     make clean                    Cleans up the installation directory
#     make distclean                Restore the installation directory to
#                                   distribution state
#     make Makefile_local           Builds Makefile_local AND
#                                   Makefile from Makefile.m4
#     make (or make all)            Compiles cif2pdb
#     make tests                    Tests the compilation of cif2pdb
#
#     The following calls are used in building release kits
#
#     make README.cif2pdb.html      Builds README.cif2pdb.html 
#                                   from README.cif2pdb.html.m4
#     make MANIFEST.html            Builds MANIFEST.html from
#                                   MANIFEST.html.m4
#     make expanded                 Creates local expanded files
#     make postshar                 Compresses files after unpacking
#                                   cif2pdb.cshar or cif2pdb.shar
#     make Makefile_expanded
#     make Makefile_sdsc
#     make Makefile_ndb
#     make Makefile_pdb
#                                   Make special versions of Makefile, 
#                                   but do not replace the default 
#                                   Makefile.  The expanded version is 
#                                   similar to the default local version, 
#                                   but expands in the current directory.
#                                   in order to use one of the special
#                                   versions, make -f Makefile_...
#     make shars                    Builds cif2pdb.cshar.Z and 
#                                   cif2pdb.shar.Z in directory
#                                   ..
#
#
#
#  Provide here a path to mkdecompln and rmdecompln
#
MKDECOMPLN	=	$(CSHELL) ./mkdecompln
RMDECOMPLN	=	$(CSHELL) ./rmdecompln
#
#  Directory for shars
#
SHARDIR		=	..
#
#  define the directory for temporary uncompressed files
#
ifelse(SITE(),`NDB',,`#')EXPAND	=	/var/tmp
ifelse(SITE(),`PDB',,`#')EXPAND =	/var/tmp
ifelse(SITE(),`SDSC',,`#')EXPAND	=	/var/tmp
ifelse(SITE(),`NODECOMP',,`#')EXPAND	=	.
ifelse(SITE(),`LOCAL',,`#')EXPAND	=	expand()
#
#  The definition of VPATH must be the httpd virtual path to
#  the directory holding this version of the directory, i.e.
#  containing the cif2pdb_n.n.n directory  MANIFEST.html
#  and README.cif2pdb.html must be rebuilt if VPATH changes
ifelse(SITE(),`NDB',,`#')VPATH	=	NDB/mmcif/software
ifelse(SITE(),`PDB',,`#')VPATH	=	~yaya/software
ifelse(SITE(),`SDSC',,`#')VPATH	=	pb/pdb2cif
ifelse(SITE(),`NODECOMP',,`#')VPATH	=	..
ifelse(SITE(),`LOCAL',,`#')VPATH	=	vpath()
#
#  The definition of DECOMP must be the httpd virtual path to
`#'  the cgi-bin script `decomp'.cgi.  If this script is not available
#  define DECOMP as NODECOMP
ifelse(SITE(),`NDB',,``#'')DECOMP	=	/cgi-bin/yaya/`decomp'.cgi
ifelse(SITE(),`PDB',,``#'')DECOMP	=	/~yaya/cgi-bin/`decomp'.cgi
ifelse(SITE(),`SDSC',,``#'')DECOMP	=	/pb/pdb2cif/cgi-bin/`decomp'.cgi
ifelse(SITE(),`NODECOMP',,``#'')DECOMP	=	NODECOMP
ifelse(SITE(),`LOCAL',,``#'')DECOMP	=	decomp()
#
#
#
ifelse(SITE(),`NODECOMP',`ZPATH =	..',`ZPATH =	$(DECOMP)/$(VPATH)')
changequote(\036,\037)

#
#  The suffixes for any files to be used or built:
.SUFFIXES:	.m4 .pl .awk .oawk .ent .acif .pcif .oacif .tcif \
		.cif .pdb .tpdb .diff .wid .twid .wdiff
#
#  The locations of any programs used

CIF2PDB		=	./cif2pdb

#
#  The flags to be used in cif builds
#  And the default header file to be used
#
CIFFLAGS	=



CIFS	=	5hvp.cif 1ace.cif 2ace.cif 1crn.cif 1cro.cif 4ins.cif 1hyh.cif \
		1cwp.cif ADH041.cif BDL001.cif BDLB13.cif DDF040.cif 4hir.cif 1zrt.cif

PDBS	=	5hvp.pdb 1ace.pdb 2ace.pdb 1crn.pdb 1cro.pdb 4ins.pdb 1hyh.pdb \
		1cwp.pdb ADH041.pdb BDL001.pdb BDLB13.pdb DDF040.pdb 4hir.pdb 1zrt.pdb

WIDS	=	5hvp.wid 1ace.wid 2ace.wid 1crn.wid 1cro.wid 4ins.wid 1hyh.wid \
		1cwp.wid ADH041.wid BDL001.wid BDLB13.wid DDF040.wid 4hir.wid 1zrt.wid

TPDBS	=	5hvp.tpdb 1ace.tpdb 2ace.tpdb 1crn.tpdb 1cro.tpdb 4ins.tpdb \
		1hyh.tpdb 1cwp.tpdb ADH041.tpdb BDL001.tpdb BDLB13.tpdb \
		DDF040.tpdb 4hir.tpdb 1zrt.tpdb

TWIDS	=	5hvp.twid 1ace.twid 2ace.twid 1crn.twid 1cro.twid 4ins.twid \
		1hyh.twid 1cwp.twid ADH041.twid BDL001.twid BDLB13.twid \
		DDF040.twid 4hir.twid 1zrt.twid

CIFZS	=	5hvp.cif.Z 1ace.cif.Z 2ace.cif.Z 1crn.cif.Z 1cro.cif.Z \
		4ins.cif.Z 1hyh.cif.Z 1cwp.cif.Z ADH041.cif.Z BDL001.cif.Z \
		BDLB13.cif.Z DDF040.cif.Z 4hir.cif.Z 1zrt.cif.Z

PDBZS	=	5hvp.pdb.Z 1ace.pdb.Z 2ace.pdb.Z 1crn.pdb.Z 1cro.pdb.Z \
		4ins.pdb.Z 1hyh.pdb.Z 1cwp.pdb.Z ADH041.pdb.Z BDL001.pdb.Z \
		BDLB13.pdb.Z DDF040.pdb.Z 4hir.pdb.Z 1zrt.pdb.Z
		
WIDZS	=	5hvp.wid.Z 1ace.wid.Z 2ace.wid.Z 1crn.wid.Z 1cro.wid.Z \
		4ins.wid.Z 1hyh.wid.Z 1cwp.wid.Z ADH041.wid.Z BDL001.wid.Z \
		BDLB13.wid.Z DDF040.wid.Z 4hir.wid.Z 1zrt.wid.Z

TPDBZS	=	5hvp.tpdb.Z 1ace.tpdb.Z 2ace.tpdb.Z 1crn.tpdb.Z 1cro.tpdb.Z \
		4ins.tpdb.Z 1hyh.tpdb.Z 1cwp.tpdb.Z ADH041.tpdb.Z \
		BDL001.tpdb.Z BDLB13.tpdb.Z DDF040.tpdb.Z 4hir.tpdb.Z 1zrt.tpdb.Z

TWIDZS	=	5hvp.twid.Z 1ace.twid.Z 2ace.twid.Z 1crn.twid.Z 1cro.twid.Z \
		4ins.twid.Z 1hyh.twid.Z 1cwp.twid.Z ADH041.twid.Z \
		BDL001.twid.Z BDLB13.twid.Z DDF040.twid.Z 4hir.twid.Z 1zrt.twid.Z

DIFFS	=	5hvp.diff 1ace.diff 2ace.diff 1crn.diff 1cro.diff \
		4ins.diff 1hyh.diff 1cwp.diff ADH041.diff \
		BDL001.diff BDLB13.diff DDF040.diff 4hir.diff 1zrt.diff \
		5hvp.wdiff 1ace.wdiff 2ace.wdiff 1crn.wdiff 1cro.wdiff \
		4ins.wdiff 1hyh.wdiff 1cwp.wdiff ADH041.wdiff \
		BDL001.wdiff BDLB13.wdiff DDF040.wdiff 4hir.wdiff 1zrt.wdiff

#
FFLAGS		= -g 
#

#
#  The path to utilities such as sed
#
BINPATH		= binpath()

.cif.pdb:
	$(TIMER) $(CIF2PDB) -d cif_mm.dic -p $* < $< | \
	compress > $*.pdb.Z
	$(MKDECOMPLN) $*.pdb $(EXPAND)

.cif.wid:
	$(TIMER) $(CIF2PDB) -w yes -d cif_mm.dic -p $* < $< | \
	compress > $*.wid.Z
	$(MKDECOMPLN) $*.wid $(EXPAND)

.cif.tpdb:
	$(TIMER) $(CIF2PDB) -d cif_mm.dic -p $* < $< | \
	compress > $*.tpdb.Z
	$(MKDECOMPLN) $*.tpdb $(EXPAND)

.cif.twid:
	$(TIMER) $(CIF2PDB) -w yes -d cif_mm.dic -p $* < $< | \
	compress > $*.twid.Z
	$(MKDECOMPLN) $*.twid $(EXPAND)

.pdb.diff:
	-@/bin/rm -f $*.diff
	-@/bin/rm -f $*.spdb
	-$(BINPATH)sed -e  "1,\$$s/ 0\./  \./g" < $*.pdb | \
		$(BINPATH)sed -e "1,\$$s/ -0\./  -\./g" >$*.spdb
	-uncompress < $*.tpdb.Z | \
		$(BINPATH)sed -e "1,\$$s/ 0\./  \./g"  | \
		$(BINPATH)sed -e "1,\$$s/ -0\./  -\./g" | \
		diff - $*.spdb > $*.diff
	-/bin/rm -f $*.spdb

.wid.wdiff:
	-@/bin/rm -f $*.wdiff
	-@/bin/rm -f $*.swid
	-$(BINPATH)sed -e  "1,\$$s/ 0\./  \./g" < $*.wid | \
		$(BINPATH)sed -e "1,\$$s/ -0\./  -\./g" >$*.swid
	-uncompress < $*.twid.Z | \
		$(BINPATH)sed -e "1,\$$s/ 0\./  \./g"  | \
		$(BINPATH)sed -e "1,\$$s/ -0\./  -\./g" | \
		diff - $*.swid > $*.wdiff
	-/bin/rm -f $*.swid

all:		cif2pdb postshar

allpdbs:	$(PDBS)

allwids:	$(WIDS)

Makefiles:	Makefile_sdsc Makefile_ndb Makefile_pdb Makefile_expanded

Makefile_local:	Makefile.m4
		m4 Makefile.m4 > Makefile_local
		cp Makefile_local Makefile

Makefile_sdsc:	Makefile.m4
		m4 -DSITE=SDSC Makefile.m4 > Makefile_sdsc

Makefile_ndb:	Makefile.m4
		m4 -DSITE=NDB Makefile.m4 > Makefile_ndb

Makefile_pdb:	Makefile.m4
		m4 -DSITE=PDB Makefile.m4 > Makefile_pdb

Makefile_expanded:	Makefile.m4
		m4 -DSITE=NODECOMP Makefile.m4 > Makefile_expanded


#
#  If any of the following fail, it means the basic ciftbx
#  installation needs to be redone

$(TBXPATH)/ciftbx.o:	$(TBXPATH)/ciftbx.f \
		$(TBXPATH)/ciftbx.sys \
		$(TBXPATH)/ciftbx.cmv
		( cd $(TBXPATH) ; make ciftbx.o )

$(TBXPATH)/hash_funcs.o:	$(TBXPATH)/hash_funcs.f
		( cd $(TBXPATH); make hash_funcs.o )

ciftbx.sys:	$(TBXPATH)/ciftbx.sys
		ln -f -s $(TBXPATH)/ciftbx.sys ciftbx.sys

ciftbx.cmn:	$(TBXPATH)/ciftbx.cmn
		ln -f -s $(TBXPATH)/ciftbx.cmn ciftbx.cmn

ciftbx.cmf:	$(TBXPATH)/ciftbx.cmf
		ln -f -s $(TBXPATH)/ciftbx.cmf ciftbx.cmf

ciftbx.cmv:	$(TBXPATH)/ciftbx.cmv
		ln -f -s $(TBXPATH)/ciftbx.cmv ciftbx.cmv


cif2pdb.o:	cif2pdb.f cif2pdb.cmn ciftbx.sys ciftbx.cmn ciftbx.cmf \
		ciftbx.cmv
		$(FC) $(FFLAGS) -c cif2pdb.f

cif2pdb:	cif2pdb.o $(TBXPATH)/hash_funcs.o $(TBXPATH)/ciftbx.o cif2pdb.o
		$(FC) $(FFLAGS) cif2pdb.o $(TBXPATH)/ciftbx.o \
		$(TBXPATH)/hash_funcs.o -o cif2pdb


postshar:
		touch $(TPDBS)
		touch $(CIFS)
		touch $(TWIDS)
		compress $(TPDBS)
		compress $(TWIDS)
		compress $(CIFS)
		ln -f -s $(MMDICPATH).Z cif_mm.dic.Z
		touch postshar


clean:		
		-@rm *.o
		-@rm *.diff
		-@rm *.wdiff
		-@rm *.BAK
		-@rm *.bak
		-@rm cif2pdb
		-@rm *.pdb.Z
		-@rm *.wid.Z
		-@$(RMDECOMPLN) *.uZ
		-@$(RMDECOMPLN) .DECOMP/*.uZ

distclean:	clean
		-@uncompress $(TPDBS)
		-@uncompress $(CIFS)
		-@uncompress $(TWIDS)
		-@rm cif_mm.dic.Z
		-@rm postshar
		-@rm ciftbx.*

tests:		postshar testcif2pdb
		-@cat *diff

testcif2pdb:	cif2pdb cif_mm.dic $(CIFS) $(PDBS) $(WIDS) \
		$(DIFFS) $(TPDBZS) $(TWIDZS)


$(CIFS):	$(CIFZS)
		$(MKDECOMPLN) $@ $(EXPAND)
$(TPDBS):	$(TPDBZS)
		$(MKDECOMPLN) $@ $(EXPAND)
$(TWIDS):	$(TWIDZS)
		$(MKDECOMPLN) $@ $(EXPAND)

$(TPDBZS):	$(CIFZS)
		make `echo $@ | $(BINPATH)sed -e "s/\.tpdb\.Z//"`.cif
	  	$(CIF2PDB) -d cif_mm.dic -p \
		  `echo $@ | $(BINPATH)sed -e "s/\.tpdb\.Z//"` < \
		  `echo $@ | $(BINPATH)sed -e "s/\.tpdb\.Z//"`.cif | \
		  compress > `echo $@ | $(BINPATH)sed -e "s/\.tpdb\.Z//"`.tpdb.Z

$(TWIDZS):	$(CIFZS)
		make `echo $@ | $(BINPATH)sed -e "s/\.twid\.Z//"`.cif
	  	$(CIF2PDB) -w yes -d cif_mm.dic -p \
		  `echo $@ | $(BINPATH)sed -e "s/\.twid\.Z//"` < \
		  `echo $@ | $(BINPATH)sed -e "s/\.twid\.Z//"`.cif | \
		  compress > `echo $@ | $(BINPATH)sed -e "s/\.twid\.Z//"`.twid.Z

cif_mm.dic:	cif_mm.dic.Z
		$(MKDECOMPLN) cif_mm.dic $(EXPAND)

MANIFEST.html:	MANIFEST.html.m4 Makefile
		-@rm MANIFEST.html.BAK
		-@mv MANIFEST.html MANIFEST.html.BAK
		m4 -DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < MANIFEST.html.m4 \
		> MANIFEST.html

README.cif2pdb.html:	README.cif2pdb.html.m4 Makefile
		-@rm README.cif2pdb.html.BAK
		-@mv README.cif2pdb.html README.cif2pdb.html.BAK
		m4 -DGRAPHICS=1 \
		-DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < \
		README.cif2pdb.html.m4 \
		> README.cif2pdb.html

expanded:	cif2pdb.cshar.Z cif2pdb.shar.Z \
		cif_mm.dic.Z $(CIFS) $(TPDBS) .DECOMP
		$(MKDECOMPLN) cif2pdb.cshar $(EXPAND)
		$(MKDECOMPLN) cif2pdb.shar $(EXPAND)
		$(MKDECOMPLN) cif_mm.dic $(EXPAND)

$(DIFFS):	$(PDBS) $(TPDBZS) $(WIDS) $(TWIDZS)
$(PDBS):	cif2pdb $(CIFZS)



shars:		$(CIFS) $(TPDBS) $(TWIDS) DISCUSS.cif2pdb.html \
		IUCR_POLICY IUCR_POLICY.html Makefile Makefile.m4 \
		NOTICE NOTICE.html \
		README README.cif2pdb.html README.cif2pdb.html.m4 \
		cif2pdb.cmn cif2pdb.f mkdecompln rmdecompln \
		MANIFEST MANIFEST.html MANIFEST.html.m4 $(SHARDIR)
		-@/bin/rm -f $(SHARDIR)/cif2pdb.cshar.Z
		-@/bin/rm -f $(SHARDIR)/cif2pdb.shar.Z
		(cd .. ; makekit -c -icif2pdb.src/MANIFEST -ocif2pdb.src/MANIFEST -h2 -p -s30000k)
		mv ../Part01 $(SHARDIR)/cif2pdb.cshar
		compress $(SHARDIR)/cif2pdb.cshar
		-@/bin/rm -f cif2pdb.cshar.Z
		ln -s $(SHARDIR)/cif2pdb.cshar.Z cif2pdb.cshar.Z
		(cd .. ; makekit -icif2pdb.src/MANIFEST -ocif2pdb.src/MANIFEST -h2 -p -s30000k)
		mv ../Part01 $(SHARDIR)/cif2pdb.shar
		compress $(SHARDIR)/cif2pdb.shar
		-@/bin/rm -f cif2pdb.shar.Z
		ln -s $(SHARDIR)/cif2pdb.shar.Z cif2pdb.shar.Z

.DECOMP:
		-@mkdir .DECOMP
$(SHARDIR):
		-@mkdir $(SHARDIR)
