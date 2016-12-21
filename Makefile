#
#  Makefile for ciftbx, cyclops, cif2cif, cif2xml web related files
#  H. J. Bernstein, Bernstein+Sons, 17 July 2000
#  The valid make calls are:
#
#     make README.html      rebuilds README.html from README.html.m4
#     make MANIFEST.html    rebuilds MANIFEST.html from MANIFEST.html.m4
#     make expanded         creates local expanded files
#     make all              makes all .src kits
#     make alltests         test all .src kits
#
#===============  The following definitions may need to be modified for
#===============  your system
#

#  The locations of any programs used
#  The most common alternative are given immediately below
M4	=	m4
MKDECOMPLN	=	./mkdecompln
RMDECOMPLN	=	./rmdecompln
EXPAND		=	/var/tmp
#
#  The definition of VPATH must be the httpd virtual path to
#  the directory holding this version of the directory, i.e.
#  containing the ciftbx_n.n.n directory  MANIFEST.html 
#  and README.html must be rebuilt if VPATH changes.
VPATH		=	NDB/mmcif/software
#VPATH		=	~yaya/software
#VPATH		=	pb/pdb2cif
#
#  The definition of DECOMP must be the httpd virtual path to
#  the cgi-bin script decomp.cgi.  If this script is not available
#  define DECOMP as NODECOMP
DECOMP		=	/cgi-bin/yaya/decomp.cgi
#DECOMP		=	/~yaya/cgi-bin/decomp.cgi
#DECOMP		=	/pb/pdb2cif/cgi-bin/decomp.cgi
#DECOMP	=	NODECOMP
#
#
#
ZPATH		=	$(DECOMP)/$(VPATH)
#
MANIFEST.html:	MANIFEST.html.m4 Makefile
		-@rm MANIFEST.html.BAK
		-@mv MANIFEST.html MANIFEST.html.BAK
		m4 -DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < MANIFEST.html.m4 \
		> MANIFEST.html

README.html:	README.html.m4 Makefile
		-@rm README.html.BAK
		-@mv README.html README.html.BAK
		m4 -DGRAPHICS=1 \
		-DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < README.html.m4 \
		> README.html

ciftbx.src:	expanded
		sh ciftbx.shar

cyclops.src:	expanded
		sh cyclops.shar

cif2cif.src:	expanded
		sh cif2cif.shar

cif2pdb.src:	expanded
		sh cif2pdb.shar

cif2xml.src:	expanded
		sh cif2xml.shar
		
all:		ciftbx.src cyclops.src cif2cif.src cif2pdb.src cif2xml.src
		(cd ciftbx.src; make all)
		(cd cyclops.src; make all)
		(cd cif2cif.src; make all)
		(cd cif2pdb.src; make all)
		(cd cif2xml.src; make all)

alltests:	all
		(cd ciftbx.src; make tests)
		(cd cyclops.src; make tests)
		(cd cif2cif.src; make tests)
		(cd cif2pdb.src; make tests)
		(cd cif2xml.src; make tests)

clean:
		-@(cd ciftbx.src; make clean)
		-@(cd cyclops.src; make clean)
		-@(cd cif2cif.src; make clean)
		-@(cd cif2pdb.src; make clean)
		-@(cd cif2xml.src; make clean)
		-@rm expanded
		-@rm *.uZ
		-@rm *shar
		

expanded:	cif2cif.shar.Z \
		dictionaries/cif_core.dic.Z dictionaries/cif_mm.dic.Z \
		ciftbx.shar.Z \
		cyclops.shar.Z \
		cif2pdb.shar.Z \
		cif2xml.shar.Z
		$(MKDECOMPLN) cif2cif.shar .
		(cd dictionaries; ../$(MKDECOMPLN) cif_core.dic . )
		(cd dictionaries; ../$(MKDECOMPLN) cif_mm.dic . )
		$(MKDECOMPLN) ciftbx.shar .
		$(MKDECOMPLN) cyclops.shar .
		$(MKDECOMPLN) cif2xml.shar .
		$(MKDECOMPLN) cif2pdb.shar .
		touch expanded
