#
#  Makefile file for cif2cif 2.0.0
#  Based on Version 0.0.6 of 6 May 1997
#  Revised 2 June 1997 for web page handling
#  Revised 29 November 2009 for CIFtbx_4.1.0 release
#  H. J. Bernstein, Bernstein+Sons
#  yaya@bernstein-plus-sons.com
#
#
#
#  provide here a path to mkdecompln and rmdecompln
#
MKDECOMPLN	= ../mkdecompln
RMDECOMPLN	= ../rmdecompln
#
#  define the directory for temporary uncompressed files
#
#EXPAND		= /var/tmp
EXPAND		= .
#
#  The definition of VPATH must be the httpd virtual path to
#  the directory holding this version of the directory, i.e.
#  containing the ciftbx_n.n.n directory  MANIFEST.html
#  and README.cif2cif.html must be rebuilt if VPATH changes.
VPATH          =       NDB/mmcif/software
#VPATH          =       ~yaya/software
#VPATH           =       pb/pdb2cif
#
#  The definition of DECOMP must be the httpd virtual path to
#  the cgi-bin script decomp.cgi.  If this script is not available
#  define DECOMP as NODECOMP
DECOMP         =       /cgi-bin/yaya/decomp.cgi
#DECOMP         =       /~yaya/cgi-bin/decomp.cgi
#DECOMP          =       /pb/pdb2cif/cgi-bin/decomp.cgi
#DECOMP	=	NODECOMP
#
#
#
ZPATH           =       $(DECOMP)/$(VPATH)
#
#  provide here a path to cif_mm.dic
#
MMDICPATH		= 	../dictionaries/cif_mm.dic
#
# provide here a path to cif_expanded_jun06.dic
CIFEXJUN06PATH		=	../dictionaries/cif_expanded_jun06.dic
#
#  provide here the path to the CIFtbx2 directory
TBXPATH		=	../ciftbx.src
#
#  provide a path to the shell sh here
SHELL		=	/bin/sh
#
#  Add any necessary compilation flags on the next line
#FFLAGS		=  -g -Wall
#FFLAGS = -g -O2 -fbounds-check -Wall -Wno-unused-variable  -Wuninitialized
FFLAGS = -g -O2
#
all:	cif2cif postshar
#
#  If any of the following fail, it means the basic ciftbx
#  installation needs to be redone

$(TBXPATH)/ciftbx.o:	$(TBXPATH)/ciftbx.f \
		$(TBXPATH)/ciftbx.sys \
		$(TBXPATH)/ciftbx.cmv
		( cd $(TBXPATH) ; make ciftbx.o )
		
$(TBXPATH)/hash_funcs.o:	$(TBXPATH)/hash_funcs.f
		( cd $(TBXPATH); make hash_funcs.o )

MANIFEST.html:	MANIFEST.html.m4 Makefile
		-@rm MANIFEST.html.BAK
		-@mv MANIFEST.html MANIFEST.html.BAK
		m4 -DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < MANIFEST.html.m4 \
		> MANIFEST.html

README.cif2cif.html:	 README.cif2cif.html.m4 Makefile
		-@rm README.cif2cif.html.BAK
		-@mv README.cif2cif.html README.cif2cif.html.BAK
		m4 -DGRAPHICS=1 \
		-DZPATH=$(ZPATH) -DDECOMP=$(DECOMP) < README.cif2cif.html.m4 \
		> README.cif2cif.html

cif2cif:	cif2cif.o $(TBXPATH)/hash_funcs.o $(TBXPATH)/ciftbx.o
		$(FC) $(FFLAGS) cif2cif.o \
		$(TBXPATH)/ciftbx.o $(TBXPATH)/hash_funcs.o \
		-o cif2cif

ciftbx.sys:	$(TBXPATH)/ciftbx.sys
		ln -f -s $(TBXPATH)/ciftbx.sys ciftbx.sys

ciftbx.cmn:	$(TBXPATH)/ciftbx.cmn
		ln -f -s $(TBXPATH)/ciftbx.cmn ciftbx.cmn

ciftbx.cmf:	$(TBXPATH)/ciftbx.cmf
		ln -f -s $(TBXPATH)/ciftbx.cmf ciftbx.cmf

ciftbx.cmv:	$(TBXPATH)/ciftbx.cmv
		ln -f -s $(TBXPATH)/ciftbx.cmv ciftbx.cmv

cif2cif.o:	cif2cif.f cif2cif.cmn ciftbx.sys ciftbx.cmn ciftbx.cmf \
		ciftbx.cmv
		$(FC) $(FFLAGS) -c cif2cif.f


clean:		
		-@rm -f *.o
		-@rm -f cif2cif
		-@rm -f *.lst
		-@rm -f *.new
		-@rm -f *.diff
		-@rm -f *.BAK
		-@rm -f *.bak
		-@$(RMDECOMPLN) *.uZ


shars:	MANIFEST Makefile README.cif2cif\
		cif2cif.cmn cif2cif.f 4ins.cif 4ins.out 4ins.prt \
		4insw.out 4insuw.out \
		xtalt2.cif xtalt2.out xttne9.out xte29.out 
		-@rm -f ../Part*
		-@rm -f ../cif2cif.cshar
		-@rm -f ../cif2cif.shar
		-@rm -f ../cif2cif.cshar.Z
		-@rm -f ../cif2cif.shar.Z
		( cd .. ; makekit -c -icif2cif.src/MANIFEST \
		 -ocif2cif.src/MANIFEST -h2 -p -s 2000k )
		mv ../Part01 ../cif2cif.cshar
		( cd .. ; makekit -icif2cif.src/MANIFEST \
		 -ocif2cif.src/MANIFEST -h2 -p -s 2000k )
		mv ../Part01 ../cif2cif.shar
		-@rm -f MANIFEST.BAK
		compress ../cif2cif.cshar
		compress ../cif2cif.shar


tests:	postshar cif_mm.dic \
		xtalt2.cif 4ins.cif 4ins.out 4ins.prt cif2cif \
		4insw.out 4insuw.out cif_expanded.dic cif_expanded.duff
		-@rm -f *.lst
		-@rm -f *.new
		-@rm -f *.diff
		./cif2cif -t y < xtalt2.cif > xtalt2.new
		./cif2cif -e 9 < xtalt2.cif > xttne9.new
		./cif2cif -t y -e 29 < xtalt2.cif > xte29.new
		./cif2cif -q qtest.req -i qtest.cif > qtest.new 2> qtest.lst
		./cif2cif -t y -d cif_mm.dic -p .._ < 4ins.cif > 4ins.new \
			2> 4ins.lst
		./cif2cif -gn -u -m2048 -i 4ins.cif > 4insuw.new 2> 4insuw.lst
		./cif2cif -d cif_mm.dic -w 120 -i 4insuw.new > 4insw.new 2> 4insw.lst 
		./cif2cif  -e 99 -B read < cif_expanded.dic > cif_expanded.new
		-@diff xtalt2.out xtalt2.new > xtalt2.diff
		-@diff xttne9.out xttne9.new > xttne9.diff
		-@diff xte29.out xte29.new > xte29.diff
		-@diff 4ins.out 4ins.new > 4ins.diff
		-@diff 4ins.prt 4ins.lst > 4insprt.diff
		-@diff 4insuw.out 4insuw.new > 4insuw.diff
		-@diff 4insw.out 4insw.new > 4insw.diff
		-@diff 4insuw.prt 4insuw.lst > 4insuwprt.diff
		-@diff 4insw.prt 4insw.lst > 4inswprt.diff
		-@diff qtest.new qtest.out > qtest.diff
		-@diff qtest.prt qtest.lst > qtestprt.diff
		-@diff -b cif_expanded.dic cif_expanded.new > cif_expanded.duff.new
		-@diff cif_expanded.duff.new cif_expanded.duff > cif_expanded.diff
		cat *.diff			

cif_mm.dic:	$(MMDICPATH).Z
		ln -f -s $(MMDICPATH).Z cif_mm.dic.Z
		$(MKDECOMPLN) cif_mm.dic $(EXPAND)

cif_expanded.dic:	$(CIFEXJUN06PATH).Z
		ln -f -s $(CIFEXJUN06PATH).Z cif_expanded.dic.Z
		$(MKDECOMPLN) cif_expanded.dic $(EXPAND)

4ins.cif:	4ins.cif.Z
		$(MKDECOMPLN) 4ins.cif $(EXPAND)

4ins.out:	4ins.out.Z
		$(MKDECOMPLN) 4ins.out $(EXPAND)

4insw.out:	4insw.out.Z
		$(MKDECOMPLN) 4insw.out $(EXPAND)

4insuw.out:	4insuw.out.Z
		$(MKDECOMPLN) 4insuw.out $(EXPAND)


postshar:	
		compress 4ins.cif
		compress 4ins.out
		compress 4insuw.out
		compress 4insw.out
		touch postshar
