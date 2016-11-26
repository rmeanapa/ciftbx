#  Makefile file for cyclops 2.1.4
#  Version of 1 November 2009
#  H. J. Bernstein, Bernstein+Sons
#  yaya@bernstein-plus-sons.com
#
#  revised from 2.1.3, 6 May 1997 to use mkdecompln
#  and to use current dictionaries
#
#  provide here a path to the CIFtbx2 directory
#
TBXPATH		=	../ciftbx.src
#
#  provide here a path to mkdecompln and rmdecompln
#
MKDECOMPLN	= ../mkdecompln
RMDECOMPLN	= ../rmdecompln
#
#  provide here a path to cif_mm.dic and cif_core.dic
#
MMDICPATH		= ../dictionaries/cif_mm.dic
COREDICPATH		= ../dictionaries/cif_core.dic
COREAUG08PATH	= ../dictionaries/cif_expanded_aug08.dic
COREOCT09PATH	= ../dictionaries/cif_expanded_oct09.dic
#
#  define SHELL as needed
#
SHELL		= /bin/sh
#
# provide flags for fortran compilations here
FFLAGS		=  -g -O0
#
all:		cyclops

cyclops:	cyclops.o \
		$(TBXPATH)/hash_funcs.o \
		$(TBXPATH)/ciftbx.o
		$(FC) $(FFLAGS) cyclops.o \
			$(TBXPATH)/ciftbx.o $(TBXPATH)/hash_funcs.o -o cyclops

cyclops.o:	cyclops.f cyclops.cmn ciftbx.sys \
		ciftbx.cmf ciftbx.cmv
		$(FC) $(FFLAGS) -c cyclops.f

$(TBXPATH)/ciftbx.o:	$(TBXPATH)/ciftbx.f \
		$(TBXPATH)/ciftbx.sys \
		$(TBXPATH)/ciftbx.cmv $(TBXPATH)/clearfp.f
		( cd $(TBXPATH) ; make ciftbx.o )

$(TBXPATH)/hash_funcs.o:	$(TBXPATH)/hash_funcs.f
		( cd $(TBXPATH) ; make hash_funcs.o )

ciftbx.sys:	$(TBXPATH)/ciftbx.sys
		ln -s -f $(TBXPATH)/ciftbx.sys ciftbx.sys

ciftbx.cmf:	$(TBXPATH)/ciftbx.cmf
		ln -s -f $(TBXPATH)/ciftbx.cmf ciftbx.cmf

ciftbx.cmv:	$(TBXPATH)/ciftbx.cmv
		ln -s -f $(TBXPATH)/ciftbx.cmv ciftbx.cmv

mtest.prt:	$(TBXPATH)/mtest.prt
		ln -s -f $(TBXPATH)/mtest.prt mtest.prt

test.prt:	$(TBXPATH)/test.prt
		ln -s -f $(TBXPATH)/test.prt test.prt

clean:		
		-@rm -f *.o
		-@rm -f *.lst
		-@rm -f cyclops
		-@rm -f STARDICT
		-@rm -f STARCHEK
		-@rm -f STARCMD
		-@rm -f *.new
		-@rm -f *.diff
		-@rm -f *.BAK
		-@rm -f *.bak
		-@$(RMDECOMPLN) *.uZ

shars:		MANIFEST Makefile README.cyclops \
	 	cyclops.f cyclops.cmn cyclops_test.prt mtest.cyc \
		../NOTICE ../COPYING
		-@rm -f ../Part*
		-@rm -f ../cyclops.cshar
		-@rm -f ../cyclops.shar
		-@rm -f ../cyclops.cshar.Z
		-@rm -f ../cyclops.shar.Z
		( cd .. ; makekit -c -icyclops.src/MANIFEST \
		  -ocyclops.src/MANIFEST -h2 -p -s 1000k )
		mv ../Part01 ../cyclops.cshar
		( cd .. ; makekit -icyclops.src/MANIFEST \
		  -ocyclops.src/MANIFEST -h2 -p -s 1000k )
		mv ../Part01 ../cyclops.shar
		-@rm -f MANIFEST.BAK
		compress ../cyclops.cshar
		compress ../cyclops.shar

		
tests:	cyclops mtest.prt cif_mm.dic cif_core.dic  mtest.cyc cyclops_test.prt
		-@rm -f cyclops_test.lst
		-@rm -f mtest_cyc.diff
		-@rm -f cyclops_prt_lst.diff
		-@rm -f STARDICT
		-@rm -f STARCMD
		-@rm -f cyclops_test2.lst
		-@rm -f STARCHEK
		-@rm -f STARCHK2
		-@rm -f cyc_cyc.diff
		-@rm -f cyclops_lst_lst.diff
		echo "#DICT" > STARDICT
		echo "cif_core.dic" >> STARDICT
		echo "cif_mm.dic" >> STARDICT
		./cyclops mtest.prt STARCHEK 2> cyclops_test.lst
		echo "-o STARCHK2 -d cif_core.dic -d cif_mm.dic" > STARCMD
		./cyclops mtest.prt -f STARCMD 2> cyclops_test2.lst
		-@diff cyclops_test.prt cyclops_test.lst > cyclops_prt_lst.diff
		-@diff mtest.cyc STARCHEK > mtest_cyc.diff
		-@diff cyclops_test.lst cyclops_test2.lst > cyclops_lst_lst.diff
		-@diff STARCHEK STARCHK2 > cyc_cyc.diff
		cat cyclops_prt_lst.diff
		cat mtest_cyc.diff			
		cat cyclops_lst_lst.diff
		cat cyc_cyc.diff

ddlmtests:  cyclops test.prt cif_core_aug08.dic cif_core_oct09.dic
		echo "-o STARCHK3 -v t -d cif_core_aug08.dic" > STARCMD
		./cyclops test.prt -f STARCMD 2> cyclops_ddlmtest3.lst
		echo "-o STARCHK4 -v t -d cif_core_oct09.dic" > STARCMD
		./cyclops test.prt -f STARCMD 2> cyclops_ddlmtest4.lst
		echo "-o STARCHK5 -v t -d cif_core_aug08.dic" > STARCMD
		./cyclops mtest.prt -f STARCMD 2> cyclops_ddlmtest5.lst
		echo "-o STARCHK6 -v t -d cif_core_oct09.dic" > STARCMD
		./cyclops mtest.prt -f STARCMD 2> cyclops_ddlmtest6.lst

cif_mm.dic:	$(MMDICPATH).Z
		ln -f -s $(MMDICPATH).Z cif_mm.dic.Z
		$(MKDECOMPLN) cif_mm.dic

cif_core.dic:	$(COREDICPATH).Z
		ln -f -s $(COREDICPATH).Z cif_core.dic.Z
		$(MKDECOMPLN) cif_core.dic
		
cif_core_aug08.dic: $(COREAUG08PATH).Z
		ln -f -s $(COREAUG08PATH).Z cif_core_aug08.dic.Z
		$(MKDECOMPLN) cif_core_aug08.dic

cif_core_oct09.dic: $(COREOCT09PATH).Z
		ln -f -s $(COREOCT09PATH).Z cif_core_oct09.dic.Z
		$(MKDECOMPLN) cif_core_oct09.dic
