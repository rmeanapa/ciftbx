<!doctype html public "-//IETF//DTD HTML 2.0//EN">
<HTML>
define(vnumb,`4.1.1')
define(concat,`$1'`$2')
ifelse(DECOMP(),`NODECOMP',
`define(`XXPATH',)',
`define(`XXPATH',`ZPATH()/concat(ciftbx_,vnumb())/cif2cif.src/')')
<HEAD>
<TITLE>
README.cif2cif
</TITLE> 
</HEAD> 
<BODY BACKGROUND="../.logos/cif2cif_wallpaper.jpg">
<a href="http://www.iucr.ac.uk/iucr-top/welcome.html">
<img alt="[IUCr Home Page]" src="../.logos/iucrhome.jpg"></a>
<a href="http://www.iucr.ac.uk/iucr-top/cif/home.html">
<img alt="[CIF Home Page]" src="../.logos/cifhome.jpg"></a>
<A HREF="../README.html"> <IMG SRC="../.logos/ciftbxHomeButton.jpg" ALT="[ciftbx home]"></A>
<A HREF="../cyclops.src/README.cyclops.html"> <IMG SRC="../.logos/cyclopsButton.jpg" 
ALT="[cyclops.src]"></A> 
<A HREF="../ciftbx.src/README.ciftbx.html"> <IMG SRC="../.logos/ciftbxButton.jpg" 
ALT="[ciftbx.src]"></A> 
<hr>
<CENTER>
<IMG SRC="../.logos/cif2cif.jpg" ALT = " ">
</CENTER>
<H2 ALIGN=CENTER>README.cif2cif</H2>
<H3 ALIGN=CENTER>
Information for cif2cif 2.0.1, 22 January 2024<BR>
</H3>
<HR>
<H3 ALIGN=CENTER>
Before using this software, please read the <BR>
<A HREF = "../NOTICE.html"> <IMG SRC="../.logos/noticeButton.jpg" ALT="NOTICE"></A><BR>
 and please read the IUCr<BR>
<A HREF="../IUCR_POLICY.html"> <IMG SRC="../.logos/policyButton.jpg" ALT="Policy"></A><BR>
 on the Use of the Crystallographic Information File (CIF)</A>
</H3> 
<HR>
<PRE>
    \ | /      \ | /
     \|/        \|/ 
    -- --&gt;&gt;&gt;&gt;&gt;&gt;-- --               c i f 2 c i f  ...... CIF COPY PROGRAM
     /|\        /|\
    / | \      / | \                          Version 2.0.1
                                            22 January 2024


 cif2cif  is a fortran program using CIFtbx2 to copy a CIF on standard
 -------  input to an equivalent CIF on standard output, while checking
          data names against dictionaries and reformating numbers with
          esd's to conform to the rule of 19.  A request list may be 
          specified

           cif2cif
                  by

                  Copyright &#169; 1997, 1998, 2000, 2005, 2009, 2024
                  Herbert J. Bernstein (yayahjb@gmail.com)
                  Bernstein + Sons
                  11 Riverside Drive Apt 5UE
                  New York, NY 10023-1317, U.S.A.

           based on suggestions by

                  Sydney R. Hall (syd@crystal.uwa.edu.au)
                  Crystallography Centre
                  University of Western Australia
                  Nedlands 6009, AUSTRALIA

           with the request list handling suggested by the program
           QUASAR by Sydney R. Hall and

                  Rolf Sievers (unc411@dbnrhrz1.bitnet)
                  Institut fuer Anorganische-Chemisches der Universitaet
                  Gerhard-Domagk-Str.
                  Bonn, GERMANY

</PRE>
<P>
<P>
In order to ensure continuing availability of source code and
documentation cif2cif and its documentation are subject to copyright.
This does not prevent you from using
the program, from making copies and changes,
but prevents the creation of &quot;closed source&quot; versions
out of the open source versions.
See <a href="NOTICE.html">NOTICE</a>.

<P>
Science is best served
when the tools we use are fully understood by those who wield those
tools and by those who make used of results obtained with those tools.
When a scientific tool exists as software, access to source code is
an important element in achieving full understanding of that tool.
As our field evolves and new versions of software are required,
access to source allows us to adapt our tools quickly and effectively.
<P>
In the early days of software development,
most scientific software source code was freely and openly shared with a
minimum of formalities.  These days, it appears that carefully drawn
legal documents are necessary to protect free access to the source
code of scientific software.  We are all deeply indebted to
<a href="http://www.stallman.org">Richard Stallman</a> for showing
us how a creative combination of copyrights and seemingly restrictive
licenses could give us truly unfettered freedom to use programs,
to read their source code and to develop new versions.  The
<a href="http://www.gnu.org">GNU project</a>, and the
<a href="http://www.linux.org">Linux project</a> have shown that
an open source approach works.  We use the GNU General Public
License
(the <a href="http://www.gnu.org/copyleft/gpl.html">&quot;GPL&quot;</a>)
for our program starting with the release of CIFtbx3.  Older versions use the license from <a
href="http://www.OpenRasMol.org">OpenRasmol</a>.
The OpenRasMol conditions
for use have correctly been called &quot;GPL-like&quot;.
<P>
If you are a user of this program, you will find that
the copyrights and notices ask little more of you than that you
avoid mistakes by others by keeping the notices with copies,
display scientific integrity by citing your sources properly and
treating this like other shared scientific developments by not
inferring a warranty.  If you are a software developer and wish to
incorporate what you
find here into new code, or to pick up bits and pieces and used them
in another context, the situation becomes more complex.  Read the
copyrights and notices carefully.  You will find that they are
&quot;infectious&quot;.  <b>Whatever you make from our Open Source code
must itself be offered as Open Source code.</b>  In addition, in
order to allow users to understand what has changed and to ensure
orderly development you have to describe your changes.
<P>
<hr>

<P>
 cif2cif reads the input CIF from the standard input device (normally
 device 5). An optional STAR data name dictionary (in DDL format) is opened.
 A reformatted copy of the input CIF is written to standard output (device 6).
 Messages are output to the standard error device (normally device 0).
 Note that the PARAMETER 'MAXBUF' should contain the maximum number of char-
 acters contained on a single text line. The default value is 2048.
<P>
 In a unix-like environment, the program is run as:
<PRE>
			  
      cif2cif [-i input_cif] [-o output_cif] [-d dictionary] [-a aliaso_] \
         [-c catck] [-e esdlim_] [-f command_file] [-g guess_type] \
         [-h htmlt] [-m maxline] [-n numasst] [-p prefix] \
         [-q request_list] [-t tabl_] [-u unfold] [-w wrap] [-B read|cif2read] \
         [input_cif [output_cif [dictionary [request_list]]]]
      where:
               input_cif defaults to $CIF2CIF_INPUT_CIF or stdin
               output_cif defaults to $CIF2CIF_OUTPUT_CIF or stdout
               dictionary defaults to $CIF2CIF_CHECK_DICTIONARY
               multiple dictionaries may be specified 
               request_list defaults to $CIF2CIF_REQUEST_LIST
               input_cif of "-" is stdin, output_cif of "-" is stdout
               request_list of "-" is stdin
               -a has values of t or 1 or y vs. f or 0 or n,
               -c has values of t or 1 or y vs. f or 0 or n,
               -e has integer values (e.g. 9, 19(default) or 29),
               -g has values of t, 1 or y vs. f, 0 or n (default t),
               -h has values of t or 1 or y vs. f or 0 or n, default f
               -m has values from 80 to 2048 for a maximum line width
               -n has values of t or 1 or y vs. f or 0 or n, default f
               -p has a string value in which "_" is replaced by blank,
               -t has values of t, 1 or y vs. f, 0 or n (default f),
               -u has values of t, 1 or y vs. f, 0 or n (default f),
               -w has values of 0 for no wrap or a column to wrap on
</PRE>
<P>
 The program can reformat CIFs with long lines to CIFs with shorter line
 by use of the '-w wrap' option, and can reformat CIFs that have been
 wrapped to a short line length to a longer line length with the '-u y'
 option.  Numbers with standard uncertainties can be reformatted to
 different rules (e.g. from rule of 9 to rule of 29), provided the
 type of the data name is known to be numeric or can be guessed from
 context.  The guessing of types for data names not specified in a
 dictionary can be suppressed by use of the '-g n' option.
<P>
 The extraction of selected fields using a request list is supported
 by the '-q request_list' option.  If a data name is requested
 that does not appear in the input CIF, a dummy entry for the
 requested name will be placed in the output CIF.
<P>
<H2>
1.  INSTALLATION
</H2>
Here is the recommended procedure for implementing and testing this
version of cif2cif.
<P>
1.0.  Before you try to install this version of cif2cif
<PRE>   
     *** ========================================================== ***
     *** ========================================================== ***
     *** ==&gt;&gt;&gt; You must have ciftbx version 4.0.1 or greater  &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; installed in a directory named ciftbx.src.     &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; The scripts mkdecompln and rmdecompln, which   &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; come with ciftbx, must be installed in the     &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; top level directory and executable.            &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; To test cif2cif, you must have a compressed    &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; copy of the dictionary cif_mm.dic in a         &lt;&lt;&lt;== ***
     *** ==&gt;&gt;&gt; directory named dictionaries.                  &lt;&lt;&lt;== ***
     *** ========================================================== ***
     *** ========================================================== ***
</PRE>
    The directory structure within which you will work is
<PRE>
                  top level directory
                  -------------------
                           |
                           |
            ------------------------------
            |              |             | 
       dictionaries   ciftbx.src     cif2cif.src
       ------------   ----------     -----------
</PRE>
            

<P>
   You may have acquired this package in one of several forms.
   The most likely are as a &quot;C-shell Archive,&quot; a &quot;Shell Archive&quot;,
   or as separate files.  The idea is to get to separate files,
   all in the same directory, named cif2cif.src, parallel to
   the directory ciftbx.src, but let's start with the possibility
   that you got the package as one big file, i.e. in one of
   the archive file formats.  Place the archive in the top level
   directory.
<PRE>
   *** ========================================================== ***
   *** ========================================================== ***
   *** ==&gt;&gt;&gt; The files in this kit will unpack into a       &lt;&lt;&lt;== ***
   *** ==&gt;&gt;&gt; directory named cif2cif.src.  It is a good idea&lt;&lt;&lt;== ***
   *** ==&gt;&gt;&gt; to save the current contents of cif2cif.src    &lt;&lt;&lt;== ***
   *** ==&gt;&gt;&gt; and then to make the directory empty           &lt;&lt;&lt;== ***
   *** ========================================================== ***
   *** ========================================================== ***
</PRE>
   If you are on a machine which does not provide a unix-like shell,
   you will need to take apart the archive by hand using a text
   editor.  We'll get to that in a moment.
<P>
   1.1.  ON A UNIX MACHINE
<P>
   If you have the shell archive on a unix machine, follow the
   instructions at the front of the archive, i.e. save the
   uncompressed archive file as &quot;file&quot;, then, if the archive 
   is a &quot;Shell Archive&quot; execute &quot;sh file&quot;.  If the archive is 
   a &quot;C-Shell Archive&quot; execute &quot;csh file&quot;.
<P>
   1.2. IF YOU DON'T HAVE UNIX
<P>
   If sh or csh are not available, then it is best to start with the
   &quot;C-Shell Archive&quot; and do the steps that follow.  If you must
   use the &quot;Shell Archive&quot; you should be aware that the lines
   you want to extract have been prefixed with &quot;X&quot;, while most of
   the lines you want to discard have not.  For a &quot;C-Shell Archive&quot;
   such prefixes are rare and the file is easier to read.  Assume
   you have a &quot;C-Shell Archive&quot;.
<P>
   Use your editor to separate the different parts of the file into
   individual files in your workspace. Each part starts with
   a lot of unixisms, then several blank lines and then two lines
   which identify the file, and most importantly, contain the
   text &quot;CUT_HERE_CUT_HERE_CUT_HERE&quot;  You can look at the line
   before and the line after to see if you are at the head or
   tail of a file.  Use your editor to search for the &quot;CUT_HERE&quot;
   lines. Each part is carefully labeled and indicates the  
   recommended filename for the separated file. On some machines
   these filenames may need to be altered to suit the OS or compiler.
<P>
   1.3.  MANIFEST
<P>
   The partitions are as follows:
<PRE>
   part  filename                   description

     1  <A HREF="../COPYING">COPYING</A>                          GPL (GNU General Public License)
     2  <A HREF="../NOTICE">NOTICE</A>                           Notices
     3  <A HREF="README.cif2cif">cif2cif.src/README.cif2cif</A>       additional information on cif2cif
     4  <A HREF="MANIFEST">cif2cif.src/MANIFEST</A>             a list of files in the kit
     5  <A HREF="Makefile">cif2cif.src/Makefile</A>             a preliminary control file for make
     6  <A HREF="XXPATH()4ins.cif">cif2cif.src/4ins.cif</A>             example mmcif file used to test cif2cif
     7  <A HREF="XXPATH()4ins.out">cif2cif.src/4ins.out</A>             example mmcif output from test of cif2cif
     8  <A HREF="4ins.prt">cif2cif.src/4ins.prt</A>             example mmcif list file from test of cif2cif
     9  <A HREF="XXPATH()4insuw.out">cif2cif.src/4insuw.out</A>           example mmcif unwrapped long output
    10  <A HREF="4insuw.prt">cif2cif.src/4insuw.prt</A>           example mmcif list unwrapped long output
    11  <A HREF="XXPATH()4insw.out">cif2cif.src/4insw.out</A>            example mmcif wrapped output
    12  <A HREF="4insw.prt">cif2cif.src/4insw.prt</A>            example mmcif list wrapped output
    13  <A HREF="cif2cif.cmn">cif2cif.src/cif2cif.cmn</A>          cif2cif common block
    14  <A HREF="cif2cif.f">cif2cif.src/cif2cif.f</A>            cif2cif fortran source
    15  <A HREF="qtest.cif">cif2cif.src/qtest.cif</A>            quasar mode test input cif
    16  <A HREF="qtest.out">cif2cif.src/qtest.out</A>            quasar mode test output cif
    17  <A HREF="qtest.req">cif2cif.src/qtest.req</A>            quasar mode test request file
    18  <A HREF="qtest.prt">cif2cif.src/qtest.prt</A>            example cif file used to test cif2cif
    19  <A HREF="xtalt2.cif">cif2cif.src/xtalt2.cif</A>           example cif file used to test cif2cif
    20  <A HREF="xtalt2.out">cif2cif.src/xtalt2.out</A>           example cif output from test of cif2cif
    21  <A HREF="xte29.out">cif2cif.src/xte29.out</A>            example cif output from test of cif2cif
    22  <A HREF="xttne9.out">cif2cif.src/xttne9.out</A>           example cif output from test of cif2cif
    23  <A HREF="cif_expanded.duff">cif2cif.src/cif_expanded.duff</A>    changes to cif_expanded.dic
</PRE>
<H3>     
2. COMPILING AND EXECUTING
</H3>
   Here are the recommended steps for a UNIX system. Vary this 
   according to the requirements of your OS and compiler.  You will 
   simplest to work if you place the cif2cif files together in a 
   common subdirectory named 'cif2cif.src'.   Be very careful if you 
   place them in a directory with other files, since some of the build 
   and test instructions remove or overwrite existing files, especially 
   with extensions such as '.o', '.lst',  or '.diff'.  On a UNIX system,
   most of what you need to do to build and test cif2cif is laid out in 
   Makefile.  *** Be sure to examine and edit this file appropriately 
   before using it.***  But, once properly edited, all you should need 
   to do is 'make clean' to remove old object files, 'make all' to 
   build new version of 'cif2cif' and 'make tests' to test what you
   have built.
<P>   
   For non-UNIX-like environments, you will have to provide 
   replacements for iargc, getarg and getenv.  The following are
   reasonable possibilities:

<PRE>
         integer function iargc(dummy)
         iargc=0
         return
         end

         subroutine getarg(narg,string)
         integer narg
         character*(*) string
         string=char(0)
         return
         end

         subroutine getenv(evar,string)
         character*(*) evar,string
         string=char(0)
         if(evar.eq.&quot;CIF2CIF_INPUT_CIF&quot;)
        *  string='INPCIF.CIF'//char(0)
         if(evar.eq.&quot;CIF2CIF_OUTPUT_CIF&quot;)
        *  string='OUTCIF.CIF'//char(0)
         if(evar.eq.&quot;CIF2CIF_CHECK_DICTIONARY&quot;) 
        *  string='CIF_CORE.DIC'//char(0)
         return
         end
</PRE>
 This combination of substitute routines would &quot;wire-in&quot; cif2cif to
 read its input cif from a file named INPCIF.CIF, write its output
 cif to a file named OUTCIF.CIF, and check names against CIF_CORE.DIC
<H3>
FILES USED
</H3>
<PRE>
       dictionary input         input   on device 2
       Reformatted CIF          output  on device 6 ('stdout')
       Input CIF                input   on device 2, if a file, 5  if 'stdin'
       Message device           output  on device 0 ('stderr')
       Direct access            in/out  on device 3
</PRE>
<H2>
TEST files
</H2>
 Three test CIFs are provided.  xtal2.cif is a test file borrowed from
 xtal_gx (file xtest2.cif at <a href="http://www.iucr.org/iucr-top/cif/software/xtal/gx/">
 www.iucr.org/iucr-top/cif/software/xtal/gx/</a>,
 provided by S. R. Hall.   4ins.cif is an mmCIF file created from the 
 PDB entry 4INS by G.G. Dodson, E. J. Dodson, D. C. Hodgkin, 
 N.W. Isaacs and M. Vijayan (1989) by the program pdb2cif 
 (P.E. Bourne, F.C. Bernstein and H.J. Bernstein, 1996,  see 
 <a href="http://www.bernstein-plus-sons.com/software/pdb2cif">
 www.bernstein-plus-sons.com/software/pdb2cif</a>).  qtest.cif 
 is the test cif from the quasar release by Hall and Sievers (1994).  
 A modified version of the request list, qtest.req, is also included.
 In addition, the June 2006 version of the expanded DDLm CIF core
 dictionary from the dictionaries directory is used as a test file.
<P>
 xtalt2.cif provides good test cases for the conversion of esd's.
 The command
<PRE>
    cif2cif -t y &lt; xtalt2.cif &gt; xtalt2.new
</PRE>
 ensures that all esd's follow the rule of 19, while
<PRE>
    cif2cif -t y -e 29 &lt; xtalt2.cif &gt; xte29.new
</PRE>
 converts esd's to the rule of 29.  The difference between the
â two rules is that for the rule of 19, all esd's lie between 2 and
 19, so that an esd of (1) has to be converted to (10), while
 for the rule of 29, all esd's lie between 3 and 29, so that
 an esd of (2) also has to be converted, in this case to (20).
 The option &quot;-t y&quot; tidies the output to tab stops.
<P>
 One last test with this file is to use the command
<PRE>
   cif2cif -e 9 &lt; xtalt2.cif &gt; xttne9.new
</PRE>
 to copy the original cif spacing and to use the rule of 9 on esd's
<P>
 4ins.cif has many comments, text fields and dense loops.  The
 test in the Makefile tests handling of these items and adds the
 additional complication of processing a prefix &quot;.._&quot; with the
 command
<PRE>
   cif2cif -t y -p .._ &lt; 4ins.cif &gt; 4ins.new
</PRE>
 The output spacing is controlled by the program.
<P>
 Two additional tests are provided.  The command
<PRE>
   cif2cif -gn -u -m2048 -i 4ins.cif &gt; 4insuw.new
</PRE>
 unwraps the 80 character version of 4ins.cif to the maximum
 line length, not using a dictionary and suppressing the guessing
 of numeric data types.  This is followed by a pass wrapping
 the resulting cif to 120 columns with the command
<PRE> 
   cif2cif -d cif_mm.dic -w 120 -i 4insuw.new > 4insw.new
</PRE>

<P>
 The quasar mode may be tested by the command
<PRE>
  cif2cif -i qtest.cif -o qtest.new -q qtest.req
</PRE>

<H3>
CHANGES
</H3>
<p>
 Version 2.0.1 (22 January 2024) adds support for optional
 html output.
<p>
 Version 2.0.0 (29 November 2009) uses CIFtbx 4 support for
 the DDLm bracketed constructs to allow DDLm files to be
 copied when -B read is specified.  -B cif2read is provided
 as a preliminary attempt to support CIF 2 files.
<P>
 Version 1.0.1 (7 April 2005) added the feature of including a
 dummy entry in the output cif for each name that appears on a
 request list,l but not on the input cif.  The new command
 line option -gn adds the ability to suppress guessing of data 
 types for names not defined in a dictionary.  Formatting
 of unwrapped output cifs was improved.  Work supported in 
 part by the International Union of Crystallography.
<P>
 Version 1.0.0 (2 April 2005) added support for line folding and
 unfolding using -w (for wrap), -u (for unfold) and -m (for
 maximum line length).  Work supported in part by the
 International Union of Crystallography.
<P>
 Version 0.1.0 (8 July 2000) fixed a parse error in the command line
 and an error in the comments.
<P>
 Version 0.0.9 (30 May 2000) changed the variable named tab to xxtab
 to avoid a conflict.
<P>
 Version 0.0.8 (2 April 1998) added code to preserve the distinction
 between esds which are numerically identical but presented with varying
 numbers of trailing zeros.   To enable this feature, the command line
 argument -e must be given a negative value large enough to span the
 desired range, e.g. "-e -19" to allow both .6870(10) and .687(1) to be
 handled.  The variables esddig_ and pesddig_ introduced with ciftbx 2.6
 are used.
 A bug in the preserving leading zeros of unlooped numeric values
 was also fixed.
<P>
 Version 0.0.7 (7 August 1997) added support for quasar-style request
 lists.  In addition, a special form of request for "data_which_contains:"
 will find the data block which contains at least one of the following
 tags.
 <P>
 Version 0.0.6 (12 May 1997) added support of negative values of
 esdlim_, indicating a range of esd's from 1 to -esdlim_.  The practical 
 use of this change is to use the command line parameter -e-9999 to
 copy most esds unchanged.  The new ciftbx variables decp_, pdecp_,
 lzero_ and plzero_ are used internally to copy decimal points and
 leading zeros more faithfully.  The 4ins example has been updated to
 use cif_mm.dic (0.9.01).  The command-line option "-c catck" was added.
 The default is not to check categories.
<P>
 Version 0.0.5 (2 December 96) adds support for copying of global_
 sections and corrects a typo in the error message issued when
 the output CIF cannot be opened.
<P> 
 Version 0.0.4 (24 September 96) adapts to the 'ciftbx.cmn' and
 'ciftbx.sys' reorganization in ciftbx 2.5.1.
<P>
 Version 0.0.3 (24 July 96) used ciftbx 2.5.0 to preserve tabs,
 fixed a case in which string quoted in the input cif may
 not have been quoted in the output cif, preserved more white
 space, and copied comments within loop headers.  The lateral
 position of &quot;loop_&quot; is preserved.  A bug in copying final comments
 which caused the a comment to be duplicated or some final comments
 to be lost was fixed.
<P>
 Version 0.0.2 (24 June 96) changed the default for -t to n 
 instead of y, and used ciftbx 2.4.6 to increase the faithfulness 
 of the copy.

<H3>
KNOWN PROBLEMS
</H3>
cif2cif does not copy white space exactly, and will reformat
some data values.  Always compare the original to the output.
<P>
<HR>
Updated 14 August 2024<P>
<address>
<A HREF="mailto:yaya@bernstein-plus-sons.com">yaya@bernstein-plus-sons.com</A>
</BODY>
</HTML>
