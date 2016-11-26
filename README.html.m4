<!doctype html public "-//IETF//DTD HTML 2.0//EN">
<HTML>
define(vnumb,`4.1.0')
define(concat,`$1'`$2')
ifelse(DECOMP(),`NODECOMP',
`define(`XXPATH',)',
`define(`XXPATH',`ZPATH()/concat(ciftbx_,vnumb())/')')
<HEAD>
<TITLE>
CIFtbx, CYCLOPS, cif2cif
</TITLE> 
</HEAD> 
ifdef(`GRAPHICS',
<BODY BACKGROUND=".logos/ciftbx_wallpaper.jpg">
,
<BODY>
)
ifdef(`GRAPHICS',
<a href="http://www.iucr.org/iucr-top/welcome.html">
<img alt="[IUCr Home Page]" src=".logos/iucrhome.jpg"></a>
<a href="http://www.iucr.org/iucr-top/cif/home.html">
<img alt="[CIF Home Page]" src=".logos/cifhome.jpg"></a>
<A HREF="ciftbx.src/README.ciftbx.html"><IMG SRC=".logos/ciftbxButton.jpg" 
ALT="[ciftbx.src]"></A>
<A HREF="cyclops.src/README.cyclops.html"><IMG SRC=".logos/cyclopsButton.jpg" 
ALT="[cyclops.src]"></A> 
<A HREF="cif2cif.src/README.cif2cif.html"><IMG SRC=".logos/cif2cifButton.jpg" 
ALT="[cif2cif.src]"></A> 
<A HREF="cif2xml.src/README.cif2xml.html"><IMG SRC=".logos/cif2xmlButton.jpg" 
ALT="[cif2xml.src]"></A> 
<A HREF="cif2pdb.src/README.cif2pdb.html"><IMG SRC=".logos/cif2pdbButton.jpg" 
ALT="[cif2pdb.src]"></A> 
<hr>
,)
<H2 ALIGN=CENTER>
CIFtbx, Cyclops, cif2cif, cif2xml, cif2pdb README FIRST Information
</H2>
<center>
<b>CIFtbx Version 4.1.0 Release, November 2009</b><br />
<font size=-1>Copyright &#169; 2009
Sydney R. Hall and Herbert J. Bernstein</font><br />
<b>Cyclops Version 2.1.5, November 2009</b><br />
<font size=-1>Copyright &#169; 2009 Sydney R. Hall and Herbert J.
Bernstein</font><br />
<b>cif2cif Version 2.0.0, November 2009</b><br />
<font size=-1>Copyright &#169; 1997, 1998, 2000, 2005, 2009 Herbert J. Bernstein</font><br />
<b>cif2xml Version 0.1.1 beta, November 2009</b><br />
<font size=-1>Copyright &#169; 2000, 2009 Herbert J. Bernstein</font><br />
<b>cif2pdb Version 2.0.3, November 2009</b><br />
<font size=-1>Copyright &#169; 2000, 2005, 2009 Herbert J. Bernstein and
Frances C. Bernstein</font>
</center>
<P>
<CENTER>
Sydney R. Hall<BR>
Crystallography Centre,
University of Western Australia, Nedlands 6009, Australia<BR>
Phone:  +61 9 380 2725, e-mail: 
<script language="javascript" type="text/javascript">
<!-- 
      var name = "syd@";
      var domain = "crystal.uwa";
      var domext = ".edu.au";
      document.write ("<a href=\"mailto:" + name + domain + domext + "\">" + name + domain + domext+"</a>"); 
// -->
</script>
<noscript>
syd&#64;crystal.uwa&#46;edu.au
</noscript>
<BR>
</CENTER>
<P>
<CENTER>
Herbert J. Bernstein<BR>
Bernstein+Sons, 5 Brewster Lane, Bellport, NY 11713-2803, USA<BR>
Phone: 1-631-286-1339, e-mail: 
<script language="javascript" type="text/javascript">
<!-- 
      var name = "yaya@";
      var domain = "bernstein-plus-sons";
      var domext = ".com";
      document.write ("<a href=\"mailto:" + name + domain + domext + "\">" + name + domain + domext+"</a>"); 
// -->
</script>
<noscript>
yaya&#64;bernstein-plus-sons&#46;com
</noscript>
<BR>
</CENTER>
<P>
<CENTER>
Frances C. Bernstein<BR>
Bernstein+Sons, 5 Brewster, Bellport, NY 11713-2803, USA<BR>
Phone: 1-631-286-1339, e-mail: 
<script language="javascript" type="text/javascript">
<!-- 
      var name = "fcb@";
      var domain = "bernstein-plus-sons";
      var domext = ".com";
      document.write ("<a href=\"mailto:" + name + domain + domext + "\">" + name + domain + domext+"</a>"); 
// -->
</script>
<noscript>
fcb&#64;bernstein-plus-sons&#46;com
</noscript>
<BR>
</CENTER>
<P>
<hr>
<b>YOU MAY REDISTRIBUTE THE CBFLIB PACKAGE UNDER THE TERMS OF THE <a href=gpl.txt>GPL</a>.
<P>
 ALTERNATIVELY YOU MAY REDISTRIBUTE THE CBFLIB API UNDER THE TERMS
 OF THE <a href=lgpl.txt>LGPL<a/>.</b>
</CENTER>
<HR>

Work on these programs has been supported by grants from various government
agencies and the IUCr.  See the documentation of the individual programs
for further information.
<B>
<P>
    CIFtbx4, Cyclops, cif2cif, cif2xml and cif2pdb are free software; 
    you can redistribute this software and/or modify
    this software under the terms of the GNU General Public License 
    as published by the Free Software Foundation; either version 2 of 
    the License, or
    (at your option) any later version.
<P>
    Alternatively you may reditribute and/or modify the CIFtbx API (but not the programs)
    under the terms of the GNU Lesser General Public 
    License as published by the Free Software Foundation; either 
    version 2.1 of the License, or (at your option) any later version.
<P>
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
<P>
    You should have received a copy of the GNU General Public License
    along with this software; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
<P>
    You should have received a copy of the GNU Lesser General Public License
    along with this software; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
<P>
<CENTER>
Before using this software, please read the <BR>
<A HREF = "NOTICE.html"> <IMG SRC=".logos/noticeButton.jpg" ALT="NOTICE"></A><BR>
 and please read the IUCr<BR>
<A HREF="IUCR_POLICY.html"> <IMG SRC=".logos/policyButton.jpg"
ALT="Policy"></A><BR>
 on the Use of the Crystallographic Information File (CIF)
</CENTER>
</B> 
<HR>
<h4 align="center">Compiled binaries for some architectures (ix86 Linux, 
Sun Sparc, SGI)
are available from the <br><a href="http://arcib.dowling.edu">Accessible 
Retired 
Computers in Biology (ARCiB) project at http://arcib.dowling.edu</a></h4>
<hr>
<P>
CIFtbx is a tool box of Fortran routines for manipulating CIF data.
CYCLOPS2, cif2cif, cif2xml and cif2pdb are open source Fortran programs based on
CIFtbx.  CYCLOPS2 checks STAR data names against data name dictionaries.
cif2cif copies a CIF to a CIF while checking data names against
dictionaries and reformating numbers with esd's to conform to the rule
of 19.  A request list may be specified.  cif2xml is a variant of
cif2cif which produces XML output.  cif2pdb converts mmCIF files
to a PDB or WPDB format.
<P>
In order to ensure continuing availability of source code and
documentation most programs and documents on this site
are subject to copyright.  This does not prevent you from using
the programs, from making copies and changes,
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
for our programs starting with the release of CIFtbx3.  Older versions use the license from 
<a
href="http://www.OpenRasMol.org">OpenRasmol</a>.
The OpenRasMol conditions
for use have correctly been called &quot;GPL-like&quot;.
<P>
If you are a user of these programs, you will find that
the copyrights and notices ask little more of you than that you
avoid mistakes by others by keeping the notices with copies,
display scientific integrity by citing your sources properly and
treating this like other shared scientific developments by not
inferring a warranty.  If you are a software developer and wish to
incorporate what you
find here into new code, or to pick up bits and pieces and use them
in another context, the situation becomes more complex.  Read the
notices carefully.  You will find that they are
&quot;infectious&quot;.  <b>Whatever you make from our Open Source code
must itself be offered as Open Source code.</b>  In addition, in
order to allow users to understand what has changed and to ensure
orderly development you have to describe your changes. 
<P><HR>
<A HREF="CIFtbx_Primer.pdf">CIFtbx Primer Plus</A> (320K in PDF format)
provides documentation on CIFtbx.
<HR>
<P>
In order to get a complete release of CIFtbx4, CYCLOPS2, cif2cif, cif2xml 
you will need either a C-shell archive or a shell archive of each
package.  That is, you need either
<P>
<A HREF="ciftbx.cshar.Z" > ciftbx.cshar.Z,</A>
<A HREF="cyclops.cshar.Z" > cyclops.cshar.Z,</A>
<A HREF="cif2cif.cshar.Z" > cif2cif.cshar.Z,</A>
<A HREF="cif2xml.cshar.Z" > cif2xml.cshar.Z </A> and
<A HREF="cif2pdb.cshar.Z" > cif2pdb.cshar.Z </A> or
<P>
<A HREF="ciftbx.shar.Z" > ciftbx.shar.Z,</A>
<A HREF="cyclops.shar.Z" > cyclops.shar.Z,</A>
<A HREF="cif2cif.shar.Z" > cif2cif.shar.Z,</A> 
<A HREF="cif2xml.shar.Z" > cif2xml.shar.Z </A> and
<A HREF="cif2pdb.shar.Z" > cif2pdb.shar.Z </A>
<P>
You do not need both sets.  If you cannot decompress files at your end, 
then you may wish to have
<A HREF="XXPATH()ciftbx.cshar">
ciftbx.cshar,</A>
<A HREF="XXPATH()cyclops.cshar">
cyclops.cshar,</A>
<A HREF="XXPATH()cif2cif.cshar">
cif2cif.cshar,</A> 
<A HREF="XXPATH()cif2xml.cshar">
cif2xml.cshar</A> and
<A HREF="XXPATH()cif2pdb.cshar">
cif2pdb.cshar </A>

<P>In order to test the programs, you will
need at least the following CIF dictionaries:
<A HREF="dictionaries/cif_core.dic.Z"> cif_core.dic.Z </A>
<A HREF="dictionaries/cif_mm.dic.Z"> cif_mm.dic.Z </A> and
<A HREF="dictionaries/cif_mm_pdbx.dic.Z"> cif_mm_pdbx.dic.Z </A>
Check the IUCr web page at 
<A HREF="http://www.iucr.org">
http://www.iucr.org </A> for latest versions.
<P>    
The directory structure within which you will work is
<PRE> 
                                 top level directory
                                 -------------------
                                         |
                                         |
     --------------------------------------------------------------------------
     |              |             |              |              |             |
<A 
HREF="dictionaries">dictionaries</A>   <A 
HREF="ciftbx.src">ciftbx.src</A>     <A 
           HREF="cyclops.src">cyclops.src</A>   <A 
                         HREF="cif2cif.src">cif2cif.src</A>   <A 
                                                     HREF="cif2xml.src">cif2xml.src</A>   <A
                                                                    HREF="cif2pdb.src">cif2pdb.src</A>
------------   ----------     -----------   -----------   -----------   -----------
</PRE>
<P>      
In order to save space and transmission time, it is possible
that each archive will be compressed before you pick it up.
Place the archives in the top level directory and the dictionaries
into the 'dictionaries' directory.  Be sure to uncompress the
archives before trying to unpack them.  If you are on a UNIX
system, leave the dictionaries compressed.
<P>
Each archive contains a detailed README file for the package
involved.  Please do read the README.  The general approach is:
<P>
     1.  unpack the archive
<P>
     2.  make all
<P>
     3.  make tests
<P>     
Since all this code depends on CIFtbx, you must install CIFtbx
before installing CYCLOPS2, cif2cif, cif2xml or cif2pdb.  If you don't need
CYCLOPS2 or cif2cif or cif2xml or cif2pdb, you may leave them out, but they are good
examples of working CIFtbx applications, so we suggest taking a
look at them if you are new to CIFtbx.
<HR>
<address>
For further information contact Syd Hall (
<script language="javascript" type="text/javascript">
<!-- 
      var name = "syd@";
      var domain = "crystal.uwa";
      var domext = ".edu.au";
      document.write ("<a href=\"mailto:" + name + domain + domext + "\">" + name + domain + domext+"</a>"); 
// -->
</script>
<noscript>
syd&#64;crystal.uwa&#46;edu.au
</noscript>
)
or Herbert Bernstein 
(
<script language="javascript" type="text/javascript">
<!-- 
      var name = "yaya@";
      var domain = "bernstein-plus-sons";
      var domext = ".com";
      document.write ("<a href=\"mailto:" + name + domain + domext + "\">" + name + domain + domext+"</a>"); 
// -->
</script>
<noscript>
yaya&#64;bernstein-plus-sons&#46;com
</noscript>
)

</address>
<HR>
<FONT FACE="Arial,Helvetica,Times" Size="2">
<center>
The CIFtbx site is provided courtesy of <a
href="http://www.bernstein-plus-sons.com">Bernstein + Sons,</a>
Information
Systems Consultants.
</center>
Portions &#169 Copyright Herbert J. Bernstein 2000, 2005, 2009.  All Rights
Reserved.  Other
copyrights apply.  License to copy most pages on this site will be
granted if certain open source and other conditions are adhered to.
See <a
`define(`XXPATH',`ZPATH()/concat(ciftbx_,vnumb())/')')

href="http://www.bernstein-plus-sons.com/software/concat(ciftbx_,vnumb())/NOTICE.html">
http://www.bernstein-plus-sons.com/software/concat(ciftbx_,vnumb())/NOTICE.html</a>.
To request permission to copy any elements from our pages, or to send
comments about our pages, please write to <a
href="mailto:webmaster@bernstein-plus-sons.com">webmaster@bernstein-plus-sons.co
m</a>.
If you establish a link to one of our pages,
please let us know!  <i>Updated 30 November 2009</i>
</FONT>

</BODY>
</HTML>
