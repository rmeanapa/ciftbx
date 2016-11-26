#!/bin/csh
# decomp.cgi -- cgi-bin script to deliver decompressed version of
#               the file given as extra path information
#               If a cached version without the .Z extension is
#               available, and the request is made without the .Z
#               the cached version will be used.  If the call is
#               made with the .Z extension, a cached version will
#               be ignored.
#
# Herbert J. Bernstein, Bernstein + Sons
#
# 13 January 1997
# Rev 24 Feb 97
# Test for .Z extension, try for expanded file if it is there
# Rev 1 June 1997
# For sgi systems with csh truncation variable names and
# losing modifiers, change PATH_TRANSLATED to PT before modifying
# Thanks to J. Westbrook of Rutgers.
#
# Call this script as a cgi-bin script as in
#   http://URL/cgi-bin/decomp.cgi/relpath
# where relpath is the relative path of the compressed version of the
# file to deliver as uncompressed text via a compress -d call
#
# To operate corectly the /bin/echo version of echo must follow
# system V conventions sufficiently to produce an empty line
# call, below
#
#set PATH_INFO=$PATH_TRANSLATED:x
set PT=$PATH_TRANSLATED
set PATH_INFO=$PT:x
/bin/echo "Content-type: text/plain"
/bin/echo
if ( $PATH_INFO:e != "Z" ) then
  if (-e $PT) then
    cat $PT
  else
    if (-e ${PT}.Z ) then
      compress -d < $PT.Z
    else
      /bin/echo "*** decomp.cgi error: file ${PT}.Z not found ***" 
    endif
  endif
else
  if (-e $PT) then
    compress -d < $PT
  else
    /bin/echo "*** decomp.cgi error: file $PT not found ***" 
  endif
endif
