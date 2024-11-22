!module hash_functions
!implicit none
!private
!public :: hash_init, hash_find, hash_fnext, hash_snext, hash_lower
!
!       hash_funcs.f -- a library of hash table management routines
!
!                                      by
!
!                              Herbert J. Bernstein
!                                Bernstein + Sons
!                    5 Brewster Lane, Bellport, NY 11713-0177, USA
!                       email: yaya@bernstein-plus-sons.com
!
!       work on these routines done in part at Brookhaven National
!       Laboratory, under contract to the U.S. Department of Energy
!
!       work on these routines as been done in part at Dowling College 
!       under contract to the International Union of Crystallography 
!       and under grants from the National Science Foundation and 
!       the U.S. Department of Energy.
!
!       Copyright (C) Herbert J. Bernstein 1997 -- 2006
!
!       hash_funcs.f90 is free software; you can redistribute this software 
!       and/or modify this software under the terms of the GNU General 
!       Public License as published by the Free Software Foundation; 
!       either version 2 of the License, or (at your option) any later version.
!
!       Alternatively you may reditribute and/or modify hash_funcs.f90
!       under the terms of the GNU Lesser General Public 
!       License as published by the Free Software Foundation; either 
!       version 2.1 of the License, or (at your option) any later version.
!
!       This software is distributed in the hope that it will be useful,
!       but WITHOUT ANY WARRANTY; without even the implied warranty of
!       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!       GNU General Public License for more details.
!
!       You should have received a copy of the GNU General Public License
!       along with this software; if not, write to the Free Software
!       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!       You should have received a copy of the GNU Lesser General Public License
!       along with this software; if not, write to the Free Software
!       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!-------------------------------------------------------------------------------
!
!       Routines
!
!       hash_init          Initializes a hash table controlled list
!                          call hash_init( data_structure_args )
!
!       hash_find          Searches for a string in a list
!                          call hash_find( name, data_structure_args, ifind )
!
!       hash_fnext         Searches for the next matching string in a list
!                          call hash_fnext( name, data_structure_args, ifind, icurr )
!
!       hash_store         Inserts a new string in a list
!                          call hash_store( name, data_structure_args, ifind )
!
!       hash_snext         Inserts a string in a list, allowing duplicates
!                          call hash_store( name, data_structure_args, ifind, icurr )
!
!       hash_value         Integer function returns index into hash_list
!                          ih = hash_value( name, hash_length )
!
!       The necessary data_structure_args for these routines are
!
!          name_list   -- an array of character strings
!                         character(len=*) name_list(list_length)
!          chain_list  -- chain pointers for searches
!                         integer chain_list(list_length)
!          list_length -- the size of the list arrays
!                         integer list_length
!          num_list    -- number of entries in the list
!                         integer num_list
!          hash_table  -- the initial hashed pointers
!                         integer hash_table
!          hash_length -- the size of the hash table
!                         integer hash_length
!
!       The three remaining arguments are
!
!          name        -- string to search for
!                         character(len=*) name
!          ifind       -- return value, 0 for not found (hash_find)
!                         or list full (hash_store), otherwise
!                         the index in name_list of the entry
!          icurr       -- the prior matching index
!
!       The relationship among the arrays used is:
!
!       hash_table is an array (preferably of a modest prime
!       dimension) which starts containing all zeros, which are
!       replaced by pointers to entries in name_list, based
!       values returned by hash_value ranging from 1 to hash_length.
!       Each name is placed in name_list.  A initial zero is placed
!       in the matching entry in chain_list, when the first entry
!       is made.  When a new entry with the same hash_value must be
!       placed a pointer is inserted into chain_list to hook the
!       values together.

! Initialization routine for a hash table controlled list
!     name_list   -- a list of character strings
!     chain_list  -- chain pointers for searches
!     list_length -- the size of the list arrays
!     num_list    -- number of entries in the list
!     hash_table  -- the initial hashed pointers
!     hash_length -- the size of the hash table
subroutine hash_init( name_list, chain_list, list_length, num_list, hash_table, hash_length )
    implicit none
    integer,          intent(in)    :: list_length, hash_length
    character(len=*), intent(inout) :: name_list(list_length)
    integer,          intent(inout) :: chain_list(list_length)
    integer,          intent(inout) :: hash_table(hash_length)
    integer   :: idummy, num_list, i
    character :: cdummy
    num_list = 0
    do i = 1, hash_length
        hash_table(i) = 0
    enddo
    cdummy = name_list(1)
    idummy = chain_list(1)
    idummy = list_length
    idummy = num_list
    return
end subroutine hash_init 

! Search routine for a hash table controlled list
!    name        -- string to find
!    name_list   -- a list of character strings
!    chain_list  -- chain pointers for searches
!    list_length -- the size of the list arrays
!    num_list    -- number of entries in the list
!    hash_table  -- the initial hashed pointers
!    hash_length -- the size of the hash table
!    ifind       -- returned index or 0
subroutine hash_find( name, name_list, chain_list, list_length, num_list, hash_table, hash_length, ifind )
    implicit none
    integer,          intent(in)  :: list_length, hash_length
    character(len=*), intent(in)  :: name, name_list(list_length)
    integer,          intent(in)  :: chain_list(list_length), num_list
    integer,          intent(in)  :: hash_table(hash_length)
    integer,          intent(out) :: ifind
    character(len=200) :: lcname
    integer            :: lenn, hash_value, ih, ip, lastnb, idummy
    ifind = 0
    lenn  = lastnb(name)
    call hash_lower(name(1:lenn),lcname(1:lenn))
    ih    = hash_value(lcname(1:lenn),hash_length)
    ip    = hash_table(ih)
    !do while
    !if( ip .eq. 0 ) return
    100 if( ip .eq. 0 ) return
    if( name_list(ip) .eq. lcname(1:lenn) )then
        ifind = ip
        return
    else
        ip = chain_list(ip)
        go to 100
    endif
    !enddo
    idummy = num_list
end subroutine hash_find
           
! Search routine for the next matching string in a list
!    name        -- string to find
!    name_list   -- a list of character strings
!    chain_list  -- chain pointers for searches
!    list_length -- the size of the list arrays
!    num_list    -- number of entries in the list
!    hash_table  -- the initial hashed pointers
!    hash_length -- the size of the hash table
!    ifind       -- returned index or 0
!    icurr       -- current match or 0
subroutine hash_fnext( name, name_list, chain_list, list_length, num_list, hash_table, hash_length, ifind, icurr )
    implicit none
    integer,          intent(in)  :: list_length, hash_length, num_list, icurr
    character(len=*), intent(in)  :: name_list(list_length)
    integer,          intent(in)  :: chain_list(list_length), hash_table(hash_length)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: ifind
    character(len=200) :: lcname
    integer            :: hash_value, ih, ip, lastnb, idummy, lenn
    lenn  = lastnb(name)
    call hash_lower(name(1:lenn),lcname(1:lenn))
    ifind = 0
    if( icurr .eq. 0 )then
        ih = hash_value(lcname(1:lenn),hash_length)
        ip = hash_table(ih)
    else
        ip = chain_list(icurr)
    endif
    ! do while
    !if( ip .eq. 0 ) return
    100 if( ip .eq. 0 ) return
    if( name_list(ip) .eq. lcname(1:lenn) )then
        ifind = ip
        return
    else
        ip = chain_list(ip)
        go to 100
    endif
    !enddo
    idummy = num_list
end subroutine hash_fnext

! Store routine for a hash table controlled list
!    name        -- string to find
!    name_list   -- a list of character strings
!    chain_list  -- chain pointers for searches
!    list_length -- the size of the list arrays
!    num_list    -- number of entries in list
!    hash_table  -- the initial hashed pointers
!    hash_length -- the size of the hash table
!    ifind       -- index of entry or 0 (table full)
subroutine hash_store( name, name_list, chain_list, list_length, num_list, hash_table, hash_length, ifind )
    implicit none
    integer,          intent(in)    :: hash_length, list_length
    character(len=*), intent(in)    :: name
    character(len=*), intent(inout) :: name_list(list_length)
    integer,          intent(inout) :: chain_list(list_length), hash_table(hash_length)
    integer,          intent(inout) :: num_list
    integer,          intent(out)   :: ifind
    character(len=200) :: lcname
    integer :: lenn, hash_value, ih, ip, iq, lastnb
    lenn  = lastnb(name)
    call hash_lower(name(1:lenn),lcname(1:lenn))
    ifind = 0
    ih    = hash_value(lcname(1:lenn),hash_length)
    ip    = hash_table(ih)
    iq    = 0
    100 if( ip .eq. 0 ) go to 200
    if( name_list(ip) .eq. lcname(1:lenn) )then
        ifind=ip
        return
    else
        iq = ip
        ip = chain_list(ip)
        go to 100
    endif
    200 if( num_list .lt. list_length )then
        num_list             = num_list+1
        name_list(num_list)  = lcname(1:lenn)
        chain_list(num_list) = 0
        if( iq .eq. 0 )then
            hash_table(ih) = num_list
        else
            chain_list(iq) = num_list
        endif
        ifind = num_list
        return
    else
        ifind = 0
        return
    endif
end subroutine hash_store

! Store routine for a hash table controlled list
!    name        -- string to find
!    name_list   -- a list of character strings
!    chain_list  -- chain pointers for searches
!    list_length -- the size of the list arrays
!    num_list    -- number of entries in list
!    hash_table  -- the initial hashed pointers
!    hash_length -- the size of the hash table
!    ifind       -- index of entry or 0 (table full)
!    icurr       -- current match or 0
subroutine hash_snext( name, name_list, chain_list, list_length, num_list, hash_table, hash_length, ifind, icurr )
    implicit none
    integer,          intent(in)    :: list_length, hash_length, icurr
    character(len=*), intent(in)    :: name
    integer,          intent(inout) :: num_list
    character(len=*), intent(inout) :: name_list(list_length)
    integer,          intent(inout) :: hash_table(hash_length), chain_list(list_length)
    integer,          intent(out)   :: ifind
    character(len=200) :: lcname
    integer :: lenn, hash_value, ih, ip, iq, lastnb
    ifind = 0
    ih    = 0
    lenn  = lastnb(name)
    call hash_lower(name(1:lenn),lcname(1:lenn))
    if( icurr .eq. 0 )then
        ih = hash_value(lcname,hash_length)
        ip = hash_table(ih)
        iq = 0
    else
        ip = chain_list(icurr)
        iq = icurr
    endif
    100 if( ip .eq. 0 ) go to 200
    if( name_list(ip) .eq. lcname(1:lenn) )then
        ifind = ip
        return
    else
        iq = ip
        ip = chain_list(ip)
        go to 100
    endif
    200 if( num_list .lt. list_length )then
        num_list             = num_list+1
        name_list(num_list)  = lcname(1:lenn)
        chain_list(num_list) = 0
        if( iq .eq. 0 )then
            hash_table(ih) = num_list
        else
            chain_list(iq) = num_list
        endif
        ifind = num_list
        return
    else
        ifind = 0
        return
    endif
end subroutine hash_snext

! Function to return a hash value of string name to fit
! a hash table of length hash_length
integer function hash_value( name, hash_length ) 
    implicit none
    character(len=*), intent(in) :: name
    integer,          intent(in) :: hash_length
    character(len=200) :: lcname
    integer :: lastnb
    integer :: id, ii, i, ic, lenn
    lenn       = lastnb(name)
    call hash_lower(name(1:lenn),lcname(1:lenn))
    hash_value = 1
    id         = 0
    do ii = 1, lenn
        i  = 1+lenn-ii
        ic = ichar(lcname(i:i))
        if( ic .ge. 32 )then
            hash_value = mod(hash_value*(ic-32),hash_length)+1
            id         = id+1
            if( id .gt. 3 ) return
        endif
    enddo
    return
end function hash_value
       
! Function to return a lower case version of string in lcstr
subroutine hash_lower( string, lcstr )
    implicit none
    character(len=*), intent(in)  :: string
    character(len=*), intent(out) :: lcstr
    character(len=1)  :: c
    integer           :: lastnb, lstring, llcstr, ipos, i, index
    character(len=26) :: UC, lc
    UC      = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    lc      = "abcdefghijklmnopqrstuvwxyz"
    lstring = lastnb(string)
    llcstr  = len(lcstr)
    if(  llcstr < lstring ) lstring                 = llcstr
    if( lstring < llcstr  ) lcstr(lstring+1:llcstr) = ' '
    do ipos = 1, lstring
        c                = string(ipos:ipos)
        i                = index(UC,c)
        if( i .gt. 0 ) c = lc(i:i)
        lcstr(ipos:ipos) = c
    enddo
    return
end subroutine hash_lower

!end module hash_functions      
