!-------------------------------------------------------------------------
! This module contains general utility subroutines and functions.
! Some of the general areas covered are:
! 1) File operations, such as 
! 2) String handling commands similar to some of those in Perl
! 3) Basic array operations such as inversions and searches
! 4) Routines for generating combinatorial quantities
! 5) Routines for determining machine specific precisions
!-------------------------------------------------------------------------

Module utils

  Use defaults, Only: strLen, lstrLen, RDbl, COMMENT_CHAR, twopi, xlstrlen, &
      one, dbgflag
#ifdef NAGCOMPILER
  Use f90_unix, Only: getenv
  Use f90_unix_proc, Only: system
#endif

  Implicit None
  Save

  Private
  Public :: allocErrDisplay, isfileopen, split, getpath, filesrchstr, &
      stripcmnt, toreal, toint, toupper, tolower, findint, findstr, &
      genfilename, getlineswpair, getinvangle, makestr, sumlogical, &
      combine, erfc, filetoStringArr, maxn, findgt, getMachPrec, isdigit, &
      filesrchwocomment, deallocErrDisplay, comb, arrnorm, minn, isreal, &
      multarrvec, multvecvec, dispvec, isblank, readblank, int2str, real2str, &
      str2seq,fileSrchStrAll, checkandstop, firstchars, digits, getext, &
      largestofseq, skiplines, cleanstring, getdepth, swap, condenseint, &
      subintlist, unionintlist, chklogical, syscall, str2intmap

  Interface minn
    Module Procedure mini
    Module Procedure minr
  End Interface

  Interface maxn
    Module Procedure maxi
    Module Procedure maxr
  End Interface

  Interface erfc
    Module Procedure erfc_regular
    Module Procedure erfc_and_exp
  End Interface

  Interface inarray
    Module Procedure utils_inarrayi
  End Interface

  Interface filesrchstr
    Module Procedure fileSrchString
    Module Procedure fileSrchStringwCase
    Module Procedure fileSrchStrings
  End Interface

  Interface swap
    Module Procedure realswap
    Module Procedure swap_integers
    Module Procedure swap_subsets
  End Interface

  Interface toreal
    Module Procedure toreal_normal
    Module Procedure toreal_or_error
  End Interface

  Interface toint
    Module Procedure toint_normal
    Module Procedure toint_or_error
  End Interface

  Interface int2str
    Module Procedure int2str_single
    Module Procedure int2str_array
  End Interface

Contains
!----------------------------------------------------------
! String Utilities
!----------------------------------------------------------
  !-----------------------------------
  ! Converts the string "s1" to upper case
  !-----------------------------------
!!$  Character(len=strLen) Function toupper(s1)
  Function toupper(s1)
    Character(len=strLen) :: toupper
    Character(*), Intent(IN) :: s1
    Character(len=strLen)    :: temp
    Character           :: c
    Integer             :: length, LowerToUpper, i
 
    temp = s1
!    temp = Adjustl(s1)    ! Remove any leading blanks
    length = len(temp)
    LowerToUpper = iachar("A") - iachar("a") 

    toupper = ""
    Do i=1, length
       c = temp(i:i)
       If ( c >= "a" .AND. c <= "z") Then
          c = achar(iachar(c) + LowerToUpper)
       Endif
       toupper(i:i) = c
    Enddo
    Return
  End Function toupper

  !-----------------------------------
  ! Converts the string "s1" to lower case
  !-----------------------------------
!!$  Character(len=strLen) Function tolower(s1)
  Function tolower(s1)
    Character(len=strLen) :: tolower
    Character(*), Intent(IN) :: s1
    Character(len=strLen)    :: temp
    Character           :: c
    Integer             :: length, LowerToUpper, i
    
    temp = Adjustl(s1)    ! Remove any leading blanks
    length = len(temp)
    LowerToUpper = iachar("A") - iachar("a") 

    tolower = ""
    Do i=1, length
       c = temp(i:i)
       If ( c >= "A" .AND. c <= "Z") Then
          c = Achar(Iachar(c) - LowerToUpper)
       Endif
       tolower(i:i) = c
    Enddo
    Return
  End Function tolower

  !-----------------------------------
  ! cleans a string by removing comments, trailing and leading blanks
  !-----------------------------------
  Function cleanstring(input)
    Character(len=strLen) :: cleanstring
    Character(*), Intent(in) :: input 
    cleanstring=Adjustl(Trim(stripcmnt(input)))
  End Function cleanstring

  !--------------------------------------------------------------
  ! Strips and returns the first non-digit character of a string
  !--------------------------------------------------------------
  Function firstchars(input)
    Character(len=strLen)        :: firstchars
    Character(*), Intent(In)     :: input

    Integer       :: i

    firstchars = ''

    Do i = 1,len(input)
      If ((input(i:i) == '1').Or. &
          (input(i:i) == '2').Or. &
          (input(i:i) == '3').Or. &
          (input(i:i) == '4').Or. &
          (input(i:i) == '5').Or. &
          (input(i:i) == '6').Or. &
          (input(i:i) == '7').Or. &
          (input(i:i) == '8').Or. &
          (input(i:i) == '9').Or. &
          (input(i:i) == '0')) Return   !** admittedly crude
      Write(firstchars,'(2a)') Trim(firstchars),input(i:i)
    End Do

  End Function firstchars

  !------------------------------------------------------------------
  ! Strips and returns the digit characters of a string AS A STRING
  !------------------------------------------------------------------
  Function digits(input)
    Character(len=strLen)        :: digits
    Character(*), Intent(In)     :: input

    Integer       :: i

    digits = ''

    Do i = 1,len(input)
      If ((input(i:i) == '1').Or. &
          (input(i:i) == '2').Or. &
          (input(i:i) == '3').Or. &
          (input(i:i) == '4').Or. &
          (input(i:i) == '5').Or. &
          (input(i:i) == '6').Or. &
          (input(i:i) == '7').Or. &
          (input(i:i) == '8').Or. &
          (input(i:i) == '9').Or. &
          (input(i:i) == '0')) Then   !** admittedly crude
        Write(digits,'(2a)') Trim(digits),input(i:i)
      End If
    End Do

  End Function digits

  !--------------------------------------------
  ! Returns the number of occurences of "delim"
  ! at the beginning of the string "line". 
  !--------------------------------------------
  Integer Function skip(line, delim)
    Character(*), Intent(in)  :: line
    Character, Intent(IN)     :: delim
    Integer    :: i, n_occur
    
    n_occur=0
    Do i=1,Len(line) 
      If (line(i:i) /= delim) Exit
     n_occur=i
    Enddo
    skip = n_occur
    Return
  End Function skip

  !-----------------------------------------------------
  ! Takes a string and returns the fields in an array
  ! where the fields are delimited using the optional 
  ! field "delim".  Default "delim" is space
  !-----------------------------------------------------
  Integer Function split(line, fields, delim)
    Character(*), Intent(in)                :: line
    Character(*), Dimension(:), Intent(out) :: fields
    Character, Intent(in), Optional         :: delim

    Character                 :: delimeter
    Character(len(line))      :: temp
    Integer     :: length, start, indx, nfields, maxfields, nskip

    ! make sure the fields are all blank
    fields(Lbound(fields,1):Ubound(fields,1)) = ""

    ! The default delimeter is a space
    If (Present(delim)) Then
      delimeter = delim
    Else
      delimeter = " "
    Endif

    ! Get the upper bound on number of fields
    maxfields = Ubound(fields,1)

    ! Remove any leading and trailing spaces
    temp = Trim(Adjustl(line))
    length = Len_trim(temp)
    If (length == 0) Then
      split = 1
      Return
    End If

    start = 1
    nfields = 0
    Do 
      ! Check if we've reached the end of the string
      If (start > length) Exit

      indx = Index(temp(start:length), delimeter)

      ! Make sure we are within bounds
      If (nfields == maxfields) Then
        Write(0,*) 'Size of the array is too small to hold all the fields'
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        split = nfields
        Return
      Endif
      nfields = nfields + 1

      If (indx /= 0) Then
        ! Get the field excluding the delimeter
        fields(nfields) = temp(start:start+indx-2)
        start = start + indx 
        If (start > length) Exit

        ! Skip over multiple delimeter occurrences
        nskip = skip(temp(start:length),delimeter)
        start = start + nskip
      Else
        fields(nfields) = temp(start:length)
        Exit
      End If
    End Do
    split = nfields
    Return
  End Function split

  !-----------------------------------------------------
  ! Takes an array and returns a string 
  ! where the fields are delimited using the optional 
  ! field "delim".  Default "delim" is space
  !-----------------------------------------------------
!!$  Character(len=lstrLen) Function combine(fields, delim)
  Function combine(fields, delim)
    Character(len=lstrLen) :: combine
    Character, Intent(in), Optional         :: delim
    Character(*), Dimension(:), Intent(in) :: fields
    Character                 :: delimeter
    Integer     :: maxfields, i, leng
    Character(len=lstrLen) :: line

    ! make sure the line is blank
    !** Initialize the return line
    combine = ""
    line = ""

    !** The default delimeter is a space
    If (Present(delim)) Then
      delimeter = delim
    Else
      delimeter = " "
    Endif

    !** Get the upper bound on number of fields
    maxfields = Ubound(fields,1)

    !** Loop through the fields and add them to the line
    line = trim(fields(1))
    Do i = 2, maxfields
      leng = len_trim(line)
      line = line(1:leng)//delimeter//Trim(fields(i))
    End Do

    !** Trim any trailing spaces
    combine = Trim(line)

    Return
  End Function combine

  !------------------------------------------------------------------
  ! This function checks to see if a character string is "blank".  
  ! It will return a false value if there is anything but spaces, 
  ! tabs or '-' characters in the string.
  !------------------------------------------------------------------
  Logical Function isblank(line)
    Character(*), Intent(In)           :: line

    Integer                            :: i

    isblank = .True.

    Do i = 1,Len(Trim(line))
      If ((line(i:i) /= ' ').And.(line(i:i) /= '-').And. &
          (Ichar(line(i:i)) /= 9)) Then
        isblank = .False.
        Return
      End If
    End Do

  End Function isblank

  !------------------------------------------------------------------
  ! This subroutine reads a line that is expected to be blank.  It
  ! Returns an error message and stops if this is not the case.
  !------------------------------------------------------------------
  Subroutine readblank(unit,filename,lineno,errormsg)
    Integer, Intent(In)                :: unit
    Character(*), Intent(In)           :: filename
    Integer, Intent(In)                :: lineno
    Character(*), Intent(In), Optional :: errormsg

    Character(len=255) :: line

    Read(unit,'(a)') line

    If (.Not.isblank(line)) Then
      Write(0,'(1x,2a,i4)') 'ERROR: expected to read a blank line at: ', &
          Trim(filename),lineno
      Write(0,'(2x,a)') 'Instead, read this: ',Trim(line)
      If (Present(errormsg)) Then
        Write(0,'(2x,a)') Trim(errormsg)
      End If
      Stop
    End If

  End Subroutine readblank

  !-------------------------------------------------------
  ! This function strips the "line" of the comment string
  ! which is defined by an optional argument "optcomment"
  ! The default is COMMENT_CHAR
  !-------------------------------------------------------
!!$  Character(len=lstrLen) Function stripcmnt(line, optcomment)
  Function stripcmnt(line, optcomment)
    Character(len=lstrLen)             :: stripcmnt
    Character(*), Intent(in)           :: line
    Character, Optional, Intent(in)    :: optcomment

    Character                             :: commentchar
    Character(len=lstrLen), Dimension(50) :: newline
    Integer                               :: nfields

    If (Present(optcomment)) Then
      Write(0,*) __FILE__,__LINE__,line
      Write(0,*) __FILE__,__LINE__,optcomment
    End If

    If (Present(optcomment)) Then
      commentchar = optcomment
    Else
      commentchar = COMMENT_CHAR
    End If

    newline = ""
    nfields = split(line, newline, commentchar)
    stripcmnt = newline(1)

  End Function stripcmnt

  !-----------------------------------------------------------------------
  ! Takes a string and converts it to a real
  ! Requires:  strfield -- the string
  !            opterror -- optional error to return if there's a problem
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function toreal_normal(strfield, opterror)
    Character(*), Intent(In)       :: strfield
    Integer, Optional, Intent(Out) :: opterror

    Integer                 :: error

    Read(strfield, *, IOSTAT=error) toreal_normal
    If (Present(opterror)) Then
      opterror = error
    Else
      If (error /= 0) Then
        Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(0,'(a)') ' Error during string to real conversion'
        Write(0,'(3a)') '  passed string was: " ',Trim(strfield), '"'
        Write(0,'(a,i4)') '  error number: ',error 
        Stop
      End If
    End If

  End Function toreal_normal

  !-----------------------------------------------------------------------
  ! Takes a string and converts it to a real, if there's an error dump 
  ! message with string to standard error
  ! Requires:  strfield -- the string
  !            opterror -- error message
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function toreal_or_error(strfield, opterror)
    Character(*), Intent(In)       :: strfield, opterror

    Integer      :: error

    Read(strfield, *, IOSTAT=error) toreal_or_error
    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)') ' Error during string to real conversion'
      Write(0,'(3a)') '  passed string was: " ',Trim(strfield), '"'
      Write(0,'(3a)') '  Conversion was attempted for "',Trim(opterror),'"'
      Stop
    End If

  End Function toreal_or_error

  !--------------------------------------------------------------
  ! Takes a string and converts it to an integer, if there's an 
  ! error, optionally return error number
  ! Requires:  strfield -- the string
  !            opterror -- integer error number
  !--------------------------------------------------------------
  Integer Function toint_normal(strfield, opterror)
    Character(*), Intent(In)       :: strfield
    Integer, Optional, Intent(Out) :: opterror

    Integer               :: error

    Read(strfield, *, IOSTAT=error) toint_normal

    If (Present(opterror)) Then
      opterror = error
    Else
      If (error /= 0) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            ' toint error'
        Write(0,'(3a)') " input string was '",Trim(strfield),"'"
        Stop
      End If
    End If

  End Function toint_normal

  !--------------------------------------------------------------
  ! Takes a string and converts it to an integer, if there's an 
  ! error dump message with string to standard error
  ! Requires:  strfield -- the string
  !            opterror -- error message
  !--------------------------------------------------------------
  Integer Function toint_or_error(strfield, opterror)
    Character(*), Intent(In)       :: strfield, opterror

    Integer      :: error

    Read(strfield, *, IOSTAT=error) toint_or_error

    If (error /= 0) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(a)') ' Error during string to integer conversion'
      Write(0,'(3a)') '  passed string was: " ',Trim(strfield), '"'
      Write(0,'(3a)') '  Conversion was attempted for "',Trim(opterror),'"'
      Stop
    End If

  End Function toint_or_error

  !-----------------------------------------------------
  ! Takes an integer and returns a string.  This is 
  ! useful because it can then be trimmed
  !-----------------------------------------------------
  Function int2str_single(int)
    Character(len=strLen)          :: int2str_single
    Integer, Intent(In)            :: int

    Character(len=strLen)          :: intstring

    Write(intstring,'(i30)') int
    intstring = Adjustl(intstring)
    int2str_single = intstring

  End Function int2str_single

  !------------------------------------------------------------
  ! Takes an integer ARRAY and returns a string.  This is 
  ! useful because it can then be trimmed
  !------------------------------------------------------------
  Function int2str_array(int)
    Character(len=strLen)                :: int2str_array
    Integer, Dimension(:), Intent(In)    :: int

    Integer                        :: i
    Character(len=lstrLen)         :: intstring

    Write(intstring,'(i30)') int(1)
    int2str_array = Adjustl(intstring)
    Do i = 2,Size(int)
      Write(intstring,'(i30)') int(i)
      intstring = Adjustl(intstring)
      Write(int2str_array,'(3a)') Trim(int2str_array),',',Trim(intstring)
    End Do

  End Function int2str_array

  !-------------------------------------------------------------------------
  ! Takes a real and returns a string.  This is useful because it can then 
  ! be trimmed.  It also switches to scientific notation if practical.  If
  ! a suggested number of character is passed, it will either restrict the
  ! length to that number of characters or increase the number of characters
  ! as necessary to make the display reasonable.  
  ! Requires:  number -- real number to convert to string
  !            sugnchar -- suggested number of characters in output, optional
  ! NOTE: Because it is a fairly complex routine, it should not be used near
  ! inner loops as it will slow down the code.  
  !-------------------------------------------------------------------------
  Function real2str(number,sugnchar)
    Character(len=strLen)          :: real2str
    Real(kind=RDbl), Intent(In)    :: number
    Integer, Intent(In), Optional  :: sugnchar

    Integer                        :: nchar,xtra,decimal_places,minchar,beforept
    Integer                        :: mindecimals,minechar,wantnchar
    Logical                        :: exp_form
    Character(len=10)              :: form
    Character(len=strLen)          :: intstring1,intstring2
    Real(kind=RDbl)                :: logn,abslogn

    !** set defaults
    exp_form = .False.
    wantnchar = 12   !** default desired number of characters
    mindecimals = 1  !** minimum number of decimal places behind the point

    !** Use the suggested number of characters if present
    If (Present(sugnchar)) wantnchar = sugnchar

    !** Get log10 of the number
    If (number == 0.0_RDbl) Then
      logn = 1
    Else
      logn = Log10(Abs(number))
    End If
    abslogn = Abs(logn)

    !** Get the minimum required number of characters to avoid overflow
    !** this applies only to the normal, non-exponential form
    xtra = 0
    If (number < 0) xtra = xtra + 1  !** add space for minus sign
    beforept = Max(0,Int(logn)) + 1
    minchar = beforept + xtra + 1 + mindecimals

    !** Get the minimum required number of characters for exponential form
    minechar = 7     !** minimum characters needed for exp form (0.0E+00)
!    If (number < 0) minechar = minechar + 1  !** add space for minus sign    

    !** Increase the 'wanted' number of characters if unreasonable
    If (minchar > wantnchar) Then !** won't work as normal notation
      If (wantnchar < minechar) Then !** need to increase wanted number
        If (minchar <= minechar) Then !** increase to use normal form
          wantnchar = minchar
        Else   !** increase to use exp form
          exp_form = .True. 
          wantnchar = minechar
        End If
      Else  !** use exp form to fit into desired nchars
        exp_form = .True.        
      End If
    End If

    !** Calculate decimal places and wwitch to exp form if number is too small
    If (.Not. exp_form) Then
      !** Calculate decimal places for normal form
      decimal_places = wantnchar - beforept - 1 - xtra

      If ((decimal_places + logn) < 0) Then
        exp_form = .True.
        wantnchar = minechar
      End If
    End If

    nchar = wantnchar

    !** set the format    
    If (exp_form) Then
      decimal_places = Max((nchar-6),mindecimals)
      intstring1 = int2str(nchar)
      intstring2 = int2str(decimal_places)
      Write(form,'(5a)') '(e',Trim(intstring1),'.', &
          Trim(intstring2),')'
    Else
      intstring1 = int2str(nchar)
      intstring2 = int2str(decimal_places)
      Write(form,'(5a)') '(f',Trim(intstring1),'.', &
          Trim(intstring2),')'
    End If

    Write(real2str,form) number

    real2str = Adjustl(real2str)

  End Function real2str

  !----------------------------------------------------------------------------
  ! Converts a string into an integer map.  Format of the string is expected 
  ! to be (for example):  1-3,5-7,4@31-37  or:  1-3@31-33;5-7@34-36;4@37
  ! integers 1-3 in the first set are mapped to 31-33, 5-7 to 34-36 and 4 to 37.
  ! Output is structured such that num(integer_set_1) = integer_set_2.  It
  ! returns the number of integer--integer mappings.
  ! Requires:  string -- string to process
  !            nums -- output array of integers in map form
  !----------------------------------------------------------------------------
  Integer Function str2intmap(string,nums)
    Character(*), Intent(In)           :: string
    Integer, Dimension(:), Intent(Out) :: nums

    Integer                    :: i,j,n,nmaps,biggestn,nset1,nset2
    Integer, Dimension(Size(nums))           :: set1,set2
    Character(len=strLen), Dimension(20)     :: set1strs,set2strs
    Character(len=strLen), Dimension(strLen) :: fields,subfields

    !** Set defaults
    str2intmap = 0
    nums = -1

    !** Split the full string into separate maps (;-separated)
    nmaps = split(string,fields,';')

    biggestn = -1
    Do i = 1,nmaps
      !** Split each map string into separate specification (@-separated)
      n = split(fields(i),subfields,'@')
      If (n > 2) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' a ";"-separated field should not have more than one "@"'
        Write(0,'(2a)') 'String is: ',Trim(fields(i))
        Stop
      End If

      !** Store the map specifications
      set1strs(i) = subfields(1)
      set2strs(i) = subfields(2)

      !** Get the largest number in the first half of the specifications
      n = largestofseq(set1strs(i))
      If (n > biggestn) biggestn = n
    End Do

    !** Check size of incoming array
    If (nmaps > Size(nums)) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          " passed array too small"
      Stop
    End If

    !** Create the mapping
    Do i = 1,nmaps
      nset1 = str2seq(set1strs(i),set1)
      nset2 = str2seq(set2strs(i),set2)

      If (nset1 /= nset2) Then
        Write(0,'(1x,2a,i4,a,i2,a)') __FILE__," : ",__LINE__, &
            ' number of integers on eiter side of map ',i,' do not match'
        Write(0,'(2a)') 'String is: ',Trim(string)
        Stop
      End If

      Do j = 1,nset1
        nums(set1(j)) = set2(j)
      End Do

      str2intmap = str2intmap + 1
    End Do

  End Function str2intmap

  !---------------------------------------------------------------------------
  ! Converts a string containing a sequence of integers into an array of 
  ! integers.  Format of the string is expected to be:  3,5-8,12 
  ! (for example).  This would be converted to: 3,5,6,7,8,12
  ! Requires: string -- string to process
  !           nums -- output array of integers
  !-----------------------------------------------------------------------------
  Integer Function str2seq(string,nums)
    Character(*), Intent(In)           :: string
    Integer, Dimension(:), Intent(Out) :: nums

    Integer                    :: i,j,n,n1,n2,nfields,nsubfields
    Character(len=strLen), Dimension(strLen) :: fields,subfields

    n = 0
    str2seq = -1
    nums = -1

    !** split the string into comma separated chunks, then process
    nfields = split(string,fields,',')
    Do i = 1,nfields
      nsubfields = split(fields(i),subfields,'-')
      Select Case (nsubfields)
      Case(0) 
        Cycle
      Case(1)
        n = n+1
        If (n > Size(nums)) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              " passed array too small"
          Stop
        End If
        nums(n) = toint(subfields(1),'subfield in str2seq utility')
      Case(2)
        n1 = toint(subfields(1),'subfield in str2seq utility')
        n2 = toint(subfields(2),'subfield in str2seq utility')
        Do j = n1,n2
          n = n + 1
          If (n > Size(nums)) Then
            Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
                " passed array too small"
            Stop
          End If
          nums(n) = j
        End Do
      Case Default
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            " Unprocessable string: ",Trim(string)
        Stop
      End Select
    End Do

    str2seq = n
    
  End Function str2seq

  !---------------------------------------------------------------------------
  ! Gets the largest number in a string sequence of numbers
  ! Format of the string is expected to be:  3,5-8,12 
  ! Requires:  string -- string to process
  !-----------------------------------------------------------------------------
  Integer Function largestofseq(string)
    Character(*), Intent(In)           :: string

    Integer                    :: i,j,n,nfields,nsubfields
    Character(len=strLen), Dimension(strLen) :: fields,subfields

    largestofseq = -1

    !** split the string into comma separated chunks, then process
    nfields = split(string,fields,',')
    Do i = 1,nfields
      nsubfields = split(fields(i),subfields,'-')
      Do j = 1,nsubfields
        n = toint(subfields(j))
        If (n > largestofseq) largestofseq = n
      End Do
    End Do
    
  End Function largestofseq

  !---------------------------------------------------------------------------
  ! Gets the smallest number in a string sequence of numbers
  ! Format of the string is expected to be:  3,5-8,12 
  ! Requires:  string -- string to process
  !-----------------------------------------------------------------------------
  Integer Function smallestofseq(string)
    Character(*), Intent(In)           :: string

    Integer                    :: i,j,n,nfields,nsubfields
    Character(len=strLen), Dimension(strLen) :: fields,subfields

    smallestofseq = 1e5

    !** split the string into comma separated chunks, then process
    nfields = split(string,fields,',')
    Do i = 1,nfields
      nsubfields = split(fields(i),subfields,'-')
      Do j = 1,nsubfields
        n = toint(subfields(j))
        If (n < smallestofseq) smallestofseq = n
      End Do
    End Do
    
  End Function smallestofseq

  !---------------------------------------------------
  ! Checks if a character is a digit
  !---------------------------------------------------
  Logical Function isdigit(char)
    Character, Intent(in) :: char
    
    If (char < '9' .And. char > '0') Then
      isdigit = .True.
    Else
      isdigit = .False.
    End If
  End Function isdigit

  !-----------------------------------------------------------------
  ! Checks to see if a character string can be converted to a real
  !-----------------------------------------------------------------
  Logical Function isreal(char)
    Character(*), Intent(In) :: char

    Integer         :: error
    Real(kind=RDbl) :: output

    isreal = .False.
    output = toreal(char,error)
    If (error == 0) isreal = .True.

  End Function isreal

  !------------------------------------------------------
  ! Generates a compact string given two numbers "num1"
  ! and "num2" and a string in between.  Useful when
  ! you want to print out mathematical operations explicitly
  ! such as "num1/num2"
  !------------------------------------------------------
!!$  Character(len=strLen) Function makestr(num1, str, num2)
  Function makestr(num1, str, num2)
    Character(len=strLen) :: makestr
    Integer, Intent(in)   :: num1, num2
    Character(*), Intent(in)  :: str

    Character(len=strLen)     :: str1, str2
    
    Write(str1,'(i12)') num1
    Write(str2,'(i12)') num2
    makestr = Trim(Adjustl(str1))//Trim(str)//Trim(Adjustl(str2))
    Return
  End Function makestr
  
!----------------------------------------------------------
! File Utilities
!----------------------------------------------------------

  !-----------------------------------------------------------------
  ! Check if a file has "tabs" in it. The ascii value of "tab" is 9.
  ! The function returns the line number where the tab was found 
  ! otherwise it returns a zero.
  !-----------------------------------------------------------------
  Integer Function check_tabs(filename)
    Character(*), Intent(in) :: filename

    Character(len=lstrLen) :: inputline
    Integer                :: i, lineno, error
    
    check_tabs = 0
    lineno = 0

    Open(unit=17, file=filename)
    Do
      Read(17, '(a)', Iostat=error) inputline
      If (error /= 0) Exit
      
      lineno = lineno + 1
      Do i=1, Len(Trim(inputline))
        If (Ichar(inputline(i:i)) == 9) Then
          check_tabs = lineno
          Exit
        End If
      End Do
    End Do

    Close(unit=17)
  End Function check_tabs

  !-----------------------------------------------------
  ! This function skips "nlines" in the file connected to
  ! "unitno"
  !-----------------------------------------------------
  Subroutine skiplines(unitno, nlines)
    Integer, Intent(in) :: unitno, nlines

    Integer        :: i

    Rewind(unitno)
    Do i=1, nlines
      Read(unitno,*)
    End Do
  End Subroutine skiplines

  !----------------------------------------------------
  ! This function returns a string with the simulation
  ! no. "simno" appended to the filename "filename"
  ! Requires: filename -- base string to add to
  !           simno -- integer to tack to the end
  !           secondno -- optional second number
  !----------------------------------------------------
!!$  Character(len=strLen) Function genfilename(filename, simno, secondno)
  Function genfilename(filename, simno, secondno)
    Character(len=strLen)          :: genfilename
    Character(*), Intent(In)       :: filename
    Integer, Intent(In)            :: simno
    Integer, Intent(In), Optional  :: secondno

    Character(len=strLen) :: i

    !** This is a hack to get it to work with the Lahey compiler
    i = int2str(simno)

    Write(genfilename, '(3a)') Trim(filename),'.',Trim(i)

    If (Present(secondno)) Then
      i = int2str(secondno)
      Write(genfilename, '(3a)') Trim(genfilename),'.',Trim(i)
    End If

  End Function genfilename

  !--------------------------------------------
  ! Get the pathname stored in the environment
  ! variable "s1"
  !--------------------------------------------
!!$  Character(len=lstrLen) Function getpath(s1)
  Function getpath(s1)
    Character(len=lstrLen) :: getpath
    Character(*), Intent(IN) :: s1
    Character(len=lstrLen)    :: path
    
    Call getenv(s1,path)
    If (path == ' ') Then
      Write(0,*)
      Write(0,'(2a,i4)') __FILE__," : ",__LINE__
      Write(0,'(3a)') 'Please set your ',s1,' environment variable'
      Stop
    Endif

    getpath = path
    If (path(len(path):len(path)) /= '/') Then
      getpath = trim(path) // '/'
    Endif
    Return
  End Function getpath

  !-----------------------------------------------------
  ! Returns the extension (.something) of a filename
  ! Requires:  filename -- filename input
  !-----------------------------------------------------
  Function getext(filename)
    Character(len=strLen) :: getext
    Character(*)          :: filename

    Integer       :: nfields
    Character(len=lstrLen), Dimension(strLen) :: fields

    nfields = split(filename,fields,'.')
    getext = fields(nfields)

  End Function getext

  !--------------------------------------------------------------------
  ! Takes a string and a file unit number and returns the line in the
  ! associated text file with that string.  The search has failed if
  ! return value is 0.  If the string is found the return value is the
  ! line number where the string was found.
  ! Requires:  unitno -- unit number to search
  !            srchstr -- string to look for
  !            line -- returned line
  !            rewindfile -- optional flag to require rewinding file
  !--------------------------------------------------------------------
  Integer Function fileSrchString(unitno, srchstr, line, rewindfile)
    Integer, Intent(In)             :: unitno
    Character(*), Intent(In)        :: srchstr
    Character(*), Intent(Out)       :: line
    Logical, Intent(In), Optional   :: rewindfile

    Character(len=2000)      :: str
    Integer                  :: ios, indx

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchString = 0    
    Do 
      fileSrchString = fileSrchString + 1
      Read(unitno,'(a)', IOSTAT=ios) str

      !** Check for end-of-file
      If (ios /= 0) Then
        fileSrchString = 0
        line = ' '
        Exit
      End If

      indx = Index(toupper(str), Trim(toupper(srchstr)), .False.)
      If (indx /= 0) Then
        line = str
        Exit
      End If
    End Do

  End Function fileSrchString

  !--------------------------------------------------------------------
  ! Same as fileSrchString, but has an option for using ToUpper or not.
  ! Eventually, this should replace fileSrchString.
  !
  ! Takes a string and a file unit number and returns the line in the
  ! associated text file with that string.  The search has failed if
  ! return value is 0.  If the string is found the return value is the
  ! line number where the string was found.
  ! Requires:  unitno -- unit number to search
  !            srchstr -- string to look for
  !            nocase -- True => use ToUpper, no case sensitivity
  !            line -- returned line
  !            rewindfile -- optional flag to require rewinding file
  !--------------------------------------------------------------------
  Integer Function fileSrchStringwCase(unitno, srchstr, nocase, line, rewindfile)
    Integer, Intent(In)             :: unitno
    Character(*), Intent(In)        :: srchstr
    Logical, Intent(In)             :: nocase
    Character(*), Intent(Out)       :: line
    Logical, Intent(In), Optional   :: rewindfile

    Character(len=2000)      :: str
    Integer                  :: ios, indx

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchStringwCase = 0    
    Do 
      fileSrchStringwCase = fileSrchStringwCase + 1
      Read(unitno,'(a)', IOSTAT=ios) str

      !** Check for end-of-file
      If (ios /= 0) Then
        fileSrchStringwCase = 0
        line = ' '
        Exit
      End If

      If (nocase) Then
        indx = Index(toupper(str), Trim(toupper(srchstr)), .False.)
      Else
        indx = Index(str, Trim(srchstr), .False.)
      End If

      If (indx /= 0) Then
        line = str
        Exit
      End If
    End Do

  End Function fileSrchStringwCase

  !---------------------------------------------------
  ! Same as fileSrchString, but first makes sure that 
  ! the line doesn't contain a comment
  !---------------------------------------------------
  Integer Function fileSrchwoComment(unitno, srchstr, line, rewindfile)
    Integer, Intent(IN)             :: unitno
    Character(*), Intent(IN)        :: srchstr
    Character(*), Intent(OUT)       :: line
    Logical, Intent(IN), Optional   :: rewindfile
    Character(len=255)              :: str
    Integer                         :: ios, indx

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchwoComment = 0    
    Do 
      fileSrchwoComment = fileSrchwoComment + 1
      Read(unitno,'(a)', IOSTAT=ios) str
      str = stripcmnt(str)  !**the only difference compared to fileSrchString
      !** Check for end-of-file
      If (ios /= 0) Then
        fileSrchwoComment = 0
        line = ' '
        Exit
      End If

      indx = Index(toupper(str), Trim(toupper(srchstr)), .False.)
      If (indx /= 0) Then
        line = str
        Exit
      End If
    End Do
  End Function fileSrchwoComment

  !----------------------------------------------------------------------------
  ! Searches file unitno for a line that contains all the strings stored in
  ! the character array srchstrs. If found, it returns the line and line no.,
  ! otherwise, it returns 0. Optionally, it will rewind to the start of the 
  ! file if rewindfile is True.
  !----------------------------------------------------------------------------
  Integer Function fileSrchStrings(unitno, srchstrs, line, rewindfile)
    Integer, Intent(IN)             :: unitno
    Character(*), Intent(IN), Dimension(:) :: srchstrs
    Character(*), Intent(OUT)       :: line
    Logical, Intent(IN), Optional   :: rewindfile

    Character(len=255)              :: str
    Character(len=strLen), Dimension(strLen) :: fields
    Integer                         :: ios, nstr, i, j, k, nfields
    Logical :: found
    Integer, Dimension(Size(srchstrs,1)) :: foundpos

    nstr = Size(srchstrs,1)

    If (Present(rewindfile)) Then
      If (rewindfile) Rewind(unitno)
    End If

    fileSrchStrings = 0    
    Do 
      fileSrchStrings = fileSrchStrings + 1
      Read(unitno,'(a)', IOSTAT=ios) str
      ! Check for end-of-file
      If (ios /= 0) Then
        fileSrchStrings = 0
        line = ' '
        Exit
      End If

      nfields = split(str,fields)
      !** Check if all the strings exist
      Do i = 1, nstr
        found = .False.
        Do j = 1, nfields
          If (Trim(tolower(srchstrs(i))) == Trim(tolower(fields(j)))) Then
            ! We want a one-to-one mapping from elements of "srchstrs" to 
            ! element of "fields".  So, make sure that the field "j" does
            ! not already have a mapping
            found = .True.
            Do k=1, i-1
              If (j == foundpos(k)) Then
                found = .False.
                Exit
              End If
            End Do
            ! If found is still true then we have truly found a match
            ! and we can move on to the next field
            If (found) Then
              foundpos(i) = j
              Exit
            End If
          End If
        End Do
        If (.Not.found) Exit
      End Do
      If (found) Then 
        line = str
        Exit
      End If
    End Do
  End Function fileSrchStrings

  !----------------------------------------------------------------------------
  ! Takes a string and a file unit number and returns the number of lines in
  ! the file that contain this string.  It also optionally returns the line
  ! numbers at which the strings occur.
  ! Requires: unitno -- unit number of an open file
  !           srchstr -- search string
  !           chkcomments -- if True, will ignore '#'-commented strings
  !           lines -- array of matching line numbers, optional
  !----------------------------------------------------------------------------
  Integer Function fileSrchStrAll(unitno,srchstr,chkcomments,lines)
    Integer, Intent(In)                           :: unitno
    Character(*), Intent(In)                      :: srchstr
    Logical, Intent(In)                           :: chkcomments
    Integer, Dimension(:), Intent(Out), Optional  :: lines

    Integer                         :: nmatch, ios, indx, length, line
    Character(len=4)                :: first_pass_str
    Character(len=255)              :: str

    length = Len(Trim(srchstr))
    If (length > 4) Then
      first_pass_str = srchstr(1:4)
    Else
      first_pass_str = srchstr
    End If

    !** rewind file
    Rewind(unitno)

    nmatch = 0 
    line = 0

    Do 
      Read(unitno,'(a)', IOSTAT=ios) str
!LC      Write(*,*) Trim(str)

      !** Check for end-of-file
      If (ios /= 0) Exit
      line = line + 1

      !** look for a shortened version first to speed (?) things up
      indx = Index(str,first_pass_str,.False.)
      If (indx /= 0) Then
        If (chkcomments) str = stripcmnt(str)
        indx = Index(str,Trim(srchstr),.False.)
        If (indx /= 0) Then
          nmatch = nmatch + 1
          If (Present(lines)) lines(nmatch) = line
        End If
      End If

    End Do

    fileSrchStrAll = nmatch

  End Function fileSrchStrAll

  !---------------------------------------------------------------------------
  ! Checks to see if the file "filename" is open.  If it is then the 
  ! associated unit number is returned else a negative number is returned.
  !---------------------------------------------------------------------------
  Integer Function isfileopen(filename)
    Character(*), Intent(in)       :: filename
    Integer      :: unitno
    Logical      :: connected

    Inquire(File=trim(filename), Opened=connected, Number=unitno)
    If (connected) Then
      isfileopen = unitno
    Else
      isfileopen = -1
    End If
    Return
  End Function isfileopen

  !---------------------------------------------------
  ! Gets the no. of characters in a file
  !---------------------------------------------------
  !---------------------------------------------------
  ! Returns a array of integer with the index holding
  ! the no. of lines in the file and the second index
  ! the no. of characters in the file "unitno"
  !---------------------------------------------------
  Function getfilestats(unitno)
    Integer, Intent(in) :: unitno
    Integer, Dimension(2) :: getfilestats
    
    Character(150)   :: inputbuffer
    Integer          :: nlines, nchars, error
    
    Rewind(unitno)
    nchars = 0
    nlines = 0
    Do 
      nlines = nlines + 1
      Read(unitno,'(a)',IOSTAT=error) inputbuffer
      Write(*,'(i8,":",a)') nlines, Trim(inputbuffer)
      If (error /= 0) Exit
      nchars = nchars + Len(Trim(inputbuffer))
    End Do
    
    getfilestats = (/nlines, nchars/)
  End Function getfilestats

  !----------------------------------------------------------------------------
  ! Picks lines out of a file that have a certain pair of names in their
  ! first two words (order irrelevant).  File must already be open!
  !----------------------------------------------------------------------------
  Subroutine getlineswpair(unit,name1,name2,delim,nlines,lines)
    Integer, Intent(IN)                      :: unit
    Character(*), Intent(IN)                 :: name1,name2
    Character(len=1), Intent(IN)             :: delim
    Integer, Intent(OUT)                     :: nlines
    Character(*), Dimension(:), Intent(OUT)  :: lines

    Integer                                  :: lineno
    Integer                                  :: nfields
    Character(len=strLen), Dimension(10)     :: fields
    Integer                                  :: name1field,name2field
    Character(len=2*lstrLen)                 :: line

    nlines = 0
    Do 
      lineno = filesrchstr(unit,name1,line,.FALSE.)
      If (lineno == 0) Then
        Exit
!        Write(0,'(1x,2a,i4, 6a)') __FILE__," : ",__LINE__, &
!            " Could not find the pair ", Trim(name1)," ",Trim(name2), &
!            " in file ", Trim(filename)
      Else
        nfields = split(line,fields,delim)
        !** Check if the name is in field1 or field2
        If (Trim(fields(1)) == Trim(name1)) Then
          name1field = 1
          name2field = 2
        Else
          name1field = 2
          name2field = 1
        End If
        If (Trim(fields(name2field)) /= Trim(name2)) Cycle
      End If

      nlines = nlines + 1
      If (nlines > Size(lines)) Then
        Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
            " passed string array too small"
        Stop
      End If
      lines(nlines) = line
      
    End Do
    Return
  End Subroutine getlineswpair

!---------------------------------------------------------------
! Array Utilities  
!---------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Searches an integer array to find the first instance of number. 
  ! Returns the position or 0 If not found. If the optional parameter 
  ! all is passed, Then the position of all occurances of number in 
  ! list are returned in all
  !----------------------------------------------------------------------------
  Integer Function findint(list,number,all)
    Integer, Dimension(:), Intent(In) :: list
    Integer, Intent(In) :: number
    Integer, Dimension(:), Intent(out), Optional :: all
    Integer :: i, j
    
    j = 0

    If (Size(list) < 1) Then
      findint = 0
      Return
    End If
    Do i = 1,Size(list)
      If (number == list(i)) Then
        findint = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
        Else
          Return
        End If
      End If
    End Do
    If (j == 0) findint = 0
  End Function findint

  !----------------------------------------------------------------------------
  ! Removes repeats from an integer array, returns the new number of integers
  ! in the complete array.
  ! Requires:  list -- integer array
  !----------------------------------------------------------------------------
  Integer Function condenseint(list)
    Integer, Dimension(:), Intent(InOut) :: list

    Integer :: i,j,n,newn

    !** run backwards through the array, for each item, check if there
    !** is a matching item between 1->i-1, if so, remove and condense
    n = Size(list)
    newn = n
    Do i = n,1,-1
      !** Check for matches in indices 1 -> i-1 and remove i if yes
      If (findint(list(1:i-1),list(i)) /= 0) Then
        newn = newn - 1
        Do j = i,n-1
          list(j) = list(j+1)
        End Do
      End If
    End Do

    condenseint = newn

  End Function condenseint

  !----------------------------------------------------------------------------
  ! Removes all integers from the 1st list that match those in the second
  ! list.  Used for extracting subsets from an integer list.  Returns the
  ! new number of integers in the 1st list.
  ! Requires:  list1 -- 1st integer array
  !            list2 -- 2nd integer array
  !----------------------------------------------------------------------------
  Integer Function subintlist(list1,list2)
    Integer, Dimension(:), Intent(InOut) :: list1
    Integer, Dimension(:), Intent(In)    :: list2

    Integer :: i,j,n,newn

    !** Check for duplication in list2, remove from list1 if so
    n = Size(list1)
    newn = n
    Do i = n,1,-1
      !** Check for matches in list2, remove and condense if present
      If (findint(list2,list1(i)) /= 0) Then
        newn = newn - 1
        Do j = i,n-1
          list1(j) = list1(j+1)
        End Do
      End If
    End Do

    subintlist = newn

  End Function subintlist

  !----------------------------------------------------------------------------
  ! Performs a union of two integer lists.  Returns the number of common
  ! integer and the resulting list.
  ! Requires:  list1 -- 1st integer array
  !            list2 -- 2nd integer array
  !            results -- resulting integer array
  !----------------------------------------------------------------------------
  Integer Function unionintlist(list1,list2,results)
    Integer, Dimension(:), Intent(In)   :: list1,list2
    Integer, Dimension(:), Intent(Out)  :: results

    Integer :: i,n,newn,maxn

    maxn = Size(results)

    !** Loop through list1, look for matches in list2 and put in results
    n = Size(list1)
    newn = 0
    Do i = 1,n
      !** Check for matches in list2, remove and condense if present
      If (findint(list2,list1(i)) /= 0) Then
        newn = newn + 1
        If (newn > maxn) Then
          Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
              ' Passed results array not large enough'
          Stop          
        End If
        results(newn) = list1(i)
      End If
    End Do

    unionintlist = newn

  End Function unionintlist

  !----------------------------------------------------------------------------
  ! Searches a STRING (len=strLen) array to find the first instance of 
  ! a given string.  Returns the position or 0 If not found. If the 
  ! optional parameter all is passed, then the position of all occurances 
  ! of this string in list are returned in all
  !----------------------------------------------------------------------------
  Integer Function findstr(list,srchstr,all)
    Character(*), Dimension(:), Intent(In)         :: list
    Character(*), Intent(In)                       :: srchstr
    Integer, Dimension(:), Intent(Out), Optional   :: all

    Integer :: i, j
        
    j = 0

    If (Size(list) < 1) Then
      findstr = 0
      Return
    End If

    Do i = 1,Size(list)
!      Write(*,*) __FILE__,i,Trim(srchstr),Trim(list(i))
      If (Trim(srchstr) == Trim(list(i))) Then
        findstr = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
        Else
!          Write(*,*) 'found it',__FILE__,i,Trim(srchstr),Trim(list(i))
          Return
        End If
      End If
    End Do

    If (j == 0) findstr = 0

  End Function findstr

  !----------------------------------------------------------------------------
  ! Returns the 'depth' of a system subset specification.  This is the level
  ! at which the specification ends.  Returned depths and examples of subsets:
  !   0 = full system      subset = (/0,0,0/)
  !   1 = species level    subset = (/a,0,0/)
  !   2 = molecule level   subset = (/a,b,0/)
  !   3 = atom level       subset = (/a,b,c/)
  ! Requires:  subset -- subset specification
  !----------------------------------------------------------------------------
  Integer Function getdepth(subset)
    Integer, Dimension(3), Intent(In)        :: subset

    Integer         :: i

    Do i = 3,0,-1
      getdepth = i
      If (i == 0) Exit
      If (subset(i) /= 0) Exit
    End Do    

  End Function getdepth

  !----------------------------------------------------------------------------
  ! Searches an integer array to find the first instance greater than number. 
  ! Returns the position or 0 If not found. If the optional parameter 
  ! all is passed, Then the position of all occurances of number in 
  ! list are returned in all and the number of instances found 
  !----------------------------------------------------------------------------
  Integer Function findgt(list,number,all)
    Integer, Dimension(:), Intent(In) :: list
    Integer, Intent(In) :: number
    Integer, Dimension(:), Intent(out), Optional :: all
    Integer :: i, j
    
    j = 0

    If (Size(list) < 1) Then
      findgt = 0
      Return
    End If
    Do i = 1,Size(list)
      If (number < list(i)) Then
        findgt = i
        If (Present(all)) Then
          j = j + 1
          all(j) = i
          findgt = j
        Else
          Return
        End If
      End If
    End Do
    If (j == 0) findgt = 0
  End Function findgt

  !---------------------------------------------
  ! find minimum in an integer or real array
  !---------------------------------------------
  Integer Function mini(list)
    Integer, Dimension(:) :: list
    Integer :: i

    mini = list(1)
    Do i = 2,Size(list)
      If (list(i) < mini) Then
        mini = list(i)
      End If
    End Do
  End Function mini

  Real(kind=RDbl) Function minr(list)
    Real(kind=RDbl), Dimension(:) :: list
    Integer :: i
    
    minr = list(1)
    Do i = 2,Size(list)
      If (list(i) < minr) Then
        minr = list(i)
      End If
    End Do
  End Function Minr

  !-----------------------------------------------
  ! finds the maximum in an array of ints or reals
  !-----------------------------------------------
  Integer Function maxi(list)
    Integer, Dimension(:) :: list
    Integer :: i

    maxi = list(1)
    Do i = 2,Size(list)
      If (list(i) > maxi) Then
        maxi = list(i)
      End If
    End Do
  End Function maxi

  Real(kind=RDbl) Function maxr(list)
    Real(kind=RDbl), Dimension(:) :: list
    Integer :: i
    
    maxr = list(1)
    Do i = 2,Size(list)
      If (list(i) > maxr) Then
        maxr = list(i)
      End If
    End Do
  End Function Maxr

  !----------------------------------------------
  ! Returns the index of the first instance of 
  ! n in the array, 0 if not found
  !----------------------------------------------
  Integer Function utils_inarrayi(n,list)

    Integer, Intent(In) :: n
    Integer, Dimension(:), Intent(In) :: list
    Integer :: i

    Do i = 1,Size(list)
      If (list(i) == n) Then
        utils_inarrayi = i
        Return
      End If
    End Do

    utils_inarrayi = 0

  End Function utils_inarrayi

  !-------------------------------------------------------------------
  ! This routine multiplies two arrays as matrices and returns the
  ! product in "prod"
  !-------------------------------------------------------------------
  Subroutine multarray(arr1, arr2, prod)
    Real(kind=RDbl), Dimension(:,:) :: arr1
    Real(kind=RDbl), Dimension(:,:) :: arr2
    Real(kind=RDbl), Dimension(:,:) :: prod

    Integer    :: dim1_1, dim1_2, dim2_1, dim2_2, dimp_1, dimp_2
    Integer    :: row, col, j
    Real(kind=RDbl) :: sum
    
    !** Make sure the dimensions of the array are consistent with matrix
    !** multiplication
    dim1_1 = Size(arr1, 1)   ! No. of rows
    dim1_2 = Size(arr1, 2)   ! No. of columns
    dim2_1 = Size(arr2, 1)
    dim2_2 = Size(arr2, 2)
    dimp_1 = Size(prod, 1)
    dimp_2 = Size(prod, 2)

    If (dim1_2 /= dim1_1) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    If (dim1_1 /= dimp_1 .Or. dim2_2 /= dimp_2) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

    Do row=1, dim1_1
      Do col=1, dim2_2
        sum = 0.0_RDbl
        Do j=1, dim1_2
          sum = sum + arr1(row, j)*arr2(j, col)
        End Do
        prod(row, col) = sum
      End Do
    End Do
  End Subroutine multarray

  !--------------------------------------------------------------
  ! This routine multiplies an array("arr") and a vector("vec")
  ! and returns the result in the array "prod"
  !-------------------------------------------------------------
  Subroutine multarrvec(arr, vec, prod)
    Real(kind=RDbl), Dimension(:,:), Intent(in) :: arr
    Real(kind=RDbl), Dimension(:), Intent(in)   :: vec
    Real(kind=RDbl), Dimension(:), Intent(out)  :: prod

    Integer      :: i, j, arrdim1, arrdim2, vecdim, proddim
    Real(kind=RDbl)   :: vecel
    
    arrdim1 = Size(arr, 1)
    arrdim2 = Size(arr, 2)
    vecdim  = Size(vec)
    proddim = Size(prod)

    If (vecdim /= proddim .Or. arrdim2 /= vecdim) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " Array dimensions don't match"
      Stop
    End If
    
    prod = 0.0_RDbl
    Do j = 1, arrdim2
      vecel = vec(j)
      Do i=1, arrdim1
        prod(i) = prod(i) + arr(i, j)*vecel
      End Do
    End Do
  End Subroutine multarrvec

  !---------------------------------------------------------------
  ! Multiplies two vectors "vec1" and "vec2" and returns the "prod"
  !---------------------------------------------------------------
  Subroutine multvecvec(vec1, vec2, prod)
    Real(kind=RDbl), Dimension(:), Intent(in) :: vec1, vec2
    Real(kind=RDbl), Intent(out)              :: prod

    Integer      :: i

    prod = 0.0_RDbl
    Do i=1, Size(vec1, 1)
      prod = prod + vec1(i)*vec2(i)
    End Do
  End Subroutine multvecvec

  !---------------------------------------------------------------
  ! Displays the contents of the array "arr" and outputs to the
  ! optional unitno "optunitno"
  !---------------------------------------------------------------
  Subroutine disparray(arr, label, fmt, optunitno)
    Real(kind=RDbl), Dimension(:, :), Intent(in) :: arr
    Character(*), Intent(in)        :: label
    Character(*), Intent(in)        :: fmt
    Integer, Optional, Intent(in)   :: optunitno

    Integer     :: i, j, unitno, dim1, dim2
    Character(len=strLen)     :: strfmt
    
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6      ! Standard output
    End If

    strfmt = "("//Trim(fmt)//")"
    dim1 = Size(arr, 1)
    dim2 = Size(arr, 2)
    Write(unitno, '(a)') label
    Do i=1, dim1
      Do j=1, dim2
        Write(unitno, strfmt, Advance='No') arr(i,j)
      End Do
      Write(unitno, *)
    End Do
  End Subroutine disparray

  !---------------------------------------------------------------
  ! Displays the contents of the array "arr" and outputs to the
  ! optional unitno "optunitno"
  !---------------------------------------------------------------
  Subroutine dispvec(vec, fmt, optunitno)
    Real(kind=RDbl), Dimension(:), Intent(in) :: vec
    Character(*), Intent(in)        :: fmt
    Integer, Optional, Intent(in)   :: optunitno

    Integer     :: i, unitno, dim1
    Character(len=strLen)     :: strfmt
    
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6      ! Standard output
    End If

    strfmt = "("//Trim(fmt)//")"
    dim1 = Size(vec, 1)
    Do i=1, dim1
      Write(unitno, strfmt, Advance='No') vec(i)
    End Do
  End Subroutine dispvec

  !-----------------------------------------------------------------
  ! This routine gets the norm of a vector "arr"
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function arrnorm(arr)
    Real(kind=RDbl), Dimension(:), Intent(in) :: arr
    
    Integer       :: i
    
    arrnorm = 0.0_RDbl
    Do i=1, Size(arr)
      arrnorm = arrnorm + arr(i)*arr(i)
    End Do

    If (arrnorm > 0.0_RDbl) Then
      arrnorm = Sqrt(arrnorm)
    Else
      arrnorm = 0.0_Rdbl
    End If
    Return
  End Function arrnorm

!-----------------------------------------------------------
! Some other random stuff  
!-----------------------------------------------------------
  !-------------------------------------------------------------
  ! Calculates the factorial of a number
  !-------------------------------------------------------------
  Integer Function factorial(n)
    Integer, Intent(in) :: n
    
    Integer    :: i

    factorial = 1
    Do i=1, n
      factorial = factorial*i
    End Do
  End Function factorial

  !------------------------------------------------------------------
  ! Gets the no. of combinations of "n" objects taken "c" at a time
  !------------------------------------------------------------------
  Integer Function comb(n, c)
    Integer, Intent(in) :: n, c

    comb = factorial(n)/(factorial(n-c)*factorial(c))
  End Function comb

  !---------------------------------------------------------------
  ! This function returns the machine precision of the variable
  ! "eps".  This algorithm has been taken from "Numerical Methods"
  ! by Kahaner, Moler and Nash.
  !---------------------------------------------------------------
  Real(kind=RDbl) Function getMachPrec(eps)
    Real(kind=RDbl), Intent(inout) :: eps
    
    Real(kind=RDbl)  :: oldeps, eps1, lasteps
        
    !** Save the value of eps
    oldeps = eps

    !** Get the machine precision.
    eps = 1.0
    eps1 = eps + 1.0
    Do 
      If (eps1 <= 1.0) Exit
      lasteps = eps
      eps = eps/2.0
      eps1 = eps + 1.0
    Enddo

    eps = oldeps
    getMachPrec =  lasteps
    Return
  End Function getMachPrec

  !---------------------------------------------------------------
  ! This function returns the machine precision of the variable
  ! "eps".  This algorithm has been taken from "Numerical Methods"
  ! by Kahaner, Moler and Nash.
  !---------------------------------------------------------------
  Real(kind=RDbl) Function getMachRange(eps)
    Real(kind=RDbl), Intent(inout) :: eps
    
    Real(kind=RDbl)  :: oldeps, lasteps
        
    !** Save the value of eps
    oldeps = eps

    !** Get the machine precision.
    eps = 1.0
    Do 
      If (eps <= 0.0) Exit
      lasteps = eps
      eps = eps/2.0
    Enddo

    !** Restore the initial value
    eps = oldeps
    getMachRange =  lasteps
    Return
  End Function getMachRange

  !------------------------------------------------------------------
  ! This routine returns an angle "theta" in the correct quadrant given
  ! values of sin(theta), cos(theta)
  !------------------------------------------------------------------
  Subroutine getinvangle(sintheta, costheta, theta, optrange)
    Real(kind=RDbl), Intent(in) :: sintheta, costheta
    Real(kind=RDbl), Intent(out)  :: theta
    Logical, Optional, Intent(in) :: optrange

    Logical    :: zero_twopi

    If (Present(optrange)) Then
      zero_twopi = optrange
    Else
      zero_twopi = .True.
    End If

    If ((costheta*costheta) > one*1.00000001) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write (*,*) "wrong sin cos values passed here"
      Write (*,*) "sin :", sintheta,  "cos :", costheta
      stop

    Endif

    If (dbgflag) Write(*,'(1x,2a,i4,e12.3)') __FILE__," : ",__LINE__, costheta
    theta = Acos(costheta)      ! theta is between 0 and pi
    If (sintheta < 0.0_RDbl) Then
      If (zero_twopi) Then
        theta = twopi - theta
      Else
        theta = -theta
      End If
    End If
  End Subroutine getinvangle

  !-----------------------------------------------------
  ! Swaps two (double) real numbers
  !-----------------------------------------------------
  Subroutine realswap(real1,real2)
    Real(Kind=RDbl), Intent(InOut)  :: real1, real2
    Real(Kind=RDbl)                 :: temp

    temp = real1
    real1 = real2
    real2 = temp

  End Subroutine realswap

  !-----------------------------------------------------
  ! Swaps two integers
  !-----------------------------------------------------
  Subroutine swap_integers(int1,int2)
    Integer, Intent(InOut)  :: int1, int2
    Integer                 :: temp

    temp = int1
    int1 = int2
    int2 = temp

  End Subroutine swap_integers

  !-----------------------------------------------------
  ! Swaps two subset arrays
  !-----------------------------------------------------
  Subroutine swap_subsets(subset1,subset2)
    Integer, Dimension(3), Intent(InOut)  :: subset1, subset2
    Integer, Dimension(3)                 :: temp

    temp = subset1
    subset1 = subset2
    subset2 = temp

  End Subroutine swap_subsets

  !----------------------------------------------------------------------------
  ! Sums the logical array with .AND. and returns the result
  !----------------------------------------------------------------------------
  Logical Function sumlogical(values)
    Logical, Dimension(:), Intent(In) :: values
    Integer :: i

    sumlogical = values(1)
    Do i = 2, Size(values,1)
      sumlogical = (sumlogical .And. values(i))
    End Do
  End Function sumlogical

  !----------------------------------------------------------------------------
  ! Compares two logical variables, if they are both of the same value, 
  ! then it returns True, otherwise it's False.  I don't understand why
  ! .Eqv. isn't work for me, so I wrote this.
  !----------------------------------------------------------------------------
  Logical Function chklogical(logical1,logical2)
    Logical, Intent(In) :: logical1,logical2

    chklogical = .False.
    If (logical1 .Or. logical2) Then  !** at least one is True
      !** See if both are True
      If (logical1 .And. logical2) chklogical = .True.  
    Else
      !** Both could still be False
      If ((.Not. logical1) .And. (.Not. logical2)) chklogical = .True.  
    End If

  End Function chklogical

  !-----------------------------------------------------------------------
  ! Check a logical, if it's .False. then stop with an error message
  ! Requires:  flag -- flag to check
  !            filename -- filename to display (use __FILE__ as input)
  !            lineno -- line number to display (use __LINE__ as input)
  !            message -- optional message
  !            integers -- optional integers to display at end of message
  !-----------------------------------------------------------------------
  Subroutine checkandstop(flag,filename,lineno,message,integers)
    Logical, Intent(In)                         :: flag
    Character(*), Intent(In)                    :: filename
    Integer, Intent(In)                         :: lineno
    Character(*), Intent(In), Optional          :: message
    Integer, Dimension(:), Intent(In), Optional :: integers

    Integer                   :: i
    Character(len=strLen)     :: string
    Character(len=xlstrLen)   :: output

    !** Skip out if there's no problem
    If (flag) Return

    !** Otherwise, write error message
    string = int2str(lineno)
    Write(0,'(1x,4a)') 'ERROR: line number ',Trim(string),' in file ', &
        Trim(filename)
    If (Present(message)) Then
      Write(output,'(1x,a)') Trim(message)

      !** Add integers to the end if they're present
      If (Present(integers)) Then
        Do i = 1,Size(integers)
          string = int2str(integers(i))
          Write(output,'(a,1x,a)') Trim(output),Trim(string)
        End Do
      End If

      !** Write to std error
      Write(0,'(1x,a)') Trim(output)
    End If

    Stop

  End Subroutine checkandstop
  
  !--------------------------------------------------------------
  ! Displays pointer allocation errors, and stops the simulation
  !--------------------------------------------------------------
  Subroutine allocErrDisplay(filename,lineno,opt_variableName)
    Character(*), Intent(In)           :: filename
    Integer, Intent(In)                :: lineno
    Character(*), Intent(In), Optional :: opt_variableName

    Character(len=strlen) :: varname,string
    
    If (Present(opt_variableName)) Then
      varname = opt_variableName
    Else
      varname = ""
    End If
    
    string = int2str(lineno)
    Write(0,'(1x,6a)') "Error while allocating variable ", &
        Trim(varname)," in '",Trim(filename),"' at line number ",Trim(string)
    Stop

  End Subroutine allocErrDisplay

  !--------------------------------------------------------------
  ! Displays pointer allocation errors, and stops the simulation
  !--------------------------------------------------------------
  Subroutine deallocErrDisplay(filename,lineno,opt_variableName)
    Character(*), Intent(In)           :: filename
    Integer, Intent(In)                :: lineno
    Character(*), Intent(In), Optional :: opt_variableName

    Character(len=strlen) :: varname,string
    
    If (Present(opt_variableName)) Then
      varname = opt_variableName
    Else
      varname = ""
    End If
    
    string = int2str(lineno)
    Write(0,'(1x,6a)') "Error while DEallocating variable ", &
        Trim(varname)," in '",Trim(filename),"' at line number ",Trim(string)
    Write(0,'(1x,a)') "This may occur when the variable was never allocated"
    Stop

  End Subroutine deallocErrDisplay
  
  !----------------------------------------------------------------------------
  !  This subroutine returns an approximation to the complementary error 
  !  Function
  !  Reference:  Abramowitz and Stegun, Handbook of Mathematical Functions,
  !              National Bureau of Standards, formula 7.1.26        
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function erfc_regular(x)
    Real(Kind=RDbl), Intent(In) :: x
    Real(Kind=RDbl) :: t,xsq,tp
    Real(Kind=RDbl), Parameter :: a1 = 0.254829592, a2 = -0.284496736, &
        a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p  =  0.3275911

    t  = 1.0/(1.0 + p*x)
    xsq = x*x
    tp = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))
    erfc_regular = tp*Exp(-xsq)
    Return
  End Function erfc_regular

  !----------------------------------------------------------------------------
  !  This subroutine returns an approximation to the complementary error 
  !  Function
  !  Reference:  Abramowitz and Stegun, Handbook of Mathematical Functions,
  !              National Bureau of Standards, formula 7.1.26        
! in addiiton returns e(-x^2)
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function erfc_and_exp(x,opt_ex2)
    Real(Kind=RDbl), Intent(In) :: x
    Real(Kind=RDbl), Intent(Out) :: opt_ex2
    Real(Kind=RDbl) :: t,xsq,tp
    Real(Kind=RDbl), Parameter :: a1 = 0.254829592, a2 = -0.284496736, &
        a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p  =  0.3275911

    t  = 1.0/(1.0 + p*x)
    xsq = x*x
    tp = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))
    opt_ex2=Exp(-xsq)
    erfc_and_exp = tp*opt_ex2
    Return
  End Function erfc_and_exp


  !--------------------------------------------------------------------------
  ! Counts the no of lines in a text file ("unitno" is that of the text file)
  ! Copies those lines into "stringArray".  This is useful while doing 
  ! post-code things control files can be copied easily and written to output 
  ! files.  no of lines in the file is passed back if "nlines" is present.
  ! Requires:  unitno -- unit number to read from
  !            stringArray -- array of strings to fill, this is a POINTER
  !            no_of_lines -- optional returned number of lines found
  !--------------------------------------------------------------------------
  Subroutine FileToStringArr(unitno,stringArray,no_of_lines)
    Integer, Intent(In)                            :: unitno
    Character(len=2*strLen), Dimension(:), Pointer :: stringArray 
    Integer, Intent(Out), Optional                 :: no_of_lines

    Integer                  :: i,error,nlines
    Character (len=2*strLen) :: inputbuffer

    !** Rewind and count the no-of-lines from beginning
    Rewind(unitno)
    nlines = 0
    Do 
      Read(unitno,'(a)',IOSTAT=error) inputbuffer
      If (error /= 0) Then
        Exit
      Endif
      nlines = nlines + 1
    End Do
    
    !** Size the string array
    Allocate(stringArray(nlines),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    !** Read values into stringArray
    Rewind(unitno)
    Do i=1,nlines
      Read(unitno,'(a)') inputbuffer
      stringArray(i)=inputbuffer
    End Do
    
    Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
    Write(*,*) stringArray

    If (present(no_of_lines)) no_of_lines=nlines

  End Subroutine FileToStringArr

  !--------------------------------------------------------------------------
  ! Make a system call.  This is a separate routine because the call may
  ! depend on the compiler.  If the optional error message is supplied, the
  ! routine will stop execution after an error and print this message.
  ! Returns True if the if there were no errors, False otherwise.
  ! Requires:  command -- the UNIX command line entry
  !            errormsg -- optional error message to print after failure
  !--------------------------------------------------------------------------
  Logical Function syscall(command,errormsg)
    Character(*), Intent(In)            :: command
    Character(*), Intent(In), Optional  :: errormsg

    Integer              :: status,errno
    Logical              :: error
#ifndef NAGCOMPILER
    Integer              :: system
#endif

    !** Execute the system call depending on the compiler type
    syscall = .False.
#ifndef NAGCOMPILER
    errno = system(command)
    If (errno == 0) syscall = .True.
#endif
#ifdef NAGCOMPILER
    Call system(command,status,errno)
    If (errno == 0) syscall = .True.
#endif

    !** Return successful if there was no error
    If (syscall) Return

    !** Write error message to stderror if error message was passed
    If (Present(errormsg)) Then
      Write(0,'(1x,a)') 'ERROR during system call.  The system command was:'
      Write(0,'(2x,3a)') '"',Trim(command),'"'
      Write(0,'(2x,a)') 'message from calling routine:'
      Write(0,'(2x,a)') Trim(errormsg)
      Stop
    End If

  End Function syscall

End Module utils

