!---------------------------------------------------------------
! This module defines the data structure of a 3-D coordinate and 
! the operations that accompany it. All operations are defined on 
! vectors with components of type real.  Vector of type integer
! is also provided including assignment operations that convert
! back and forth between the real and integer types
!----------------------------------------------------------------
Module vector

  Use defaults, Only: RDbl, strLen, Pi, zero, zeroTolerance, lstrLen, degTorad
  Use utils, Only: real2str

  Implicit None
  Save

  Private
  Public :: VecType, Assignment(=), Operator(*), Operator(+), &
      Operator(-), Operator(/), vector_display,vector_filedisplay, &
      mag , unitvec, IntVecType, vector_tripledotprod, vector_getnormsq, &
      vector_getcomp, vector_iscollinear, vector_crossprod, &
      vector_getnorm, vector_getdistsq, vector_sumVecs, vector_bondangle, &
      vector_angle, vector_getcom, vector_getunitvec, vector_getdist, &
      vector_getplanenorm, swap, vector_ptonplane, vector_zerovec

  ! When allocatable arrays are allowed in
  ! structures we can make this module more
  ! general.
  Integer, Parameter  :: vec_size = 3

  Type IntVecType
    Integer    :: comp(vec_size)
  End Type IntVecType

  Type VecType
    Real(kind=RDbl)    :: comp(vec_size)
  End Type VecType

  Interface Assignment(=)
    Module Procedure vector_vectype_eq_intvectype
    Module Procedure vector_vec_eq_scalar
    Module Procedure vector_intvectype_eq_scalar
    Module Procedure vector_vectype_eq_rarray
    Module Procedure vector_vectype_eq_intarray
    Module Procedure vector_rarray_eq_vectype
  End Interface

  Interface Operator(+)
    Module Procedure vector_add
    Module Procedure vector_vec_plus_array
  End Interface

  Interface Operator(-)
    Module Procedure vector_subtract
  End Interface

  Interface Operator(*)
    Module Procedure vector_dotprod
    Module Procedure vector_vecxscalar
    Module Procedure vector_scalarxvec
    Module Procedure vector_vecxinteger
    Module Procedure vector_integerxvec
    Module Procedure vector_dotprod_AryxVec
  End Interface

  Interface Operator(/)
    Module Procedure vector_intdivide
    Module Procedure vector_realdivide
  End Interface

  !** Vector display function-interface
  Interface vector_display
    Module Procedure vector_intstrdisplay
    Module Procedure vector_realstrdisplay
  End Interface

  Interface unitvec
    Module Procedure vector_getunitvec
  End Interface

  Interface mag
    Module Procedure vector_getnorm
  End Interface

  Interface isinplane
    Module Procedure vector_isin3dplane
  End Interface

  Interface getplanedist
    Module Procedure vector_get3dplanedist
    Module Procedure vector_getoutofplanedist
  End Interface

  Interface swap
    Module Procedure vector_swap
  End Interface

Contains
  !------------------------------------------------------------
  ! Assign an array to a real vector type
  !------------------------------------------------------------
  Subroutine vector_vectype_eq_rarray(vec1, array1)
    Type(VecType), Intent(out) :: vec1
    Real(kind=RDbl), Dimension(:), Intent(in) :: array1
    Integer   :: arraysize, i
    
    arraysize = Size(array1, 1)
    If (arraysize > vec_size) Then
      Write(0,'(a,i4)') 'Array is too big.  Vector size is fixed at :',vec_size
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    Do i=1, arraysize
      vec1%comp(i) = array1(i)
    Enddo
    Return
  End Subroutine vector_vectype_eq_rarray

  !------------------------------------------------------------
  ! Assign an integer array to a integer vector type
  !------------------------------------------------------------
  Subroutine vector_vectype_eq_intarray(vec1, array1)
    Type(IntVecType), Intent(out) :: vec1
    Integer, Dimension(:), Intent(in) :: array1
    Integer   :: arraysize, i
    
    arraysize = Size(array1, 1)
    If (arraysize > vec_size) Then
      Write(0,'(a,i4)') 'Array is too big.  Vector size is fixed at :', vec_size
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    Do i=1, arraysize
      vec1%comp(i) = array1(i)
    Enddo
    Return
  End Subroutine vector_vectype_eq_intarray

  !----------------------------------------------
  ! Assign a real vector type to an array 
  !----------------------------------------------
  Subroutine vector_rarray_eq_vectype(array1, vec1)
    Real(kind=RDbl), Dimension(:), Intent(out) :: array1
    Type(VecType), Intent(in) :: vec1

    Integer   :: arraysize, i
    
    arraysize = Size(array1, 1)
    If (arraysize < vec_size) Then
      Write(0,'(a,i4)') 'Array is too small.  Vector size is :', vec_size
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If
    Do i=1, vec_size
      array1(i) = vec1%comp(i)
    Enddo
    Return
  End Subroutine vector_rarray_eq_vectype
  

  !----------------------------------------------------------
  ! Assign an integer vector type to a real vector type
  !----------------------------------------------------------
  Subroutine vector_vectype_eq_intvectype(vec1, intvec2)
    Type(VecType), Intent(out)   :: vec1
    Type(IntVecType), Intent(in) :: intvec2
    Integer           :: i

    Do i = 1, vec_size
      vec1%comp(i) = intvec2%comp(i)
    End Do
    Return
  End Subroutine vector_vectype_eq_intvectype


  !-------------------------------------------
  ! Initialize a vector's components to "num1"
  !-------------------------------------------
  Subroutine vector_vec_eq_scalar(vec1, num1)
    Type(VecType), Intent(out) :: vec1
    Real(kind=RDbl), Intent(in) :: num1
    Integer   :: i
    
    Do i=1, vec_size
      vec1%comp(i) = num1
    End Do
    Return
  End Subroutine vector_vec_eq_scalar

  !-------------------------------------------
  ! Initialize a vector's components to "num1"
  !-------------------------------------------
  Subroutine vector_intvectype_eq_scalar(intvec1, num1)
    Type(IntVecType), Intent(out) :: intvec1
    Integer, Intent(in) :: num1
    Integer   :: i, nsize
    
    nsize = Size(intvec1%comp, 1)
    Do i=1, nsize
      intvec1%comp(i) = num1
    End Do
    Return
  End Subroutine vector_intvectype_eq_scalar

  !--------------------------------------------------
  ! returns a zero vector
  !--------------------------------------------------
  Type(VecType) Function vector_zerovec()
      vector_zerovec = zero 
    Return
  End Function vector_zerovec

  !--------------------------------------------------
  ! Add two vector structures
  !--------------------------------------------------
  Type(VecType) Function vector_add(vec1, vec2)
    Type(VecType), Intent(in)  :: vec1, vec2
    Integer               :: i
    
    Do i=1, vec_size
      vector_add%comp(i) = vec1%comp(i) + vec2%comp(i)
    End Do
    Return
  End Function vector_add

  !----------------------------------------------------------------------------
  ! Add an array to a vector
  !----------------------------------------------------------------------------
  Type(VecType) Function vector_vec_plus_array(vec1,arr1)
    Type(VecType), Intent(in) :: vec1
    Real(Kind=RDbl), Dimension(:), Intent(in) :: arr1
    Integer :: i

    Do i = 1, vec_size
      vector_vec_plus_array%comp(i) = vec1%comp(i) + arr1(i)
    End Do

  End Function vector_vec_plus_array

  !--------------------------------------------------
  ! Subtract two vector structures
  !--------------------------------------------------
  Type(VecType) Function vector_subtract(vec1, vec2)
    Type(VecType), Intent(in)  :: vec1, vec2
    Integer               :: i
    
    Do i=1, vec_size
      vector_subtract%comp(i) = vec1%comp(i) - vec2%comp(i)
    Enddo
    Return
  End Function vector_subtract

  !--------------------------------------------------
  ! Get the product of a scalar with a vector
  !--------------------------------------------------
  Type(VecType) Function vector_vecxscalar(vec1, num1)
    Type(VecType), Intent(in)  :: vec1
    Real(kind=RDbl), Intent(in)  :: num1
    Integer      :: i
    
    Do i=1, vec_size
      vector_vecxscalar%comp(i) = vec1%comp(i)*num1
    End Do
    Return
  End Function vector_vecxscalar

  !--------------------------------------------------
  ! Get the product of a scalar with a vector
  !--------------------------------------------------
  Type(VecType) Function vector_scalarxvec(num1,vec1)
    Real(kind=RDbl), Intent(in)  :: num1
    Type(VecType), Intent(in)  :: vec1
    Integer      :: i
    
    Do i=1, vec_size
      vector_scalarxvec%comp(i) = vec1%comp(i)*num1
    End Do
    Return
  End Function vector_scalarxvec

  !--------------------------------------------------
  ! Get the product of an integer with a vector
  !--------------------------------------------------
  Type(VecType) Function vector_vecxinteger(vec1, num1)
    Type(VecType), Intent(In)  :: vec1
    Integer, Intent(In)        :: num1

    Integer      :: i
    
    Do i = 1,vec_size
      vector_vecxinteger%comp(i) = vec1%comp(i)*num1
    End Do

  End Function vector_vecxinteger

  !--------------------------------------------------
  ! Get the product of a integer with a vector
  !--------------------------------------------------
  Type(VecType) Function vector_integerxvec(num1,vec1)
    Integer, Intent(In)        :: num1
    Type(VecType), Intent(in)  :: vec1

    Integer      :: i
    
    Do i = 1,vec_size
      vector_integerxvec%comp(i) = vec1%comp(i)*num1
    End Do

  End Function vector_integerxvec

  !--------------------------------------------------
  ! Get the dot product of two vectors
  !--------------------------------------------------
  Real(kind=RDbl) Function vector_dotprod(vec1, vec2)
    Type(VecType), Intent(in)  :: vec1, vec2
    Real(kind=RDbl)  :: dotprod
    Integer               :: i

    dotprod = 0.0_RDbl

    Do i=1, vec_size
      dotprod = dotprod + vec1%comp(i)*vec2%comp(i)
    Enddo
    vector_dotprod = dotprod
    Return
  End Function vector_dotprod

  !----------------------------------------------------------------------------
  ! Returns the dot product of a Fortran array-type vector and VecType vector
  !----------------------------------------------------------------------------
  Real(Kind=RDbl) Function vector_dotprod_AryxVec(arr,vec)
    Type(VecType), Intent(In) :: vec
    Real(Kind=RDbl), Dimension(:), Intent(In) :: arr
    
    !** Verify they are the same length. If not, exit
    If (Size(arr,1) /= vec_size) Then
      Write(0,'(3(a,i2),a)') __FILE__//":", __LINE__,&
          " Array length of ", Size(arr,1), " is longer than the default &
          &vec length ", vec_size, " as is required by dotprod."
      Stop
    End If

    !** Get the answer
    vector_dotprod_AryxVec = arr(1)*vec%comp(1) + arr(2)*vec%comp(2) + &
        arr(3)*vec%comp(3)

  End Function vector_dotprod_AryxVec

  !-----------------------------------------------------
  ! Swaps two vectors
  !-----------------------------------------------------
  Subroutine vector_swap(vec1,vec2)
    Type(VecType), Intent(InOut)  :: vec1,vec2

    Type(VecType)                 :: temp

    temp = vec1
    vec1 = vec2
    vec2 = temp

  End Subroutine vector_swap

  !----------------------------------------------------------------------
  ! Get the angle between two bond vectors.  Make sure the vectors are 
  ! defined correctly or this routine will, of course, return junk.
  ! Requires:  vec12 = (atom1 - atom2) vector
  !            vec32 = (atom3 - atom2) vector
  !            adotb = dot product between vec12 and vec32
  !            costheta = cosine of the angle
  !            theta = the angle itself (in radians)
  !----------------------------------------------------------------------
  Subroutine vector_bondangle(vec12,vec32,adotb,costheta,theta)
    Type(VecType), Intent(In)    :: vec12,vec32
    Real(kind=RDbl), Intent(Out) :: theta,costheta,adotb

    Real(kind=RDbl)              :: norm1,norm2

    norm1 = vector_getnorm(vec12)
    norm2 = vector_getnorm(vec32)

    adotb = vector_dotprod(vec12,vec32)
    costheta = adotb/(norm1*norm2)
    theta = dacos(costheta)

  End Subroutine vector_bondangle

  !-----------------------------------------------------------------
  ! Get the angle between three position vectors in degrees
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function vector_angle(vec1,vec2,vec3)
    Type(VecType), Intent(in)  :: vec1,vec2,vec3
    Real(kind=RDbl)            :: dotprod
    Type(VecType)              :: vec21,vec23

    vec21 = vector_getunitvec(vec2 - vec1)
    vec23 = vector_getunitvec(vec2 - vec3)

    dotprod = vector_dotprod(vec21,vec23)
    vector_angle = Acos(dotprod)/degTorad

  End Function vector_angle
 
  !--------------------------------------------------
  ! Get the cross product of two vectors
  !--------------------------------------------------
  Type(VecType) Function vector_crossprod(vec1, vec2)
    Type(VecType), Intent(in)  :: vec1, vec2
    Type(VecType)    :: cprod
    
    cprod%comp(1) =   vec1%comp(2)*vec2%comp(3) - vec1%comp(3)*vec2%comp(2)
    cprod%comp(2) = -(vec1%comp(1)*vec2%comp(3) - vec1%comp(3)*vec2%comp(1))
    cprod%comp(3) =   vec1%comp(1)*vec2%comp(2) - vec1%comp(2)*vec2%comp(1)
    vector_crossprod = cprod
  End Function vector_crossprod

  !--------------------------------------------------
  ! Get the scalar triple product of 3 vectors
  !--------------------------------------------------
  Real(kind=RDbl) Function vector_tripledotprod(v1, v2, v3)
    Type(VecType), Intent(in)  :: v1, v2, v3
    Type(VecType)    :: temp

    ! Triple product = (v1 X v2).v3 
    
    ! Take the cross product of v1, v2
    temp = vector_crossprod(v1, v2)
    ! Finally take the dot product with v3
    vector_tripledotprod = temp*v3
    Return
  End Function vector_tripledotprod

  !--------------------------------------------------
  ! Divide a vector by an int
  !--------------------------------------------------
  Type(VecType) Function vector_intdivide(vec1, int1)
    Type(VecType), Intent(in)  :: vec1
    Integer, Intent(in)        :: int1
    Integer   :: i
    
    Do i=1, vec_size
      vector_intdivide%comp(i) = vec1%comp(i)/int1
    End Do
    Return
  End Function vector_intdivide

  !--------------------------------------------------
  ! Divide a vector by a real
  !--------------------------------------------------
  Type(VecType) Function vector_realdivide(vec1, real1)
    Type(VecType), Intent(in)  :: vec1
    Real(kind=RDbl), Intent(in) :: real1
    Integer      :: i
    
    Do i=1, vec_size
      vector_realdivide%comp(i) = vec1%comp(i)/real1
    Enddo
    Return
  End Function vector_realdivide

  !----------------------------------------------------------------------
  ! Get the component number "cmp" of the vector "vec1"
  !----------------------------------------------------------------------
  Real(kind=RDbl) Function vector_getcomp(vec1, cmp)
    Type(VecType), Intent(in)   :: vec1
    Integer, Intent(in)         :: cmp
    
    vector_getcomp = vec1%comp(cmp)
  End Function vector_getcomp

  !---------------------------------------------------------------------
  ! Check if the vectors vec1, vec2 are collinear
  !---------------------------------------------------------------------
  Logical Function vector_iscollinear(vec1, vec2)
    Type(VecType), Intent(In) :: vec1, vec2
    
    Integer         :: i
    Real(kind=RDbl) :: factor
    Type(VecType)   :: vec3

    !** Get the multiplication factor which when multiplied by vec1 
    !** would result in vec2 if the vectors are collinear
    factor = 0.0_RDbl
    Do i = 1,vec_size
      If (Abs(vec1%comp(i)) > zeroTolerance) Then
        factor = vec2%comp(i)/vec1%comp(i)
        Exit
      End If
    End Do
    If (factor == 0.0_Rdbl) Then
      !** vec1 is a null vector and hence vectors are collinear
      vector_iscollinear = .True.
      Return
    End If

    !** Now check if we can generate vec2
    vec3 = vec1*factor
    If (mag(vec2 - vec3) < zeroTolerance) Then
      vector_iscollinear = .True.
      Return
    End If
    vector_iscollinear = .False.

  End Function vector_iscollinear

  !---------------------------------------------------------------------
  ! Performs Gram-Schmidt orthogonalization on the array of 
  ! vectors "vecs" to generate "orthovecs".  See 
  ! "Linear Algebra and its Applications" by Strang, Pg. 172, 3rd Ed.
  !---------------------------------------------------------------------
  Subroutine vector_gsortho(vecs, orthovecs)
    Type(VecType), Dimension(:), Intent(in) :: vecs
    Type(VecType), Dimension(:), Intent(out):: orthovecs

    Integer         :: i, j
    Real(kind=RDbl) :: dotprod

    Do i=1, vec_size
      orthovecs(i) = vecs(i)
      Do j=1, i-1
        dotprod = orthovecs(i)*orthovecs(j)
        orthovecs(i) = orthovecs(i) - orthovecs(j)*dotprod
      End Do
      If (vector_getnormsq(orthovecs(i)) < 1.0e-8_RDbl) Then
        Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
            " Initial Vectors are not linearly independent "
        Stop
      Else
        orthovecs(i) = vector_getunitvec(orthovecs(i))        
      End If
    End Do
  End Subroutine vector_gsortho

  !----------------------------------------------------------
  ! Gets the unit vector in the direction of "vec"
  !----------------------------------------------------------
  Type(VecType) Function vector_getunitvec(vec)
    Type(VecType), Intent(In) :: vec

    Real(kind=RDbl) :: magnitude

    magnitude = vector_getnorm(vec)
    vector_getunitvec = vector_realdivide(vec,magnitude)

  End Function vector_getunitvec

  !-----------------------------------------------------
  ! Get the norm square of a vector
  !-----------------------------------------------------    
  Real(kind=RDbl) Function vector_getnormsq(vec1)
    Type(VecType), Intent(in)  :: vec1
    Integer             :: i

    vector_getnormsq = 0.0_RDbl
    
    Do i=1, vec_size
      vector_getnormsq = vector_getnormsq + vec1%comp(i)**2
    End Do
    Return
  End Function vector_getnormsq

  !-----------------------------------------------------
  ! Get the norm of a vector "vec1"
  !-----------------------------------------------------
  Real(kind=RDbl) Function vector_getnorm(vec1)
    Type(VecType), Intent(in) :: vec1

    vector_getnorm =  sqrt(vector_getnormsq(vec1))
  End Function vector_getnorm

  !--------------------------------------------------------------
  ! Get the distance between two coordinates "vec1" and "vec2"
  !--------------------------------------------------------------
  Real(kind=RDbl) Function vector_getdist(vec1, vec2)
    Type(VecType), Intent(IN) :: vec1, vec2
    Type(VecType)   :: temp
    
    temp = vec1 - vec2
    vector_getdist = vector_getnorm(temp)

  End Function vector_getdist

  !-----------------------------------------------------
  ! Get the squared distance between two coordinates
  !-----------------------------------------------------
  Real(kind=RDbl) Function vector_getdistsq(vec1, vec2)
    Type(VecType), Intent(IN) :: vec1, vec2
    Type(VecType)             :: temp
    
    temp = vec1 - vec2
    vector_getdistsq = vector_getnormsq(temp)

    Return
  End Function vector_getdistsq

  !----------------------------------------------------------------------------
  ! Returns the sum of all the vectors of the dimensioned vecs
  !----------------------------------------------------------------------------
  Type(VecType) Function vector_sumVecs(vecs)
    Type(VecType), Dimension(:), Intent(In) :: vecs
    Integer :: i

    Do i = 1, vec_size
      vector_sumVecs%comp(i) = Sum(vecs%comp(i))
    End Do

  End Function vector_sumVecs

  !-----------------------------------------------------
  ! Gets the center of mass of an array of vectors
  !-----------------------------------------------------
  Type(VecType) Function vector_getcom(veclist)
    Type(VecType), Dimension(:), Intent(in) :: veclist
    Integer     :: nsize, i
    
    vector_getcom = 0.0_RDbl
    nsize = Size(veclist,1)
    Do i=1, nsize
      vector_getcom = vector_getcom + veclist(i)
    End Do
    vector_getcom = vector_getcom/(nsize*1.0_RDbl)
  End Function vector_getcom

  !--------------------------------------------------------------
  ! This function checks if the points in the array "xyzcoords"
  ! lie in the plane normal to "normal" going throught the point
  ! "pt"
  !--------------------------------------------------------------
  Logical Function vector_isin3dplane(xyzcoords, pt, normal, opttol)
    Type(VecType), Intent(in)      :: xyzcoords
    Type(VecType), Intent(in)      :: pt, normal
    Real(kind=RDbl), Optional, Intent(in)  :: opttol

    Type(VecType)          :: inplanevec
    Real(kind=Rdbl)        :: dotprod, tol

    If (Present(opttol)) Then
      tol = opttol
    Else
      tol = 1.0e-5_RDbl
    End If

    inplanevec = xyzcoords - pt
      
    ! If inplanevec is actually in the plane then its
    ! dot product with the normal should be zero
    dotprod = inplanevec*normal
    
    If (Abs(dotprod) > tol) Then
      vector_isin3dplane = .False.
      Return
    End If
    vector_isin3dplane = .True.
    Return

  End Function vector_isin3dplane

  !-----------------------------------------------------------------
  ! This function gets the distance of the point "pt1" to the
  ! plane defined by the point "planept" and normal "planenormal"
  ! along the normal "normal1". 
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function vector_get3dplanedist(pt1,normal1,planept,planenormal)
    Type(VecType), Intent(in) :: pt1, normal1, planept, planenormal
    
    vector_get3dplanedist = planenormal*(planept - pt1)/(planenormal*normal1)
  End Function vector_get3dplanedist

  !-----------------------------------------------------------------
  ! This function returns the normal to a plane specified by three 
  ! points.  The returned vector is normalized
  ! normvec = (p2 - p1) .crossprod. (p3 - p1)  (and normalized)
  ! consequently, the normal is the one given by the right-hand rule
  !-----------------------------------------------------------------
  Type(VecType) Function vector_getplanenorm(pt)
    Type(VecType), Dimension(3), Intent(In) :: pt
    Type(VecType)                           :: a,b
    
    a = pt(2) - pt(1)
    b = pt(3) - pt(1)
    vector_getplanenorm = unitvec(vector_crossprod(a,b))

  End Function vector_getplanenorm

  !-----------------------------------------------------------------
  ! This Subroutine gets the distances from each point of a set of 
  ! points to a plane defined by a given three points.  The first 
  ! distance is returned as the value of the function and all the
  ! distances are returned in the array "distance"
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function vector_getoutofplanedist(planept,testpt,distance)
    Type(VecType), Dimension(3), Intent(In)               :: planept
    Type(VecType), Dimension(:), Intent(In)               :: testpt
    Real(kind=RDbl), Dimension(Size(testpt)), Intent(Out) :: distance

    Integer                               :: i
    Type(VecType)                         :: normal,vec

    normal = vector_getplanenorm(planept)

    Do i = 1,Size(testpt)
      vec = testpt(i) - planept(1)
      distance(i) = vec*normal
    End Do

    vector_getoutofplanedist = distance(1)

  End Function vector_getoutofplanedist

  !-----------------------------------------------------------------
  ! Given a point, find the point on the plane that connects the two
  ! points along the normal vector.  This could perhaps be viewed
  ! as a projection of a point onto a plane.
  ! Requires:  normal -- normal vector of plane, MUST BE UNIT VECTOR
  !            knownpt -- known point on the plane
  !            otherpt -- given point outside the plane
  !-----------------------------------------------------------------
  Type(VecType) Function vector_ptonplane(normal,knownpt,otherpt)
    Type(VecType), Intent(In)               :: normal,knownpt,otherpt

    Real(kind=RDbl)   :: dist
    Type(VecType)     :: vec1

    !** Make a vector point from knownpt -> otherpt
    vec1 = otherpt - knownpt

    !** Get the distance from the otherpt to the plane
    dist = vec1*normal

    !** Use distance and normal vector to find new plane point
    vector_ptonplane = otherpt - (dist*normal)
    
    If (.Not. vector_isin3dplane(vector_ptonplane, &
        knownpt,normal,1.0e-6_RDbl)) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
    End If

  End Function vector_ptonplane

  !----------------------------------------------------------
  ! Write out the components of the vector to unit "unitno"
  !----------------------------------------------------------
  Subroutine vector_filedisplay(vec1, unitno)

    Type(VecType), Intent(in) :: vec1
    Integer, Optional, Intent(in) :: unitno
    Integer        :: i, unit
    Character(len=strLen) :: strformat

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If

    If (vec_size < 10) Then
      Write(strformat, '(a,i1,2a)') "(1x, ",vec_size, "f8.3)"
    Else
      Write(strformat, '(a,i2,a)') "(1x, ",vec_size, "f8.3)"
    End If

    Write(unit,strformat) (vec1%comp(i), i=1, vec_size)
  End Subroutine vector_filedisplay

  !--------------------------------------------------------
  ! This function returns a string with the components of
  ! the vector formatted according to "fmt".
  ! Requires:  vec1 -- a vector type
  !            fmt -- optional format for number
  !--------------------------------------------------------
!!$  Character(len=lstrLen) Function vector_realstrdisplay(vec1, fmt)
  Function vector_realstrdisplay(vec1, fmt)
    Character(len=lstrLen)                :: vector_realstrdisplay
    Type(VecType), Intent(In)             :: vec1
    Character(*), Intent(In), Optional    :: fmt

    Character(len=strLen)      :: strformat,string1,string2,string3
    Integer                    :: i

    If (Present(fmt)) Then
      strformat = fmt
    
      If (vec_size < 10) Then
        Write(strformat, '(a,i1,2a)') "(1x, ",vec_size, Trim(fmt),")"
      Else
        Write(strformat, '(a,i2,2a)') "(1x, ",vec_size, Trim(fmt),")"
      End If
      Write(vector_realstrdisplay, strformat) (vec1%comp(i), i=1, vec_size)
    Else

      !** This piece of code was giving problems while ran on 
      !** contumacy With lahey compiler
 !SDEBUG
!!$      Write(vector_realstrdisplay,'(2(a,3x),a)') &
!!$          (Trim(real2str(vec1%comp(i),10)), i=1, vec_size)
!!$    End If
!!$!!    vector_realstrdisplay = Adjustl(vector_realstrdisplay)
 !SDEBUG

      !** Replacement for above lines
 !SDEBUG
      string1 = real2str(vec1%comp(1),10)
      string2 = real2str(vec1%comp(2),10)
      string3 = real2str(vec1%comp(3),10)
      Write(vector_realstrdisplay,'(a25,2x,a25,2x,a25)') &
          Trim(string1), Trim(string2), Trim(string3)
    End If
    vector_realstrdisplay = Adjustl(vector_realstrdisplay)
 !SDEBUG


  End Function vector_realstrdisplay

  !--------------------------------------------------------
  ! This function returns a string with the components of
  ! the vector formatted according to "fmt".
  !--------------------------------------------------------
!!$  Character(len=strLen) Function vector_intstrdisplay(vec1, fmt)
  Function vector_intstrdisplay(vec1, fmt)
    Character(len=strLen)      :: vector_intstrdisplay
    Type(IntVecType), Intent(in)  :: vec1
    Character(*), Intent(in)   :: fmt
    Character(len=strLen)      :: strformat
    Integer                    :: i

    strformat = fmt
    
    If (vec_size < 10) Then
      Write(strformat, '(a,i1,2a)') "(1x, ",vec_size, Trim(fmt),")"
    Else
      Write(strformat, '(a,i2,2a)') "(1x, ",vec_size, Trim(fmt),")"
    End If
    Write(vector_intstrdisplay, strformat) (vec1%comp(i), i=1, vec_size)
    vector_intstrdisplay = Adjustl(vector_intstrdisplay)

  End Function vector_intstrdisplay

End Module vector
