!------------------------------------------------------------------------------
! This module contains the shape and transformation information for a 
! fundamental cell.  The fundamental cell type contains shape and size
! information.  Angles are stored in radians and lengths in Angstroms.
!------------------------------------------------------------------------------

Module fundcell

  Use defaults, Only: RDbl, strLen, RSgl, degTorad, Pi, lstrLen, dbgflag
  Use utils, Only: toreal,minn,split,real2str
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_display, vector_crossprod, &
      IntVecType, vector_tripledotprod, vector_getnormsq, vector_getunitvec, &
      vector_getnorm
  Use matrix, Only: MatrixType, matrix_inverse, Operator(*), Assignment(=)

  Implicit None
  Save

  Private
  Public :: Fundamental_Cell, fundcell_slant, fundcell_getell, &
      fundcell_unslant, fundcell_getvolume, fundcell_init, &
      fundcell_simpleinit, fundcell_construct, fundcell_isortho, &
      fundcell_geteff, fundcell_idstring, fundcell_angles, &
      fundcell_minwidth, fundcell_width, fundcell_latticevec, &
      fundcell_snglMinImage, fundcell_multMinImage, fundcell_center, &
      fundcell_maptocell, fundcell_display, fundcell_cellstring, &
      fundcell_shift, fundcell_changevolume, fundcell_singleArrayMinImage, &
      fundcell_multminimage2, fundcell_listminimage, fundcell_oldbox

  Character(len=strLen), Parameter :: fundcell_idstring = "Fundcell_Info"

  Type Fundamental_Cell
    Logical                     :: orthoflag
    Real(kind=RDbl)             :: ell(3)                  ! Edge Lengths
    Real(kind=RDbl)             :: half_ell(3)             ! half Edge Lengths
    Real(kind=RDbl)             :: ell_inv(3)           ! inverse Edge Lengths
    Real(kind=RDbl)             :: anglebc,angleac,angleab ! cell angles
    Type(VecType)               :: origin     ! origin of slanted system
    Real(kind=RDbl)             :: eff(3)     ! effective size of unit cell
    Real(kind=RDbl)             :: volume     ! volume of cell
    Type(MatrixType)            :: slantmatrix,unslantmatrix ! xform matrices 
    Type(VecType), Dimension(3) :: lvec                      ! Lattice vectors
    Real(kind=RDbl)             :: width(3)     ! width of unit cell in a,b,c
    Real(kind=RDbl)             :: minwidth     ! minimum width of box
    Real(kind=RDbl)             :: mx, my, mz
  End Type Fundamental_Cell

! these will be stored during volume moves and copied back if move is not 
! accepted
Real(rdbl), Dimension(3) :: old_ell, old_eff

Contains
  !----------------------------------------------------------------------------
  ! Initialize the fundamental cell information
  ! It expects the following information and format in a file or passed:
  ! a, b, c                   # edge lengths (RDbl)
  ! anglebc, angleac, angleab # cell angles  (RDbl)  (alpha,beta,gamma)
  ! x0, y0, z0                # origin coordinates (RDbl)
  ! Requires: fcell -- fundamental cell type to be initialized
  !----------------------------------------------------------------------------
  Subroutine fundcell_init(fcell)
    Type(Fundamental_Cell), Intent(OUT)   :: fcell

    Integer                   :: i
    Real(kind=RDbl)           :: anglebc,angleac,angleab,err


    !** Process the information
    fcell%half_ell = 0.5_RDbl*fcell%ell
    fcell%ell_inv = 1.0_RDbl/fcell%ell

    !** just checking how far does it deviate from a 90 Deg fundcell
    err = Abs(fcell%anglebc - 90.000) + Abs(fcell%angleac - 90.000) + Abs(fcell%angleab - 90.000)
    If (err > 0.001) Then
      fcell%orthoflag = .False.
    Else
      fcell%orthoflag = .True.
    End If


    !** convert to radians
    fcell%anglebc = fcell%anglebc*degTorad
    fcell%angleac = fcell%angleac*degTorad
    fcell%angleab = fcell%angleab*degTorad

    !** Generate transformation matrices, lattice vectors and volume
    Call fundcell_initxform(fcell)
    Call fundcell_genlatticevects(fcell)
    Call fundcell_genwidths(fcell)
    Call fundcell_genvolume(fcell)

  End Subroutine fundcell_init


  !----------------------------------------------------------------------------
  ! scale by given factro
  !----------------------------------------------------------------------------
  Subroutine fundcell_changevolume(fcell,scale)
    Type(Fundamental_Cell), Intent(OUT)   :: fcell
    Real(kind=RDbl), Intent(In)           :: scale

    Integer                   :: i
    Real(kind=RDbl)           :: anglebc,angleac,angleab,err

    ! store for later use if required
old_ell=fcell%ell(1:3) 
old_eff=fcell%eff(1:3)

    fcell%ell(1:3) =       fcell%ell(1:3) * scale
    fcell%eff(1:3) =       fcell%eff(1:3) * scale

    !** Process the information
    fcell%half_ell = 0.5_RDbl*fcell%ell
    fcell%ell_inv = 1.0_RDbl/fcell%ell


    !** Generate transformation matrices, lattice vectors and volume
    !    Call fundcell_initxform(fcell)
    Call fundcell_genlatticevects(fcell)
    Call fundcell_genwidths(fcell)
    Call fundcell_genvolume(fcell)

  End Subroutine fundcell_changevolume


  !----------------------------------------------------------------------------
  ! back to the old box
  !----------------------------------------------------------------------------
  Subroutine fundcell_oldbox(fcell)
    Type(Fundamental_Cell), Intent(InOUT)   :: fcell

    Integer                   :: i
    Real(kind=RDbl)           :: anglebc,angleac,angleab,err

    ! store for later use if required
fcell%ell(1:3)=old_ell
fcell%eff(1:3)=old_eff

    !** Process the information
    fcell%half_ell = 0.5_RDbl*fcell%ell
    fcell%ell_inv = 1.0_RDbl/fcell%ell


    !** Generate transformation matrices, lattice vectors and volume
    !    Call fundcell_initxform(fcell)
    Call fundcell_genlatticevects(fcell)
    Call fundcell_genwidths(fcell)
    Call fundcell_genvolume(fcell)

  End Subroutine fundcell_oldbox


  !----------------------------------------------------------------------------
  ! Initialize the fundamental cell information from a string containing this
  ! information:  a_size b_size c_size alpha_angle beta_angle gamma_angle
  ! Requires:  fcell -- fundamental cell type to be initialized
  !            info -- string containing above information
  !----------------------------------------------------------------------------
  Subroutine fundcell_simpleinit(fcell,info)
    Type(Fundamental_Cell), Intent(OUT)   :: fcell
    Character(*), Intent(In)              :: info

    Integer                   :: nfields
    Character(len=lstrLen), Dimension(20) :: fields

    !** Split the information string and assign variables
    nfields = split(info,fields)
    fcell%ell(1) = toreal(fields(1),'a edge length in fundcell_simpleinit')
    fcell%ell(2) = toreal(fields(2),'b edge length in fundcell_simpleinit')
    fcell%ell(3) = toreal(fields(3),'c edge length in fundcell_simpleinit')

    !** Convert angles to radians
    fcell%anglebc = toreal(fields(4))*degTorad
    fcell%angleac = toreal(fields(5))*degTorad
    fcell%angleab = toreal(fields(6))*degTorad

    !** Set the orthorhombic field
    fcell%orthoflag = .False.

    !** Make assumptions about this stuff
    fcell%origin = (/0.0_RDbl,0.0_RDbl,0.0_RDbl/)
    fcell%eff = (/0.0_RDbl,0.0_RDbl,0.0_RDbl/)

    !** Process the information
    fcell%half_ell = 0.5_RDbl*fcell%ell
    fcell%ell_inv = 1.0_RDbl/fcell%ell

    !** Generate transformation matrices, lattice vectors and volume
    Call fundcell_initxform(fcell)
    Call fundcell_genlatticevects(fcell)
    Call fundcell_genwidths(fcell)
    Call fundcell_genvolume(fcell)

  End Subroutine fundcell_simpleinit

  !----------------------------------------------------------
  ! Calculates the volume of the fundamental cell
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------
  Subroutine fundcell_genvolume(fcell)
    Type(Fundamental_Cell), Intent(InOut)      :: fcell

    fcell%volume = vector_tripledotprod &
        (fcell%lvec(1), fcell%lvec(2), fcell%lvec(3))

  End Subroutine fundcell_genvolume

  !----------------------------------------------------------------------------
  ! Generate lattice vectors given the lattice constants and
  ! the angles between them.  The coordinate system is aligned
  ! such that "ell(1)" is along the x-axis, ell(2) is in the x-y
  ! plane ell(3) is set to complete the right-handed system.
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------------------------
  Subroutine fundcell_genlatticevects(fcell)
    Type(Fundamental_Cell), Intent(inout) :: fcell

    Real(kind=RDbl)                :: magvec3sq, angleab, anglebc, angleac
    Type(VecType)                  :: vec3
    Real(kind=RDbl), Dimension(3)  :: comp

    !** 1st-lattice vector
    comp(1) = fcell%ell(1)
    comp(2) = 0.0_RDbl
    comp(3) = 0.0_RDbl
    fcell%lvec(1) = comp

    !** 2nd-lattice vector
    angleab = fcell%angleab
    comp(1) = fcell%ell(2)*cos(angleab)
    comp(2) = fcell%ell(2)*sin(angleab)
    comp(3) = 0.0_RDbl
    fcell%lvec(2) = comp

    !** 3rd-lattice vector
    angleac = fcell%angleac
    anglebc = fcell%anglebc
    comp(1) = fcell%ell(3)*cos(angleac)
    comp(2) = fcell%ell(3)* & 
        (Cos(anglebc) - Cos(angleac)*Cos(angleab))/Sin(angleab)
    ! Figure out what component 3 should be so that the magnitude
    ! of the vector is ell(3)
    comp(3) = 0.0_RDbl
    vec3 = comp
    magvec3sq = vector_getnormsq(vec3)
    comp(3) = Sqrt(fcell%ell(3)**2 - magvec3sq)
    fcell%lvec(3) = comp

  End Subroutine fundcell_genlatticevects

  !----------------------------------------------------------------------------
  ! Returns a specific one of the three lattice vectors
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_latticevec(fcell,vecno)
    Type(Fundamental_Cell), Intent(InOut) :: fcell
    Integer, Intent(In)                   :: vecno

    fundcell_latticevec = fcell%lvec(vecno)

  End Function fundcell_latticevec

  !-----------------------------------------------------------
  ! Build the cell widths in each of the a,b,c directions
  ! Dot each lattice vector with the cross product of the
  ! other two.  This projects the lattice vector along the
  ! vector normal to the other two and gives the "width"
  ! Requires: fcell -- fundamental cell data structure
  !-----------------------------------------------------------
  Subroutine fundcell_genwidths(fcell)
    Type(Fundamental_Cell), Intent(InOut) :: fcell

    Integer         :: i,j,k
    Type(VecType)   :: crossprod

    Do i = 1,3
      j = Mod(i+1,3)
      k = Mod(i+2,3)
      If (j == 0) j = 3
      If (k == 0) k = 3
      crossprod = vector_crossprod(fcell%lvec(j),fcell%lvec(k))
      fcell%width(i) = (fcell%lvec(i)*crossprod)/vector_getnorm(crossprod)
      !      Write(*,*) i,j,k,fcell%width(i)
    End Do

    fcell%minwidth = minn(fcell%width)

  End Subroutine fundcell_genwidths

  !----------------------------------------------------------------------------
  ! Returns the cell width in a specified direction (a=1,b=2,c=3)
  ! Requires:  fcell -- fundamental cell data structure
  !            vecno -- width direction (cooresponds to lattice vector idx)
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function fundcell_width(fcell,vecno)
    Type(Fundamental_Cell), Intent(In) :: fcell
    Integer, Intent(In)                :: vecno

    fundcell_width = fcell%width(vecno)

  End Function fundcell_width

  !----------------------------------------------------------------------------
  ! Returns the minimum cell width
  ! Requires:  fcell -- fundamental cell data structure
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function fundcell_minwidth(fcell)
    Type(Fundamental_Cell), Intent(In) :: fcell

    fundcell_minwidth = fcell%minwidth

  End Function fundcell_minwidth

  !----------------------------------------------------------
  ! Gets the edge lengths of the fundamental cell
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------
  Function fundcell_getell(fcell)
    Type(Fundamental_Cell), Intent(in)      :: fcell
    Real(kind=RDbl), Dimension(3)           :: fundcell_getell

    fundcell_getell = fcell%ell

  End Function fundcell_getell

  !--------------------------------------------------------------------
  ! Gets the angles of the fundamental cell
  ! Requires: fcell -- fundamental cell data structure
  !           deg_flag -- flags if angles should be in degrees
  !--------------------------------------------------------------------
  Function fundcell_angles(fcell,deg_flag)
    Type(Fundamental_Cell), Intent(In)      :: fcell
    Logical, Intent(In), Optional           :: deg_flag
    Real(kind=RDbl), Dimension(3)           :: fundcell_angles

    fundcell_angles(1) = fcell%anglebc
    fundcell_angles(2) = fcell%angleac
    fundcell_angles(3) = fcell%angleab

    If (Present(deg_flag)) Then
      If (deg_flag) Then
        fundcell_angles = fundcell_angles/Pi*180.0_RDbl
      End If
    End If

  End Function fundcell_angles

  !----------------------------------------------------------
  ! Gets the effective box size of the fundamental cell
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------
  Function fundcell_geteff(fcell)
    Type(Fundamental_Cell), Intent(In) :: fcell
    Real(Kind=RDbl), Dimension(3)      :: fundcell_geteff

    fundcell_geteff = fcell%eff

  End Function fundcell_geteff

  !----------------------------------------------------------
  ! Get the volume of the fundamental unit cell
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------
  Real(kind=RDbl) Function fundcell_getvolume(fcell)
    Type(Fundamental_Cell), Intent(in)      :: fcell

    fundcell_getvolume = fcell%volume

  End Function fundcell_getvolume

  !----------------------------------------------------------
  ! Returns TRUE if the cell is orthorhombic
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------
  Logical Function fundcell_isortho(fcell)
    Type(Fundamental_Cell), Intent(in)      :: fcell
    fundcell_isortho = fcell%orthoflag
  End Function fundcell_isortho

  !----------------------------------------------------------------------
  ! Construct a larger cell from a more fundamental unit
  ! Requires: smallcell -- fundamental cell to build from
  !           size -- integer vector given multiples in each direction
  !           largecell -- larger fundamental cell to create
  !----------------------------------------------------------------------
  Subroutine fundcell_construct(smallcell,size,largecell)
    Type(Fundamental_Cell), Intent(In)   :: smallcell
    Type(IntVecType), Intent(In)         :: size
    Type(Fundamental_Cell), Intent(Out)  :: largecell

    Integer         :: i
    Real(Kind=RDbl) :: rsize  ! Real values of the int size

    !** Set the parameters for the new cell
    largecell = smallcell

    !** change edge lengths
    Do i = 1,3
      rsize = Real(size%comp(i),Kind=RDbl)
      largecell%ell(i) = smallcell%ell(i)*rsize
      largecell%half_ell(i) = smallcell%half_ell(i)*rsize
      largecell%ell_inv(i) = 1.0_RDbl/largecell%ell(i)
      largecell%eff(i) = smallcell%eff(i)*rsize
    End Do

    !** generate the lattice vectors, widths and volume
    Call fundcell_genlatticevects(largecell)
    Call fundcell_genwidths(largecell)
    Call fundcell_genvolume(largecell)

    !LC    Call fundcell_display(largecell,2,6)

  End Subroutine fundcell_construct

  !----------------------------------------------------------------------
  ! Initialize the transformation matrices
  ! Requires: fcell -- fundamental cell data structure
  !----------------------------------------------------------------------
  Subroutine fundcell_initxform(fcell)
    Type(Fundamental_Cell), Intent(InOut)   :: fcell    

    Real(kind=RSgl)                         :: alpha,beta,gamma
    Real(kind=RDbl)                         :: sign,sign2

    !** Calculate the slant to cartesian transformation matrix (unslant mtx)
    alpha = fcell%anglebc
    beta  = fcell%angleac
    gamma = fcell%angleab

    !** This calculation assumes that the (a) vector is co-linear with the
    !** cartesian x-axis and that the (b) vector is in the x-y plane
    sign = 1.0_RDbl         !+/- possible
    sign2 = 1.0_RDbl        !+/- possible

    fcell%unslantmatrix%comp(1,1) = 1.0_RDbl
    fcell%unslantmatrix%comp(1,2) = Cos(gamma)
    fcell%unslantmatrix%comp(1,3) = Cos(beta)

    fcell%unslantmatrix%comp(2,1) = 0.0_RDbl
    fcell%unslantmatrix%comp(2,2) = sign*Sqrt(1.0_RDbl - Cos(gamma)**2)
    fcell%unslantmatrix%comp(2,3) = (Cos(alpha) - Cos(gamma)*Cos(beta))/&
        fcell%unslantmatrix%comp(2,2)

    fcell%unslantmatrix%comp(3,1) = 0.0_RDbl
    fcell%unslantmatrix%comp(3,2) = 0.0_RDbl
    fcell%unslantmatrix%comp(3,3) = sign2*Sqrt(1.0_RDbl - &
        fcell%unslantmatrix%comp(1,3)**2 - &
        fcell%unslantmatrix%comp(2,3)**2)

    !** Invert to the get the slantmatrix
    fcell%slantmatrix = matrix_inverse(fcell%unslantmatrix)

  End Subroutine fundcell_initxform

  !----------------------------------------------------------------------------
  ! Transform a vector from the unslanted to slanted basis
  ! Requires: fcell -- fundamental cell data structure
  !           vec -- vector to transform
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_slant(fcell, vec)
    Type(Fundamental_Cell), Intent(in)   :: fcell
    Type(VecType), Intent(In)            :: vec

    If (fcell%orthoflag) Then
      fundcell_slant = vec
      Return
    End If

    fundcell_slant = fcell%slantmatrix * vec

  End Function fundcell_slant

  !----------------------------------------------------------------------------
  ! Transform a vector from the slanted to unslanted basis
  ! Requires: fcell -- fundamental cell data structure
  !           vec -- vector to transform
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_unslant(fcell, vec)
    Type(Fundamental_Cell), Intent(In)   :: fcell
    Type(VecType), Intent(In)            :: vec

    If (fcell%orthoflag) Then
      fundcell_unslant = vec
      Return
    End If

    fundcell_unslant = fcell%unslantmatrix * vec
  End Function fundcell_unslant

  !----------------------------------------------------------------------------
  ! Translate a vector to the primary fundamental cell
  ! Requires:  fcell -- the fundamental cell
  !            vec -- the vector to be transformed
  !            cvec -- optional (output only) continuation vector
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_maptocell(fcell,vec,cvec)
    Type(Fundamental_Cell), Intent(In)       :: fcell
    Type(VecType), Intent(In)                :: vec
    Type(IntVecType), Intent(Out), Optional  :: cvec

    Integer                            :: i
    Type(VecType)                      :: dummy

    dummy = vec - fcell%origin
    dummy = fundcell_slant(fcell,dummy)
    dummy%comp = dummy%comp * fcell%ell_inv

    If (Present(cvec)) Then
      !** also get the continuation vectors
      Do i = 1,3
        cvec%comp(i) = floor(dummy%comp(i))
        dummy%comp(i) = dummy%comp(i) - cvec%comp(i)
      End Do
    Else
      !** simply translate back to the primary unit cell
      Do i = 1,3
        If (dummy%comp(i) >= 0.0d0) Then
          dummy%comp(i) = mod(dummy%comp(i), 1.0d0)
        Else
          dummy%comp(i) = mod(dummy%comp(i), 1.0d0) + 1.0d0
        End If
      End Do
    End If

    dummy%comp = dummy%comp * fcell%ell
    dummy = fundcell_unslant(fcell,dummy)
    dummy = dummy + fcell%origin    

    fundcell_maptocell = dummy

  End Function fundcell_maptocell

  !----------------------------------------------------------------------------
  ! Get the coordinates of the center of the fundamental cell
  ! Requires:  fcell -- fundamental cell data structure
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_center(fcell)
    Type(Fundamental_Cell), Intent(In)   :: fcell

    fundcell_center = fcell%half_ell

    If (.Not. fcell%orthoflag) Then
      fundcell_center = fundcell_unslant(fcell,fundcell_center)
    End If

  End Function fundcell_center

  !----------------------------------------------------------------------------
  ! Returns the quadrant in which the point lies
  !
  !                  *-----------+-----------*
  !                 /           /           /|
  !                /     7     /     8     / |
  !               /           /           /  |
  !              +-----------+-----------+   |
  !             /           /           /|   |
  !            /     5     /     6     / |   +
  !           /           /           /  |  /|
  !          +-----------+-----------+   | / |
  !          |           |           |   |/  |
  !          |           |           |   +   |
  !          |           |           |  /| 4 |
  !     >    |           |           | / |   *
  !  3 /     |           |           |/  |  /
  !          +-----------+-----------+   | /
  !          |           |           |   |/
  !          |           |           |   +
  !          |     1     |     2     |  /
  !          |           |           | /
  !          |           |           |/
  !          *-----------+-----------*    pretty picture, eh?  Amit's an artist
  ! Requires: fcell -- fundamental cell data structure
  !           point -- coordinate vector to place in quadrant
  !----------------------------------------------------------------------------
  Integer Function fundcell_getquadrant(fcell,point)
    Type(Fundamental_Cell), Intent(In) :: fcell
    Type(VecType), Intent(In)          :: point

    Integer :: upperSlice, backSlice, rightSlice

    upperSlice = 0
    backSlice = 0
    rightSlice = 0

    If (fcell%orthoflag.EQV..False.) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' fundcell is not orthorhombic. Cannot getquadrant yet.'
      Stop              
    End If

    If (point%comp(1) >= fcell%ell(3)/2.0_RDbl) upperSlice = 4
    If (point%comp(2) >= fcell%ell(2)/2.0_RDbl) backSlice = 2
    If (point%comp(3) >= fcell%ell(1)/2.0_RDbl) rightSlice = 1

    fundcell_getquadrant = 1 + rightSlice + backSlice + upperSlice

  End Function fundcell_getquadrant

  !----------------------------------------------------------------------------
  ! Returns the passed position vector shifted into an image fundamental cell
  ! Requires:  fcell -- fundamental cell data structure
  !            vec -- input position vector 
  !            shift -- integer array specifying number of cells to shift away
  !----------------------------------------------------------------------------
  Type(VecType) Function fundcell_shift(fcell,vec,shift)
    Type(Fundamental_Cell), Intent(In)     :: fcell
    Type(VecType), Intent(In)              :: vec
    Integer, Dimension(3), Intent(In)      :: shift

    fundcell_shift = vec + shift(1)*fcell%lvec(1) + &
        fcell%lvec(2)*shift(2) + fcell%lvec(3)*shift(3)

  End Function fundcell_shift

  !------------------------------------------------------------------
  ! Get the minimum image distance between two position vectors
  ! Requires:  fcell -- fundamental cell data structure
  !            atvec1 -- 1st position vector
  !            atvec2 -- 2nd position vector
  !            sepvec -- minimum image separation vector
  !            dist2 -- square of minimum image distance 
  !------------------------------------------------------------------
  Subroutine fundcell_snglMinImage(fcell,atvec1,atvec2,sepvec,dist2)
    Type(Fundamental_Cell), Intent(In)     :: fcell    
    Type(VecType), Intent(In)              :: atvec1,atvec2
    Type(VecType), Intent(Out)             :: sepvec
    Real(kind=RDbl), Optional, Intent(Out) :: dist2

    Integer                :: i
    Logical                :: is_slanted
    Real(kind=RDbl)        :: d1, d2, d3, XL, YL, ZL
    Type(VecType)          :: vec1,vec2

    !** check whether fundamental cell is orthorhombic
    is_slanted = .True.
    If (fundcell_isortho(fcell)) is_slanted = .False.


    XL=fcell%ell(1)   ! full edge lengths
    YL=fcell%ell(2)
    ZL=fcell%ell(3)



    If (is_slanted) Then

      !vec1 = fundcell_slant(fcell,atvec1)
      !vec2 = fundcell_slant(fcell,atvec2)
      vec1 = atvec1
      vec2 = atvec2
        d1 = vec1%comp(1) - vec2%comp(1)
        d2 = vec1%comp(2) - vec2%comp(2)
        d3 = vec1%comp(3) - vec2%comp(3)
        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))

        sepvec%comp(1:3) = (/d1,d2,d3/)
        sepvec =  fundcell_unslant(fcell,sepvec)

        If (Present(dist2)) Then
          dist2 = vector_getnormsq(sepvec)
        End If

    else

        d1 = atvec1%comp(1) - atvec2%comp(1)
        d2 = atvec1%comp(2) - atvec2%comp(2)
        d3 = atvec1%comp(3) - atvec2%comp(3)

        ! I have tested "If loops" vs. the use of ANINT
        ! This seems slightly faster. See discusion in Allen and Tildesely
        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))
        sepvec%comp(1:3) = (/d1,d2,d3/)


        !** pass distance-squared if needed
        If (Present(dist2)) Then
          dist2 = d1*d1+d2*d2+d3*d3
        End If
!If (dbgflag) Write(*,*) d1, d2, d3, dist2

    Endif


  End Subroutine fundcell_snglMinImage

  !----------------------------------------------------------------------------
  ! Returns the minimum image distances between all coords in atvecs1 and 
  ! atvecs2.  Separation vectors are calculated as (atom1 - atom2).  
  ! sepvecs is indexed so (atom1,atom2) is the distance between point 
  ! atvecs1(atom1) and atvecs2(atom2).  Optionally returns the distance
  ! squared in dists2.  
  ! Requires:  fcell -- fundamental cell data structure
  !            atvecs1 -- 1st position vectors
  !            atvecs2 -- 2nd position vectors
  !            sepvecs -- minimum image separation vectors
  !            dists2 -- square of minimum image distances 
  !----------------------------------------------------------------------------
  Subroutine fundcell_multMinImage(fcell,atvecs1,atvecs2,sepvecs,dists2)
    Type(Fundamental_Cell), Intent(In)                     :: fcell    
    Type(VecType), Intent(In), Dimension(:)                :: atvecs1
    Type(VecType), Intent(In), Dimension(:)                :: atvecs2
    Type(VecType), Intent(Out), Dimension(:,:)             :: sepvecs
    Real(kind=RDbl), Dimension(:,:), Intent(Out), Optional :: dists2    

    Integer                  :: a1, a2, na1, na2
    Logical                  :: is_slanted
    Real(kind=RDbl)          :: d1,d2,d3,XL,YL,ZL
    Type(VecType)            :: vec1,vec2

    !** checks whether orthorhombic or not
    is_slanted = .True.
    If (fundcell_isortho(fcell)) is_slanted=.False.

    na1=Size(atvecs1,1)
    na2=Size(atvecs2,1)

    XL=fcell%ell(1)   ! full edge lengths
    YL=fcell%ell(2)
    ZL=fcell%ell(3)

    If (is_slanted) Then


      Do a1 = 1, na1
        vec1 = fundcell_slant(fcell,atvecs1(a1))
        Do a2 = 1, na2
          vec2 = fundcell_slant(fcell,atvecs2(a2))
          d1 = vec1%comp(1) - vec2%comp(1)
          d2 = vec1%comp(2) - vec2%comp(2)
          d3 = vec1%comp(3) - vec2%comp(3)
          d1=d1-XL*(Anint(d1/XL))
          d2=d2-YL*(Anint(d2/YL))
          d3=d3-ZL*(Anint(d3/ZL))

          sepvecs(a1,a2)%comp(1:3) = (/d1,d2,d3/)
          sepvecs(a1,a2) =  fundcell_unslant(fcell,sepvecs(a1,a2))

          If (Present(dists2)) Then
            dists2(a1,a2) = vector_getnormsq(sepvecs(a1,a2))
          End If

        End Do
      End Do

    else


      Do a1 = 1, na1
        Do a2 = 1, na2
          d1 = atvecs1(a1)%comp(1) - atvecs2(a2)%comp(1)
          d2 = atvecs1(a1)%comp(2) - atvecs2(a2)%comp(2)
          d3 = atvecs1(a1)%comp(3) - atvecs2(a2)%comp(3)

          ! I have tested "If loops" vs. the use of ANINT
          ! This seems slightly faster. See discusion in Allen and Tildesely
          d1=d1-XL*(Anint(d1/XL))
          d2=d2-YL*(Anint(d2/YL))
          d3=d3-ZL*(Anint(d3/ZL))

          sepvecs(a1,a2)%comp(1:3) = (/d1,d2,d3/)


          !** pass distance-squared if needed
          If (Present(dists2)) Then
            dists2(a1,a2) = d1*d1+d2*d2+d3*d3
          End If

        End Do
      End Do

    endif

  End Subroutine fundcell_multMinImage


  !----------------------------------------------------------------------------
  ! Returns the minimum image distances between one atom and a set of atoms
  ! the number of atoms is passed in nvecs
  !----------------------------------------------------------------------------
  Subroutine fundcell_multMinImage2(fcell,atvec1,atvecs2,sepvecs,nvecs,dists2)
    Type(Fundamental_Cell), Intent(In)                     :: fcell    
    Type(VecType), Intent(In)               :: atvec1
    Type(VecType), Intent(In), Dimension(:)                :: atvecs2
    Type(VecType), Intent(Out), Dimension(:)             :: sepvecs
    Integer, Intent(in) :: nvecs
    Real(kind=RDbl), Dimension(:), Intent(Out), Optional :: dists2    

    Integer                  :: a2
    Logical                  :: is_slanted
    Real(kind=RDbl)          :: d1,d2,d3,XL,YL,ZL
    Type(VecType)            :: vec1,vec2

    !** checks whether orthorhombic or not
    is_slanted = .True.
    If (fundcell_isortho(fcell)) is_slanted=.False.

    XL=fcell%ell(1)   ! full edge lengths
    YL=fcell%ell(2)
    ZL=fcell%ell(3)

    If (is_slanted) Then

      vec1 = fundcell_slant(fcell,atvec1)
      Do a2 = 1, nvecs

        vec2 = fundcell_slant(fcell,atvecs2(a2))
        d1 = vec1%comp(1) - vec2%comp(1)
        d2 = vec1%comp(2) - vec2%comp(2)
        d3 = vec1%comp(3) - vec2%comp(3)
        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))

        sepvecs(a2)%comp(1:3) = (/d1,d2,d3/)
        sepvecs(a2) =  fundcell_unslant(fcell,sepvecs(a2))

        If (Present(dists2)) Then
          dists2(a2) = vector_getnormsq(sepvecs(a2))
        End If

      End Do

    else


      Do a2 = 1, nvecs
        d1 = atvec1%comp(1) - atvecs2(a2)%comp(1)
        d2 = atvec1%comp(2) - atvecs2(a2)%comp(2)
        d3 = atvec1%comp(3) - atvecs2(a2)%comp(3)

        ! I have tested "If loops" vs. the use of ANINT
        ! This seems slightly faster. See discusion in Allen and Tildesely
        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))

        sepvecs(a2)%comp(1:3) = (/d1,d2,d3/)


        !** pass distance-squared if needed
        If (Present(dists2)) Then
          dists2(a2) = d1*d1+d2*d2+d3*d3
        End If

      End Do

    Endif

  End Subroutine fundcell_multMinImage2



  !----------------------------------------------------------------------------
  ! Returns the minimum image distances between one atom and a set of atoms
  ! the number of atoms is passed in nvecs
  ! list contains the index of atoms from atvecs2 for which calcualtions
  ! should be done
  !----------------------------------------------------------------------------
  Subroutine fundcell_listMinImage(fcell,atvec1,atvecs2,sepvecs,nvecs,&
      dists2,list)
    Type(Fundamental_Cell), Intent(In)                     :: fcell    
    Type(VecType), Intent(In)               :: atvec1
    Type(VecType), Intent(In), Dimension(:)                :: atvecs2
    Type(VecType), Intent(Out), Dimension(:)             :: sepvecs
    Integer, Intent(in) :: nvecs
    Real(kind=RDbl), Dimension(:), Intent(Out) :: dists2    
    Integer, Dimension(:), Intent(in) :: list

    Integer                  :: a2, l
    Logical                  :: is_slanted
    Real(kind=RDbl)          :: d1,d2,d3,XL,YL,ZL
    Type(VecType)            :: vec1,vec2

    !** checks whether orthorhombic or not
    is_slanted = .True.
    If (fundcell_isortho(fcell)) is_slanted=.False.

    XL=fcell%ell(1)   ! full edge lengths
    YL=fcell%ell(2)
    ZL=fcell%ell(3)

    If (is_slanted) Then
      vec1 = fundcell_slant(fcell,atvec1)
      Do l = 1, nvecs
        a2=list(l)
        vec2 = fundcell_slant(fcell,atvecs2(a2))
        d1 = vec1%comp(1) - vec2%comp(1)
        d2 = vec1%comp(2) - vec2%comp(2)
        d3 = vec1%comp(3) - vec2%comp(3)
        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))

        sepvecs(l)%comp(1:3) = (/d1,d2,d3/)
        sepvecs(l) =  fundcell_unslant(fcell,sepvecs(l))

        dists2(l) = vector_getnormsq(sepvecs(l))

      End Do

    else
      Do l = 1, nvecs
        a2=list(l)
        d1 = atvec1%comp(1) - atvecs2(a2)%comp(1)
        d2 = atvec1%comp(2) - atvecs2(a2)%comp(2)
        d3 = atvec1%comp(3) - atvecs2(a2)%comp(3)

        d1=d1-XL*(Anint(d1/XL))
        d2=d2-YL*(Anint(d2/YL))
        d3=d3-ZL*(Anint(d3/ZL))

        sepvecs(l)%comp(1:3) = (/d1,d2,d3/)
        dists2(l) = d1*d1+d2*d2+d3*d3
      End Do

    Endif

  End Subroutine fundcell_listMinImage





  !----------------------------------------------------------------------------
  ! Returns the minimum image distances between all coords in atvecs1
  ! only uper triangel of utput is filled
  !----------------------------------------------------------------------------
  Subroutine fundcell_singleArrayMinImage(fcell,atvecs,sepvecs,dists2,nvecs)
    Type(Fundamental_Cell), Intent(In)                     :: fcell    
    Type(VecType), Intent(In), Dimension(:)                :: atvecs
    Type(VecType), Intent(Out), Dimension(:,:)             :: sepvecs
    Real(kind=RDbl), Dimension(:,:), Intent(Out) :: dists2    
    Integer,Intent(in) :: nvecs

    Integer                  :: a1, a2
    Logical                  :: is_slanted
    Real(kind=RDbl)          :: d1,d2,d3
    Type(VecType)            :: vec1,vec2

    !** checks whether orthorhombic or not
    is_slanted = .True.
    If (fundcell_isortho(fcell)) is_slanted=.False.


    Do a1 = 1, nvecs
      Do a2 = a1+1, nvecs


        !** Apply the minimum image convention
        If (is_slanted) Then
          !** Transform to the slanted basis
          vec1 = fundcell_slant(fcell,atvecs(a1))
          vec2 = fundcell_slant(fcell,atvecs(a2))
          d1 = vec1%comp(1) - vec2%comp(1)
          d2 = vec1%comp(2) - vec2%comp(2)
          d3 = vec1%comp(3) - vec2%comp(3)
Write(*,*) "is dists2 correct in this routine?"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop

        Else
          d1 = atvecs(a1)%comp(1) - atvecs(a2)%comp(1)
          d2 = atvecs(a1)%comp(2) - atvecs(a2)%comp(2)
          d3 = atvecs(a1)%comp(3) - atvecs(a2)%comp(3)
        End If

        If (d1 > fcell%half_ell(1)) Then
          d1 = d1 - fcell%ell(1)
        Else if (d1 < -fcell%half_ell(1)) Then
          d1 = d1 + fcell%ell(1)
        End If
        If (d2 > fcell%half_ell(2)) Then
          d2 = d2 - fcell%ell(2)
        Else If (d2 < -fcell%half_ell(2)) Then
          d2 = d2 + fcell%ell(2)
        End If
        If (d3 > fcell%half_ell(3)) Then
          d3 = d3 - fcell%ell(3)
        Else If (d3 < -fcell%half_ell(3)) Then
          d3 = d3 + fcell%ell(3)
        End If

        sepvecs(a1,a2)%comp(1) = d1
        sepvecs(a1,a2)%comp(2) = d2
        sepvecs(a1,a2)%comp(3) = d3

        !** Transform back to the cartesian system to get the distance
        !** not needed if fundamental cell is orthogonal
        If (is_slanted) sepvecs(a1,a2) = &
            fundcell_unslant(fcell,sepvecs(a1,a2))

        !** pass distance-squared if needed
        dists2(a1,a2) = d1*d1+d2*d2+d3*d3

      End Do
    End Do

  End Subroutine fundcell_singleArrayMinImage


  !----------------------------------------------------------------------------
  ! Dump the edge lengths and angles into a string
  ! Requires:  fcell -- fundamental cell data structure
  !----------------------------------------------------------------------------
  Function fundcell_cellstring(fcell)
    Character(len=lstrLen)                  :: fundcell_cellstring
    Type(Fundamental_Cell), Intent(In)      :: fcell    

    Write(fundcell_cellstring,'(6f10.5)') fcell%ell,fcell%anglebc/degTorad, &
        fcell%angleac/degTorad,fcell%angleab/degTorad

  End Function fundcell_cellstring

  !----------------------------------------------------------------------------
  ! Display the fundamental cell information
  ! Requires: fcell -- fundamental cell data structure
  !           indent -- no. of spaces from the left margin
  !           unit -- unit to dump into
  !----------------------------------------------------------------------------
  Subroutine fundcell_display(fcell,indent,unit)
    Type(Fundamental_Cell), Intent(In)      :: fcell    
    Integer, Intent(In)                     :: indent
    Integer, Intent(In)                     :: unit

    Integer                :: i
    Character(len=indent)  :: blank
    Character(len=lstrLen) :: string

    blank = Repeat(' ',indent)

    If (fcell%orthoflag) Then
      Write(unit,'(2a)') blank,'ORTHORHOMBIC'
    Else
      Write(unit,'(2a)') blank,'NON-ORTHORHOMBIC'
    End If

    Write(unit,'(2a,3f10.4)') blank,'edge lengths (Ang)    : ', &
        (fcell%ell(i),i=1,3)
    Write(unit,'(2a,3f10.4)') blank,'cell angles (degrees) : ', &
        fcell%anglebc/Pi*180, fcell%angleac/Pi*180, &
        fcell%angleab/Pi*180
    !    Write(unit,'(2a,3f10.4)') blank,'bound box size        : ', &
    !        (fcell%eff(i), i=1,3)
    string = vector_display(fcell%origin,'f9.4')
    Write(unit,'(3a)')        blank,'origin shift (Ang)    : ', &
        Trim(string)
    string = vector_display(fcell%lvec(1),'f9.4')
    Write(unit,'(3a)')        blank,'cell vector a         : ', &
        Trim(string)
    string = vector_display(fcell%lvec(2),'f9.4')
    Write(unit,'(3a)')        blank,'cell vector b         : ', &
        Trim(string)
    string = vector_display(fcell%lvec(3),'f9.4')
    Write(unit,'(3a)')        blank,'cell vector c         : ', &
        Trim(string)
    Write(unit,'(2a,3f10.4)') blank,'box width in (a,b,c)  : ', &
        (fcell%width(i), i=1,3)
    Write(unit,'(2a,f10.4)')  blank,'minimum width (Ang)   : ', &
        fcell%minwidth
    string = real2str(fcell%volume,8)
    Write(unit,'(3a)')  blank,'volume (Angstroms^3)  : ',Trim(string)

  End Subroutine fundcell_display

End Module fundcell



