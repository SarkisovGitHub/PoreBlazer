!----------------------------------------------------------
! This module defines various matrix operations.  Although
! I have tried to keep the matrix dimensions general a lot
! of routines work on only 3x3 matrices.
!----------------------------------------------------------
Module matrix

  Use defaults, Only: RDbl, strLen, one, zero, pi, twopi
  Use utils, Only: int2str,real2str
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb
  Use random, Only:   rranf


  Implicit None
  Save

  Private
  Public :: MatrixType, Assignment(=), Operator(*), matrix_genrotcx, &
      matrix_genrotcy, matrix_genrotcz, matrix_display, matrix_inverse, &
      matrix_getinv, matrix_Db2fDPsi, matrix_b2feuler, matrix_f2beuler, &
      matrix_Db2fDTheta, &
      matrix_Db2fDPhi, matrix_transpose, matrix_identity, invertarray, &
      matrix_rotAroundVec, matrix_rotate_random

  Integer, Parameter   :: nrows = 3
  Integer, Parameter   :: ncols = 3
  
  Type MatrixType
    Real(kind=RDbl), Dimension(nrows, ncols)  :: comp
  End Type MatrixType

  Interface Assignment(=)
    Module Procedure matrix_m_eq_real
    Module Procedure matrix_m_eq_2darray
    Module Procedure matrix_2darray_eq_m
  End Interface

  Interface Operator(*)
    Module Procedure matrix_m_mult_v
    Module Procedure matrix_v_mult_m
    Module Procedure matrix_m_mult_m
  End Interface

  Interface matrix_genrotcx
    Module Procedure matrix_genrotcx_trig
    Module Procedure matrix_genrotcx_angle
  End Interface

  Interface matrix_genrotcy
    Module Procedure matrix_genrotcy_trig
    Module Procedure matrix_genrotcy_angle
  End Interface

  Interface matrix_genrotcz
    Module Procedure matrix_genrotcz_trig
    Module Procedure matrix_genrotcz_angle
  End Interface

  Interface display
    Module Procedure matrix_display
  End Interface

Contains

  !-------------------------------------------------------------------------
  ! Generates an identity matrix
  !-------------------------------------------------------------------------
  Type(MatrixType) Function matrix_identity()
    
    matrix_identity = 0.0_RDbl
    matrix_identity%comp(1,1) = 1.0_RDbl
    matrix_identity%comp(2,2) = 1.0_RDbl
    matrix_identity%comp(3,3) = 1.0_RDbl

  End Function matrix_identity

  !-------------------------------------------------------
  ! This function initializes the matrix from a 2-D array
  !-------------------------------------------------------
  Subroutine matrix_m_eq_2darray(m1, array)
    Type(MatrixType), Intent(out) :: m1
    Real(kind=RDbl), Dimension(:,:), Intent(in)   :: array
    Integer  :: r, c, dim1, dim2

    dim1 = Size(array,1)
    dim2 = Size(array,2)
    Do r=1, dim1
      Do c=1, dim2
        m1%comp(r, c) = array(r, c)
      End Do
    End Do
    Return
  End Subroutine matrix_m_eq_2darray

  !-------------------------------------------------------
  ! This function sets a 2-D array from a matrix
  !-------------------------------------------------------
  Subroutine matrix_2darray_eq_m(array, m1)
    Real(kind=RDbl), Dimension(:,:), Intent(Out)  :: array
    Type(MatrixType), Intent(In)                  :: m1

    Integer  :: r, c, dim1, dim2

    dim1 = Size(m1%comp,1)
    dim2 = Size(m1%comp,2)
    Do r = 1,dim1
      Do c = 1,dim2
        array(r,c) = m1%comp(r,c)
      End Do
    End Do

  End Subroutine matrix_2darray_eq_m

  !-------------------------------------------------------
  ! This function assigns the scalar "r1" to each element 
  ! of the matrix "m1" 
  !-------------------------------------------------------
  Subroutine matrix_m_eq_real(m1, r1)
    Type(MatrixType), Intent(out) :: m1
    Real(kind=RDbl), Intent(in):: r1
    Integer     :: r, c
    
    Do c=1, ncols
      Do r=1, nrows
        m1%comp(r, c) = r1
      End Do
    End Do
  End Subroutine matrix_m_eq_real

  !--------------------------------------------------
  ! This function multiplies the matrix "m1" by the
  ! vector "vec1".  Could have written in terms of
  ! dot products but vectors and matrices are so
  ! intimate that they have access to each others
  ! private members.
  !--------------------------------------------------
  Type(VecType) Function matrix_m_mult_v(m1, vec1)
    Type(MatrixType), Intent(in) :: m1
    Type(VecType), Intent(in)    :: vec1
    Integer            :: r, c
    
    matrix_m_mult_v = 0.0_RDbl
    Do r=1, nrows
      Do c=1, ncols
        matrix_m_mult_v%comp(r) = matrix_m_mult_v%comp(r) + &
            m1%comp(r,c)*vec1%comp(c)
      End Do
    End Do
    Return
  End Function matrix_m_mult_v
  
  !---------------------------------------------------
  ! This function multiplies the vector "vec1" by the
  ! matrix "m1"
  !---------------------------------------------------
  Type(VecType) Function matrix_v_mult_m(vec1, m1)
    Type(MatrixType), Intent(in) :: m1
    Type(VecType), Intent(in)    :: vec1
    Type(VecType)      :: temp
    Integer            :: r, c
    
    temp = 0.0_RDbl
    Do c=1, ncols
      Do r=1, nrows
        temp%comp(c) = temp%comp(c) + vec1%comp(r)*m1%comp(r,c)
      End Do
    End Do
    matrix_v_mult_m = temp
    Return
  End Function matrix_v_mult_m

  !-------------------------------------------------
  ! Multiplies matrices "m1" and "m2"
  !-------------------------------------------------
  Type(MatrixType) Function matrix_m_mult_m(m1, m2)
    Type(MatrixType), Intent(in) :: m1, m2
    Type(MatrixType)   :: temp
    Integer            :: r, c, i
    
    temp = 0.0_RDbl
    Do r=1, nrows
      Do c=1, ncols
        Do i=1, ncols
          temp%comp(r,c) = temp%comp(r,c) + m1%comp(r,i)*m2%comp(i,c)
        End Do
      End Do
    End Do
    matrix_m_mult_m = temp
    Return
  End Function matrix_m_mult_m

  !----------------------------------------------------------------------------
  ! Inverts a 3x3 matrix numerically
  !----------------------------------------------------------------------------
  Type(MatrixType) Function matrix_inverse(m1)
    Implicit None
    Type(MatrixType), Intent(in) :: m1

    Integer                      :: i,j
    !variables for the NumRec subroutines
    Integer                      :: np,n
    Integer, Dimension(3)        :: indx
    Real(kind=RDbl)              :: a(3,3),y(3,3),d

    !** Invert the matrix using techniques from Numerical Recipes
    np = 3
    n = 3
    d = 1.0d0
    
    !setup an identity matrix
    Do i = 1,n
       Do j = 1,n
          y(i,j) = 0
       EndDo
       y(i,i) = 1
    EndDo
    
    !** Perform LU decomposition 
    a = m1%comp
    Call matrixops_ludcmp(a,n,np,indx,d)
    
    !** find the inverse by solving an equation for each column
    Do j = 1,n
       Call matrixops_lubksb(a,n,np,indx,y(1:n,j))
    EndDo
    
    matrix_inverse%comp = y

  End Function matrix_inverse

  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM by an angle "theta" about the X-axis.  This matrix
  ! when multiplied by a vector gives the components of the
  ! vector in the rotated coordinate system.
  ! See Notes on coordinate rotation at the end of this Module
  ! L_theta
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcx_angle(theta)
    Real(kind=RDbl), Intent(in) :: theta
    Real(kind=RDbl)  :: sintheta, costheta
    
    sintheta = Sin(theta)
    costheta = Cos(theta)
    matrix_genrotcx_angle = matrix_genrotcx_trig(sintheta, costheta)
    Return
  End Function matrix_genrotcx_angle
  
  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM Anti-CW by an angle "theta" about the X-axis given
  ! sintheta and costheta.  This matrix when multiplied by a 
  ! vector gives the components of the vector in the rotated 
  ! coordinate system.  That is, the matrix expresses coordinates
  ! from the fixed to the body coordinate system
  ! See note at end of this module- Shaji.
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcx_trig(sintheta, costheta)
    Real(kind=RDbl), Intent(in)  :: sintheta, costheta
    
    matrix_genrotcx_trig%comp(1, 1:3) =  (/1.0_RDbl, 0.0_RDbl, 0.0_RDbl/)
    matrix_genrotcx_trig%comp(2, 1:3) =  (/0.0_RDbl, costheta, sintheta/)
    matrix_genrotcx_trig%comp(3, 1:3) =  (/0.0_RDbl, -sintheta, costheta/)
    Return
  End Function matrix_genrotcx_trig

  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM by an angle "theta" about the Y-axis.  This matrix
  ! when multiplied by a vector gives the components of the
  ! vector in the rotated coordinate system.
  ! See note at end of this module- Shaji.
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcy_angle(theta)
    Real(kind=RDbl), Intent(in) :: theta
    Real(kind=RDbl)  :: sintheta, costheta
    
    sintheta = Sin(theta)
    costheta = Cos(theta)
    matrix_genrotcy_angle = matrix_genrotcy_trig(sintheta, costheta)
    Return
  End Function matrix_genrotcy_angle
  
  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM Anti-CW by an angle "theta" about the Y-axis given
  ! sintheta and costheta.  This matrix when multiplied by a 
  ! vector gives the components of the vector in the rotated 
  ! coordinate system.
  ! See note at end of this module- Shaji.
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcy_trig(sintheta, costheta)
    Real(kind=RDbl), Intent(in)  :: sintheta, costheta

    matrix_genrotcy_trig%comp(1, 1:3) =  (/costheta,  0.0_RDbl,-sintheta/)
    matrix_genrotcy_trig%comp(2, 1:3) =  (/0.0_RDbl,1.0_RDbl,0.0_RDbl/)
    matrix_genrotcy_trig%comp(3, 1:3) =  (/sintheta,  0.0_RDbl, costheta/)
    Return
  End Function matrix_genrotcy_trig

  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM by an angle "theta" about the Z-axis.  This matrix
  ! when multiplied by a vector gives the components of the
  ! vector in the rotated coordinate system.
  ! See note at end of this module- Shaji.
  ! L_phi, L_psi
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcz_angle(theta)
    Real(kind=RDbl), Intent(in) :: theta
    Real(kind=RDbl)  :: sintheta, costheta
    
    sintheta = Sin(theta)
    costheta = Cos(theta)
    matrix_genrotcz_angle = matrix_genrotcz_trig(sintheta, costheta)
    Return
  End Function matrix_genrotcz_angle
  
  !-----------------------------------------------------------
  ! This function generates a matrix to rotate a COORDINATE
  ! SYSTEM Anti-CW by an angle "theta" about the Z-axis given
  ! sintheta and costheta.  This matrix when multiplied by a 
  ! vector gives the components of the vector in the rotated 
  ! coordinate system.
  ! See note at end of this module- Shaji.
  !-----------------------------------------------------------
  Type(MatrixType) Function matrix_genrotcz_trig(sintheta, costheta)
    Real(kind=RDbl), Intent(in)  :: sintheta, costheta

    matrix_genrotcz_trig%comp(1, 1:3) =  (/costheta,  sintheta, 0.0_RDbl /)
    matrix_genrotcz_trig%comp(2, 1:3) =  (/-sintheta, costheta, 0.0_RDbl /)
    matrix_genrotcz_trig%comp(3, 1:3) =  (/0.0_RDbl,0.0_RDbl,1.0_RDbl/)
    Return
  End Function matrix_genrotcz_trig

  !--------------------------------------------------------------------
  ! Generates the matrix that takes coordinates from the FIXED
  ! COORDINATE system to the BODY COORDINATE system given the
  ! 3 eulerian angles (phi, theta, psi) describing the orientation
  ! of the body coordinate sytem.  The convention is that of
  ! "Classical Dynamics of Particles and Systems by Marion & Thornton"
  ! 4th ed., pg. 432. Rotate by phi about z, theta about x, psi about z 
  ! See note at end of this module- Shaji.
  !  Result : L_psi * L_theta * L_phi
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_f2beuler(phi, theta, psi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, psi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = sin(psi)
    cospsi = cos(psi)
    matrix_f2beuler%comp(1, 1) = cospsi*cosphi - costheta*sinphi*sinpsi
    matrix_f2beuler%comp(1, 2) = cospsi*sinphi + costheta*cosphi*sinpsi
    matrix_f2beuler%comp(1, 3) = sinpsi*sintheta
    matrix_f2beuler%comp(2, 1) = -sinpsi*cosphi - costheta*sinphi*cospsi
    matrix_f2beuler%comp(2, 2) = -sinpsi*sinphi + costheta*cosphi*cospsi
    matrix_f2beuler%comp(2, 3) = cospsi*sintheta
    matrix_f2beuler%comp(3, 1) = sintheta*sinphi
    matrix_f2beuler%comp(3, 2) = -sintheta*cosphi
    matrix_f2beuler%comp(3, 3) = costheta

  End Function matrix_f2beuler




  !--------------------------------------------------------------------
  ! Generates the matrix that takes coordinates from the BODY
  ! COORDINATE system to the FIXED COORDINATE system given the
  ! 3 eulerian angles (phi, theta, psi) describing the orientation
  ! of the BODY coordinate sytem.  The convention is that of
  ! "Classical Dynamics of Particles and Systems by Marion & Thornton"
  ! 4th ed., pg. 432. Rotate by phi about z, theta about x, psi about z. 
  ! This matrix is essentially the transpose of the "b2feuler" matrix.
  !
  ! Multiply a vector by this matrix to rotate it by phi, theta, psi. 
  ! See note at end of this module- Shaji.
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_b2feuler(phi, theta, psi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, psi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = sin(psi)
    cospsi = cos(psi)
    matrix_b2feuler%comp(1, 1) = cospsi*cosphi  - costheta*sinphi*sinpsi
    matrix_b2feuler%comp(1, 2) = -sinpsi*cosphi - costheta*sinphi*cospsi
    matrix_b2feuler%comp(1, 3) = sintheta*sinphi
    matrix_b2feuler%comp(2, 1) = cospsi*sinphi  + costheta*cosphi*sinpsi
    matrix_b2feuler%comp(2, 2) = -sinpsi*sinphi + costheta*cosphi*cospsi
    matrix_b2feuler%comp(2, 3) = -sintheta*cosphi
    matrix_b2feuler%comp(3, 1) = sinpsi*sintheta
    matrix_b2feuler%comp(3, 2) = cospsi*sintheta
    matrix_b2feuler%comp(3, 3) = costheta
    Return
  End Function matrix_b2feuler

  !--------------------------------------------------------------------
  ! Generates the matrix that rotates any vector by "dpsi" around unit 
  ! vector defined by phi, theta ( all angles are measure anti clockwise 
  ! See note at end of this module- Shaji.
  !  Result : [L_(-phi) L_(-theta)] * [ L_(-dpsi) * L_theta * L_phi ]
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_rotAroundVec(phi, theta, dpsi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, dpsi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi
    Type(MatrixType) :: M1, M2

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = Sin(-dpsi)  !** NOTE -sign
    cospsi = Cos(-dpsi)  !** NOTE -sign
    ! same as f2beuler
    M1%comp(1, 1) = cospsi*cosphi - costheta*sinphi*sinpsi
    M1%comp(1, 2) = cospsi*sinphi + costheta*cosphi*sinpsi
    M1%comp(1, 3) = sinpsi*sintheta
    M1%comp(2, 1) = -sinpsi*cosphi - costheta*sinphi*cospsi
    M1%comp(2, 2) = -sinpsi*sinphi + costheta*cosphi*cospsi
    M1%comp(2, 3) = cospsi*sintheta
    M1%comp(3, 1) = sintheta*sinphi
    M1%comp(3, 2) = -sintheta*cosphi
    M1%comp(3, 3) = costheta

    ! L(-phi) * L(-theta)
    M2%comp(1, 1) =  cosphi
    M2%comp(1, 2) = -sinphi * costheta
    M2%comp(1, 3) =  sinphi * sintheta
    M2%comp(2, 1) =  sinphi
    M2%comp(2, 2) =  cosphi * costheta
    M2%comp(2, 3) = -cosphi * sintheta
    M2%comp(3, 1) =  zero 
    M2%comp(3, 2) =  sintheta
    M2%comp(3, 3) =  costheta

    matrix_rotAroundVec = M2 * M1 

  End Function matrix_rotAroundVec



  !--------------------------------------------------------------------
  ! Generates the matrix that gives the derivative of the b2f matrix 
  ! wrt to theta.  The b2f matrix takes the coordinates in the BODY
  ! COORDINATE system to the FIXED COORDINATE system given the
  ! 3 eulerian angles (phi, theta, psi) describing the orientation
  ! of the BODY coordinate sytem.  The convention is that of
  ! "Classical Dynamics of Particles and Systems by Marion & Thornton"
  ! 4th ed., pg. 432. Rotate by phi about z, theta about x, psi about z. 
  ! This matrix is essentially the transpose of the "b2feuler" matrix.
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_Db2fDTheta(phi, theta, psi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, psi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = sin(psi)
    cospsi = cos(psi)
    matrix_Db2fDTheta%comp(1, 1) = sintheta*sinphi*sinpsi
    matrix_Db2fDTheta%comp(1, 2) = sintheta*sinphi*cospsi
    matrix_Db2fDTheta%comp(1, 3) = costheta*sinphi
    matrix_Db2fDTheta%comp(2, 1) = -sintheta*cosphi*sinpsi
    matrix_Db2fDTheta%comp(2, 2) = -sintheta*cosphi*cospsi
    matrix_Db2fDTheta%comp(2, 3) = -costheta*cosphi
    matrix_Db2fDTheta%comp(3, 1) = sinpsi*costheta
    matrix_Db2fDTheta%comp(3, 2) = cospsi*costheta
    matrix_Db2fDTheta%comp(3, 3) = -sintheta
    Return
  End Function matrix_Db2fDTheta

  !--------------------------------------------------------------------
  ! Generates the matrix that gives the derivative of the b2f matrix 
  ! wrt to phi.  The b2f matrix takes the coordinates in the BODY
  ! COORDINATE system to the FIXED COORDINATE system given the
  ! 3 eulerian angles (phi, theta, psi) describing the orientation
  ! of the BODY coordinate sytem.  The convention is that of
  ! "Classical Dynamics of Particles and Systems by Marion & Thornton"
  ! 4th ed., pg. 432. Rotate by phi about z, theta about x, psi about z. 
  ! This matrix is essentially the transpose of the "b2feuler" matrix.
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_Db2fDPhi(phi, theta, psi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, psi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = sin(psi)
    cospsi = cos(psi)
    matrix_Db2fDPhi%comp(1, 1) = -cospsi*sinphi - costheta*cosphi*sinpsi
    matrix_Db2fDPhi%comp(1, 2) = sinpsi*sinphi  - costheta*cosphi*cospsi
    matrix_Db2fDPhi%comp(1, 3) = sintheta*cosphi
    matrix_Db2fDPhi%comp(2, 1) = cospsi*cosphi  - costheta*sinphi*sinpsi
    matrix_Db2fDPhi%comp(2, 2) = -sinpsi*cosphi  - costheta*sinphi*cospsi
    matrix_Db2fDPhi%comp(2, 3) = sintheta*sinphi
    matrix_Db2fDPhi%comp(3, 1) = 0.0_RDbl
    matrix_Db2fDPhi%comp(3, 2) = 0.0_RDbl
    matrix_Db2fDPhi%comp(3, 3) = 0.0_RDbl
    Return
  End Function matrix_Db2fDPhi

  !--------------------------------------------------------------------
  ! Generates the matrix that gives the derivative of the b2f matrix 
  ! wrt to psi.  The b2f matrix takes the coordinates in the BODY
  ! COORDINATE system to the FIXED COORDINATE system given the
  ! 3 eulerian angles (phi, theta, psi) describing the orientation
  ! of the BODY coordinate sytem.  The convention is that of
  ! "Classical Dynamics of Particles and Systems by Marion & Thornton"
  ! 4th ed., pg. 432. Rotate by phi about z, theta about x, psi about z. 
  ! This matrix is essentially the transpose of the "b2feuler" matrix.
  !--------------------------------------------------------------------
  Type(MatrixType) Function matrix_Db2fDPsi(phi, theta, psi)
    Real(kind=RDbl), Intent(in)  :: phi, theta, psi
    Real(kind=RDbl)    :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi

    sinphi = sin(phi)
    cosphi = cos(phi)
    sintheta = sin(theta)
    costheta = cos(theta)
    sinpsi = sin(psi)
    cospsi = cos(psi)
    matrix_Db2fDPsi%comp(1, 1) = -sinpsi*cosphi  - costheta*sinphi*cospsi
    matrix_Db2fDPsi%comp(1, 2) = -cospsi*cosphi + costheta*sinphi*sinpsi
    matrix_Db2fDPsi%comp(1, 3) = 0.0_RDbl
    matrix_Db2fDPsi%comp(2, 1) = -sinpsi*sinphi + costheta*cosphi*cospsi
    matrix_Db2fDPsi%comp(2, 2) = -cospsi*sinphi - costheta*cosphi*sinpsi
    matrix_Db2fDPsi%comp(2, 3) = 0.0_RDbl
    matrix_Db2fDPsi%comp(3, 1) = cospsi*sintheta
    matrix_Db2fDPsi%comp(3, 2) = -sinpsi*sintheta
    matrix_Db2fDPsi%comp(3, 3) = 0.0_RDbl
    Return
  End Function matrix_Db2fDPsi
  
  !----------------------------------------------------------
  ! Gets the transpose of a matrix "m1"
  !----------------------------------------------------------
  Type(MatrixType) Function matrix_transpose(m1)
    Type(MatrixType), Intent(in)  :: m1
    Integer     :: nr, nc, r, c
    
    nr = Size(m1%comp, 1)
    nc = Size(m1%comp, 2)
    Do r=1, nr
      Do c=1, nc
        matrix_transpose%comp(c,r) = m1%comp(r,c)
      Enddo
    End Do
    Return
  End Function matrix_transpose

  !----------------------------------------------------------
  ! Gets the determinant of a matrix of dimension "size X size"
  ! The default size is 3.
  !----------------------------------------------------------
  Recursive Function matrix_getdet(m1, size) Result(det)
    Type(MatrixType), Intent(in)  :: m1
    Integer, Optional, Intent(in) :: size
    Real(kind=RDbl)      :: det   
    Type(MatrixType)     :: minor
    Integer              :: dim, c, sign
    
    If (Present(size)) Then
      dim = size
    Else
      dim = 3
    End If

    If (dim == 1) Then
      det = m1%comp(1,1)
      Return
    Endif

    det = 0.0_RDbl
    Do c=1, dim
      minor = matrix_getminor(m1, 1, c, dim)
      sign = (-1)**(1+c)
      det = det + m1%comp(1,c)*sign*matrix_getdet(minor, dim-1)
    End Do
    Return
  End Function matrix_getdet

  !-----------------------------------------------------
  ! Takes a matrix "m1" of dimension "dim X dim" and returns 
  ! the minor of element (r,c).
  ! The minor is matrix resulting from removing row "r" and
  ! and column "c"
  !-----------------------------------------------------
  Type(MatrixType) Function matrix_getminor(m1, row, col, dim)
    Type(MatrixType), Intent(in)  :: m1
    Integer, Intent(in)           :: row, col
    Integer, Intent(in)           :: dim
    Integer    :: r, c, rowno, colno

    rowno = 0
    matrix_getminor = 0.0_RDbl
    Do r=1, dim
      If (r == row) Cycle
      rowno = rowno + 1
      colno = 0
      Do c=1, dim
        If (c == col) Cycle
        colno = colno + 1
        matrix_getminor%comp(rowno, colno) = m1%comp(r, c)
      End Do
    End Do
    Return
  End Function matrix_getminor

  !----------------------------------------------------------
  ! Gets the inverse of a matrix "m1"
  !----------------------------------------------------------
  Type(MatrixType) Function matrix_getinv(m1)
    Type(MatrixType), Intent(in)  :: m1
    Type(MatrixType) :: minor
    Real(kind=RDbl)  :: det, minordet
    Integer     :: r, c, dim, sign

    matrix_getinv = 0.0_RDbl
    dim = Size(m1%comp, 1)
    det = matrix_getdet(m1, dim)
    Do r=1, dim
      Do c=1, dim
        minor = matrix_getminor(m1, r, c, dim)
        minordet = matrix_getdet(minor, dim-1)
        sign = (-1)**(r+c)
        matrix_getinv%comp(c,r) = sign*minordet/det
      End Do
    End Do
    Return
  End Function matrix_getinv

  !----------------------------------------------------------------------------
  ! Dumps the contents of a matrix to an optional unit using a specified format
  ! Requires: m1 -- matrix to display
  !           fmt -- format of entries, if 'a', will pick the shortest 
  !           unitno -- unit to dump into
  !----------------------------------------------------------------------------
  Subroutine matrix_display(m1,fmt,unitno)
    Type(MatrixType), Intent(In)        :: m1
    Character(*), Intent(In)            :: fmt
    Integer, Optional, Intent(In)       :: unitno

    Integer                   :: r, c, unit
    Character(len=strLen)     :: strformat,string

    If (Present(unitno)) Then
      unit = unitno
    Else
      unit = 6
    End If
    string = int2str(ncols)
    Write(strformat, '(4a)') '(1x,',Trim(string),Trim(fmt),')'

    Do r = 1, nrows
      If (Trim(fmt) == 'a') Then
        string = real2str(m1%comp(r,c),1)
        Write(unit, strformat) (Trim(string)//' ', c=1, ncols)
      Else
        Write(unit, strformat) (m1%comp(r,c), c=1, ncols)
      End If
    End Do

  End Subroutine matrix_display

  !----------------------------------------------------------------
  ! This routine takes a square 2-dimensional array and inverts it
  ! The inversion is done by calling the lubksb and ludcmp routines
  ! of Press et al.
  !----------------------------------------------------------------
  Subroutine invertarray(arr, invarr)
    Real(kind=RDbl), Dimension(:,:), Intent(in) :: arr
    Real(kind=RDbl), Dimension(:,:), Intent(out):: invarr

    Integer     :: i, dim1, dim2
    Real(kind=RDbl), Dimension(Size(arr,1), Size(arr,2)) :: temparr
    Real(kind=RDbl), Dimension(Size(arr,1), Size(arr,2)) :: unitarr
    Real(kind=RDbl), Dimension(Size(arr,1)) :: rowperms
    Real(kind=RDbl)      :: nperms
    
    !** Make sure the array is square
    dim1 = Size(arr, 1)
    dim2 = Size(arr, 2)
    If (dim1 /= dim2) Then
      Write(0,'(1x,2a,i4, a)') __FILE__," : ",__LINE__, &
          " The array is not square"
      Stop
    End If

    !** Set up a unitary array
    unitarr = 0.0_RDbl
    Do i=1, dim1
      unitarr(i,i) = 1.0_RDbl
    End Do

    !** Do the decomposition to the LU form
    temparr = arr
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    !problem with NAG compiler
    Stop
!    Call ludcmp(temparr, dim1, dim1, rowperms, nperms)

    !** Do the back substitution for each column of the unitary array
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    !problem with NAG compiler
    Stop
    Do i=1, dim1
!      Call lubksb(temparr, dim1, dim1, rowperms, unitarr(1:dim1, i))
      invarr(1:dim1, i) = unitarr(1:dim1, i)
    End Do
  End Subroutine invertarray

  !-------------------------------------------------------------------------
  ! Multyplies two matrices and returns a resulting matrix
  !-------------------------------------------------------------------------

  Function matrix_mult(m1,m2) result(mout)
  Real(kind=RDbl), Dimension(:,:), Intent(In)         :: m1,m2
  Real(kind=RDbl), Dimension(size(m1,1),size(m2,2))   :: mout
  Integer            :: j, k, i, n, n11, n12, n21, n22
    
    !** Check if the matrices are linked)
    If(size(m1,2)/=size(m2,1)) Then
    Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
        ' Matrices are not linked'
        print*, size(m1,2), size(m2,1)
    Stop
    End If

    n11=size(m1,1)
    n12=size(m1,2)
    n21=size(m2,1)
    n22=size(m2,2)

    mout = 0.0_RDbl
    
    Do i=1, n11
      Do j=1, n22
          Do k=1, n12
             mout(i,j) = mout(i,j)+m1(i, k)*m2(k,j)
          End Do
       End Do
     End Do

  Return
  End Function matrix_mult

  !-------------------------------------------------------------------------
  ! Takes an array of coordinates in system of coordinates 1, randomly rotates
  ! this system of coordinates and returns this array of coordinates in 
  ! the rotated system
  !-------------------------------------------------------------------------

  Subroutine matrix_rotate_random(coord_in) 
  Type(VecType), Dimension(:), Intent(InOut)       :: coord_in 
  Real(kind=Rdbl)                                  :: costheta, theta, psi, phi
  
  Real(kind=Rdbl), Dimension(3,3)               :: rotation_matrix
  Real(kind=Rdbl), Dimension(3,size(coord_in))  :: matrix_in, matrix_out
  Real(kind=Rdbl)                               :: xc, yc, zc

  Integer                                       :: i

  !** Centering the geometric center of the molecule at 0, 0, 0

  xc = zero
  yc = zero
  zc = zero

  If(size(coord_in)<3) Return
  !! Rotation around geometric center of the molecule
  Do i = 1, size(coord_in)
  xc = xc + coord_in(i)%comp(1)
  yc = yc + coord_in(i)%comp(2)
  zc = zc + coord_in(i)%comp(3)
  End Do

  xc = xc/size(coord_in)
  yc = yc/size(coord_in)
  zc = zc/size(coord_in)

  Do i = 1, size(coord_in)
  coord_in(i)%comp(1) = coord_in(i)%comp(1) - xc
  coord_in(i)%comp(2) = coord_in(i)%comp(2) - yc
  coord_in(i)%comp(3) = coord_in(i)%comp(3) - zc
  End Do

  costheta = 1.0 - rranf() * 2.0
  theta = Acos(costheta) 
  phi = pi-rranf()*twopi
  psi = pi-rranf()*twopi

  Do i = 1, size(coord_in)
  matrix_in(1,i) = coord_in(i)%comp(1)
  matrix_in(2,i) = coord_in(i)%comp(2)
  matrix_in(3,i) = coord_in(i)%comp(3)
  End Do

  rotation_matrix(1,1) =   dcos(psi)*dcos(phi) - dcos(theta)*dsin(phi)*dsin(psi)
  rotation_matrix(1,2) =   dcos(psi)*dsin(phi) + dcos(theta)*dcos(phi)*dsin(psi)
  rotation_matrix(1,3) =   dsin(psi)*dsin(theta)
  rotation_matrix(2,1) =  -dsin(psi)*dcos(phi) - dcos(theta)*dsin(phi)*dcos(psi)
  rotation_matrix(2,2) =  -dsin(psi)*dsin(phi) + dcos(theta)*dcos(phi)*dcos(psi)
  rotation_matrix(2,3) =   dcos(psi)*dsin(theta)
  rotation_matrix(3,1) =   dsin(theta)*dsin(phi)
  rotation_matrix(3,2) =  -dsin(theta)*dcos(phi)
  rotation_matrix(3,3) =   dcos(theta)

  matrix_out = matrix_mult(rotation_matrix, matrix_in)
  
  Do i = 1, size(coord_in)
  coord_in(i)%comp(1) = matrix_out(1,i)
  coord_in(i)%comp(2) = matrix_out(2,i)
  coord_in(i)%comp(3) = matrix_out(3,i)
  End Do

  Do i = 1, size(coord_in)
  coord_in(i)%comp(1) = coord_in(i)%comp(1) + xc
  coord_in(i)%comp(2) = coord_in(i)%comp(2) + yc
  coord_in(i)%comp(3) = coord_in(i)%comp(3) + zc
  End Do

  
  Return
  End Subroutine matrix_rotate_random


  
End Module matrix

!------------------------------------------------------------------------
! Notes on coordinate rotation : Marion and Thornton gives the matrixes
! for rotation of coordinate systems.  But here we are more bothered
! about how to rotate a vector.  To rotate a coord system by angle phi
! around z axis in a-clockwise we use matrix L_phi.  Lets call the fixed
! (lab) coord system values as x and the values in the new coord system
! as x'. Then
! 
! x' = L_phi x
! 
! Where L_phi =  cos(phi)  sin(phi) 0
!               -sin(phi)  cos(phi) 0
!                0         0        1
!
! similarly for rotation by angle theta around x in a-clockwise
! 
! x' = L_theta x 
! 
! Where L_theta= 1         0           0
!                0         cos(theta)  sin(theta) 
!                0        -sin(theta)  cos(theta)
!
! 
! similarly for rotation by angle psi around z in a-clockwise  
!                
! x' = L_psi x 
! 
! Where L_psi =  cos(psi)  sin(psi) 0
!               -sin(psi)  cos(psi) 0
!                0         0        1 
! note that L_phi and L_psi have same forms
! 
! If we do rotations by phi, theta, psi in that order then
!  x' = L_psi L_theta L_phi x
! [This is equation 12.65 in Marion and thornton] where x' is values in new 
!  coord and x is the original value.
! 

! Assume a rigid body is rotated by phi, theta, psi. Originally we
! have two frames. The lab(fixed frame) and the body frame. Let a point
! on the body be denoted by x (in the fixed frame). Originaly these
! systems coincide. So x=x' in numerical value.  
! 
! Now we rotated the body frame by phi theta psi. But the value of x'
! does not change because the body is atteched the body frame. But he
! values of that point in the fixed co-ord has changed. Lets call the
! new coords as x1 and x1'
! 
! x1' = L_psi L_theta L_phi x1
! 
! we know x1' = x' =x
! 
! x = Lpsi L_theta L_phi x1
! 
! 
! where x1 is the value of the point after the body is rotated [in fixed frame]
!       x  is the value of the point before rotation [in fixed frame]
! 
! To start with we know x but our aim is to find x1
! 
! so,
! 
! x1=  L_phi^-1 L_theta^-1 L_psi^-1 x
! 
! So if we multiply a vector by the matrix given 12.66 then we
! ( matrix_f2beular)
! effectively rotate the body by -psi, -theta, -phi in that order. You
! might notice that phi, psi are interchangeable.
! ( Or rotate by phi, theta, psi in clockwise direction)
! 
! 1) To rotate a vector by phi, theta, psi :
!    use matrix_b2feular (L_phi^-1 L_theta^-1 L_psi^-1)
!  
! 2) To rotate a vector by dpsi around unit vector pointing along phi, theta
!    multiply by  - [L(-phi) L(-theta)] * [L(phi) L(theta) L(-dpsi)]
!                 - matrix_rotAroundVec(phi, theta, dpsi)             
!--------------------------------------------------------------------------
