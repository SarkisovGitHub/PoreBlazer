!------------------------------------------------------------------------------
! This module contains various generic matrix operations, such as L-U 
! decomposition, back substitution, SVD decomposition, and more.
! See Numerical Recipies for more information.
!------------------------------------------------------------------------------
Module matrixops

  Use defaults, Only: RDbl

  Implicit None

  Private
  Public :: matrixops_ludcmp, matrixops_lubksb

  Interface matrixops_ludcmp
    Module Procedure matrixops_ludcmpReal
    Module Procedure matrixops_ludcmpInt
  End Interface

  Interface matrixops_lubksb
    Module Procedure matrixops_lubksbReal
    Module Procedure matrixops_lubksbInt
  End Interface

Contains

  !----------------------------------------------------------------------------
  ! Performs a L-U decomposition on matrix A of size N. In this case, 
  ! the INDX parameter is a real
  !----------------------------------------------------------------------------
  Subroutine matrixops_ludcmpReal(A,N,NP,INDX,D)
    Integer, Parameter :: NMAX=100
    Real(Kind=RDbl) :: TINY=1.0E-20
    Real(Kind=RDbl) :: D, AAMAX, SUM, DUM
    Integer         :: NP, N
    Real(Kind=RDbl), Dimension(:,:) :: A
    Real(Kind=RDbl), Dimension(NMAX)  :: VV
    Real(Kind=RDbl), Dimension(:)     :: INDX
    Integer    :: I, J, K, IMAX
    
    
    D=1.0_RDbl
    
    Do I=1,N
      AAMAX=0.0_RDbl
      Do J=1,N
        If (Abs(A(I,J)).Gt.AAMAX) AAMAX=Abs(A(I,J))
      End Do
      If (AAMAX.Eq.0.) Write(0,'(a,i5,a)') __FILE__,__LINE__, &
          ': Singular matrix in LUDCMP.'
      VV(I)=1./AAMAX
    End Do
    
    Do J=1,N
      If (J.Gt.1) Then
        Do I=1,J-1
          SUM=A(I,J)
          If (I.Gt.1)Then
            Do K=1,I-1
              SUM=SUM-A(I,K)*A(K,J)
            End Do
            A(I,J)=SUM
          End If
        End Do
      End If
      AAMAX=0.
      Do I=J,N
        SUM=A(I,J)
        If (J.Gt.1)Then
          Do K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          End Do
          A(I,J)=SUM
        End If
        DUM=VV(I)*Abs(SUM)
        If (DUM.Ge.AAMAX) Then
          IMAX=I
          AAMAX=DUM
        End If
      End Do
      If (J.Ne.IMAX)Then
        Do K=1,N
          DUM=A(IMAX,K)
          A(IMAX,K)=A(J,K)
          A(J,K)=DUM
        End Do
        D=-D
        VV(IMAX)=VV(J)
      End If
      INDX(J)=Real(IMAX,Kind=RDbl)
      If(J.Ne.N)Then
        If(A(J,J).Eq.0.)A(J,J)=TINY
        DUM=1./A(J,J)
        Do I=J+1,N
          A(I,J)=A(I,J)*DUM
        End Do
      End If
    End Do
    
    If(A(N,N).Eq.0.)A(N,N)=TINY
    Return
  End Subroutine matrixops_ludcmpReal

  !----------------------------------------------------------------------------
  ! Performs a L-U decomposition on matrix A of size N. In this case, 
  ! the INDX parameter is an integer
  !----------------------------------------------------------------------------
  Subroutine matrixops_ludcmpInt(A,N,NP,INDX,D)
    Integer, Parameter :: NMAX=100
    Real(Kind=RDbl) :: D
    Integer         :: NP, N
    Real(Kind=RDbl), Dimension(:,:) :: A
    Integer, Dimension(:)     :: INDX
    Real(Kind=RDbl), Dimension(Size(INDX,1)) :: RINDX

    !** Convert to reals
    RINDX = Real(INDX,kind=RDbl)
    
    !** Call matrixops_ludcmpReal
    Call matrixops_ludcmpReal(A,N,NP,RINDX,D)
    
    !** Convert back to INT
    INDX = Int(RINDX)
  End Subroutine matrixops_ludcmpInt

  !----------------------------------------------------------------------------
  ! Subroutine for performing L-U backsubtitution
  !----------------------------------------------------------------------------
  Subroutine matrixops_lubksbReal(A,N,NP,INDX,B)
    Integer         ::    NP, N, I, II, LL, J
    Real(Kind=RDbl) ::    SUM
    Real(Kind=RDbl), Dimension(NP,NP) :: A
    Real(Kind=RDbl), Dimension(NP)    :: B
    Real(Kind=RDbl), Dimension(N)     :: INDX

    II=0
    Do I=1,N
      LL=Int(INDX(I))
      SUM=B(LL)
      B(LL)=B(I)
      If (II.Ne.0)Then
        Do J=II,I-1
          SUM=SUM-A(I,J)*B(J)
        End Do
      Else If (SUM.Ne.0.) Then
        II=I
      End If
      B(I)=SUM
    End Do

    Do I=N,1,-1
      SUM=B(I)
      If(I.Lt.N)Then
        Do J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        End Do
      End If
      B(I)=SUM/A(I,I)
    End Do

    Return
  End Subroutine matrixops_lubksbReal

  !----------------------------------------------------------------------------
  ! Subroutine for performing L-U backsubtitution
  !----------------------------------------------------------------------------
  Subroutine matrixops_lubksbInt(A,N,NP,INDX,B)
    Integer         ::    NP, N
    Real(Kind=RDbl), Dimension(NP,NP) :: A
    Real(Kind=RDbl), Dimension(NP)    :: B
    Integer, Dimension(N)     :: INDX

    Real(kind=Rdbl), Dimension(N) :: RINDX
    
    RINDX = Real(INDX,Kind=RDbl)
    Call matrixops_lubksbReal(A,N,NP,RINDX,B)
    INDX = Int(RINDX)
  End Subroutine matrixops_lubksbInt

End Module matrixops


