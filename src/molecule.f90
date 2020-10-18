!--------------------------------------------------------
! This module defines molecular data type
!--------------------------------------------------------

 Module molecule

  Use defaults, Only: RDbl, strLen, lstrLen, one, zero
  Use vector, Only: VecType
  Implicit None
  Save

  Private
  Public :: MoleculeType


        Type MoleculeType
        Character(len=strLen)                         :: MoleculeName ! Molecule(or species name)
        Integer                                       :: natoms       ! Number of atoms
        Real(kind=Rdbl)                               :: molmass      ! Molecular mass
        Character(len=strLen), Dimension(:), Pointer  :: AtomName     ! Labels of individual atoms
        Type(VecType), Dimension(:), Pointer          :: coord        ! Atom coodinates within molecule
        Real(kind=Rdbl), Dimension(:), Pointer        :: sigma        ! Sigma
        Real(kind=Rdbl), Dimension(:), Pointer        :: eps          ! Epsilon
        Real(kind=Rdbl), Dimension(:), Pointer        :: q            ! Charge
        Real(kind=Rdbl), Dimension(:), Pointer        :: mass         ! Charge
        End Type MoleculeType

 End Module molecule
    
    
    
    
    
