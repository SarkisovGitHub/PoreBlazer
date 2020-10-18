!------------------------------------------------------------------------------
! This module contains upper limits on arrays difficult to allocate on the fly.
! It also contains the filename tag definitions
!------------------------------------------------------------------------------
Module defaults
  Implicit None
  Save

  Public 

  Integer :: pcalls

  !-------------------------------------------------
  ! Input file comment character
  !-------------------------------------------------
  Character(len=1), Parameter :: COMMENT_CHAR = "#"

  !--------------------------------------------
  ! Type Definitions
  !--------------------------------------------
  Integer, Parameter   :: RDbl = Selected_real_kind(10, 50)
  Integer, Parameter   :: RSgl = Selected_real_kind(5,10)
  Integer, Parameter   :: IDbl = Selected_int_kind(10)
  Integer, Parameter   :: ISgl = Selected_int_kind(4)
  Integer, Parameter   :: strLen = 48 ! Default character string length
  Integer, Parameter   :: lstrLen = 2*strLen ! Long string length
  Integer, Parameter   :: xlstrLen = 3*strLen ! Extra long string length
  Real(kind=RDbl), Parameter  :: zeroTolerance = 1.0e-6_Rdbl ! For comparing
                                                             ! Reals

  !--------------------------------------------
  ! Mathematical Constants
  !--------------------------------------------
  Real(kind=RDbl), Parameter  :: pi    = 3.1415927_RDbl
  Real(kind=RDbl), Parameter  :: twopi = 2.0_RDbl*pi
  Real(kind=RDbl), Parameter  :: degTorad = pi/180.0_RDbl
  Real(kind=RDbl), Parameter  :: radTodeg = 1.0_RDbl/degTorad
  Real(kind=RDbl), Parameter  :: zero = 0.0_RDbl
  Real(kind=Rdbl), Parameter  :: one = 1.0_RDbl

  !--------------------------------------------
  ! Physical Constants
  ! Note that the internal units are as follows:
  ! length = Ang
  ! time   = ps
  ! mass   = amu
  ! temperature = K
  ! velocity = Ang/ps
  ! accelerations = Ang/ps^2
  ! energy = kcal (when calculated)
  ! energy = kJ (when stored and reported)
  !--------------------------------------------
  Character(len=strLen), Parameter :: lengthUnit = "Ang"
  Character(len=strLen), Parameter :: timeUnit = "ps"
  Character(len=strLen), Parameter :: massUnit = "amu"
  Character(len=strLen), Parameter :: TUnit = "K"
  Character(len=strLen), Parameter :: nrgUnit = "kJ/mol"
  Character(len=strLen), Parameter :: velUnit = Trim(lengthUnit)//"/"//Trim(timeUnit)

  Real(kind=RDbl), Parameter  :: Echarge = 1.60E-19 
  ! electron chrage in coulombs

  Real(kind=RDbl), Parameter  :: Nav = 6.0221367E+23_RDbl   ! Avogadro's Number
  Real(kind=RDbl), Parameter  :: calToj = 4.184_RDbl  
  Real(kind=RDbl), Parameter  :: bohrToAng  = 0.529177_RDbl
  Real(kind=RDbl), Parameter  :: HTokJmol = 4.3598E-21_RDbl*Nav !Hartrees->kJ/mol
  Real(kind=RDbl), Parameter  :: HTokcalmol = HTokJmol/calToj   !H->kcal/mol
  Real(kind=RDbl), Parameter  :: eVTokcalmol = 23.0605_RDbl  

  Real(kind=RDbl), Parameter :: Rgas    = 8.314_Rdbl ! J/mol/K
  Real(kind=RDbl), Parameter :: hplanck = 6.626E-34_RDbl   ! J.s
  Real(kind=RDbl), Parameter :: scalepe = 4.184_RDbl ! (k)cal -> (k)J Unit conv
  Real(kind=RDbl), Parameter :: scalef  = 418.4_RDbl 
  ! Actual conversion : kcal -> amu. ang^2 / ps^2      (energy)
  !  OR                 kcal/ang -> amu. ang/ps^2      (force)
  Real(kind=RDbl), Parameter :: scaleke = 0.01_RDbl   !(Ang/ps)^2*amu -> KJ/mol
  Real(kind=RDbl), Parameter :: jmolec_kb    = Rgas/Nav       ! J/molecule/K
  Real(kind=RDbl), Parameter :: kjmole_kb    = 8.314E-3_RDbl ! R, KJ/mole/K
  Real(kind=RDbl), Parameter :: jmole_kb     = 8.314_RDbl         ! J/mole/K
  Real(kind=RDbl), Parameter :: kcalmole_kb  =  kjmole_kb/calToj! kcal/mole/K
  Real(kind=RDbl), Parameter :: calmole_kb   = jmole_kb/calToj         ! cal/mole/K



  !** Converts from e^2/Ang to kcal/mol. Remember that the charges'
  !** values are in elementry charge units e = 1.602177E-19 C. 
  !** Usual derivations (this one included) drop the 4*Pi*eps0 term
  !** from the equations, where eps0 = 8.85419E-12 C^2/(J*m).
  !** This yields units of e^2/Ang. The conversion factor, therefore, is 
  !  Real(Kind=RDbl), Parameter :: tokcal = 332.064_RDbl
  Real(Kind=RDbl), Parameter :: e2kcal = 331.5 
  ! we used 331.5 in other parts of music, need to verify this
  !** conversion constant(from atomic unit to kcal/mol) 
  !** = 9.0*10^9*(1.6*10^-19)^2/10^(-10) * 6.02*10^23 /4.184 = 331 kcal/mol
  
  !------------------------------------------------------------------
  ! Just a flag that can make debugging easier
  !------------------------------------------------------------------
  Logical            :: dbgflag = .False.


End Module defaults




