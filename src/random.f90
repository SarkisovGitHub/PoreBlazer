!-------------------------------------------------------------------------
! This module contains the random number generator and related subroutines
! The main generator is expexted to give values from 0.0000.. to 
! 0.9999..
! There are two random number generators coded in here
! 1.) rranf_old : the one used till recently, but gives value 1.0000000000
!                 roughly after 800,000,000 calls. So not good for use in music
! 2.) rranf_NR  : "Taken from numerical recipes in fortran90: the art of 
!                 parallel  scientific computing" W.H Press et al. . 
!                 See pages 1141-1143 in chapter B7. 
!                 The book is online and free. The code is available at 
!                 "/fiesta/software/numerical.recipes/f90/recipes/ran.f90"
!                 It was modified to fit inot music's way of calling rranf
! Only one will be used for now. rranf_old can be later removed
! Shaji - Wed Oct 03 2001.
!-------------------------------------------------------------------------
Module random

  Use defaults, Only: RDbl, IDbl, Pi, lstrlen

  Implicit None
  Save

  Private
  Public :: random_init, rranf, random_gaussian, random_getiseed, & 
            random_gettrial, random_getnewseed

  !** Initialized to a negative value
  Integer       ::  iseed =-1

  !** K4B = 64bit Integer ( I guess!), These are specified in the ran.f90 
  !** program in the book. Safe to leave it unchanged. 
  !** variables IA, IM, IQ and IR has been rename by adding "CONST_" prefix
  Integer, Parameter :: K4B=Selected_int_kind(9)
  Integer(K4B), Parameter :: CONST_IA=16807,CONST_IM=2147483647
  Integer(K4B), Parameter :: CONST_IQ=127773, CONST_IR=2836

  !** these are the changing seeds for the new generator
  Integer(K4B)  :: ix=-1, iy=-1

  Interface rranf
!!$    Module Procedure rranf_old !** Only one of these can be used at a time
    Module Procedure rranf_NR
  End Interface

Contains

  !---------------------------------------------
  ! Set the seed for the random no. generator
  !---------------------------------------------
  Subroutine random_init(seed)
    Integer   :: seed
    iseed = seed
  End Subroutine random_init

  !-----------------------------------------------------------------
  ! Random Number Generator                                         
  ! Uniform pseudo-random number generator; mixed congruential      
  ! Type. The Function value is a Real number on the half-Open      
  ! interval [0,1) (i.e., including 0, but not 1).                  
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function rranf_old()
    Implicit None
    Integer(kind=IDbl) :: ib, ia, ibc, ida, isum, iff, ie, ix, iy, ix2, ix1
    
    ib   = ISeed/65536
    ia   = ISeed - ib*65536
    ibc  = ib*63253
    ida  = ia*24301
    isum = ibc - 2147483647 + ida
    
    If( isum .Le. 0 ) Then
      isum= isum + 2147483647
    Else
      isum= isum - 1
    End If
    iff= isum/32768
    
    ie= isum - iff*32768
    ix= ie + ia
    iy= 453816691 - 2283*ia
    ix2= ix/32768
    ix1= ix - 32768*ix2
    ISeed= ix1*65536 - 2147483647 + iy
    
    If( ISeed .Le. 0 ) Then
      ISeed= ISeed + 2147483647
    Else
      ISeed= ISeed - 1
    End If
    
    rranf_old= ISeed*4.65661287524579692410575D-10
    
    !   Note on the final Assignment of random number.
    !   
    !   The value of ISeed upon entering routine:  1349134852
    !   The value of ISeed upon leaving routine :  2147483591
    !   
    !   rranf= ISeed*4.65661287308E-10
    !        = 1.0000000000000
    !    
    !   rranf= ISeed*4.65661287308D-10
    !        = 0.99999997345787
    !   
    !   rranf= ISeed/Dble(2147483647)
    !        = ISeed*4.6566128752458D-10
    !        = 0.99999997392297
    !   
    !   rranf= ISeed*4.65661287524579692410575D-10
    !        = 0.99999997392297
    
    Return
  End Function rranf_old

  !-----------------------------------------------------------------
  ! Uniform Random Number Generator
  ! see header of the Module for more details
  ! The Function value is a Real number on the Open      
  ! interval (0,1) (0 and 1 not included)
  !-----------------------------------------------------------------
  Real(kind=RDbl) Function rranf_NR()
    Implicit None
    Real(kind=RDbl), Save :: am
    Integer(K4B)  :: k
    Logical, Save :: first_time=.True.

    !** Initialize
    If (first_time) Then
      first_time=.False.
      am=Nearest(1.0,-1.0)/CONST_IM
      iy=Ior(Ieor(888889999,Abs(iseed)),1)
      ix=Ieor(777755555,Abs(iseed))
      iseed=Abs(iseed)+1
    End If

    !** The Marsaglia shift sequence with period of 2**32 -1
    ix=Ieor(ix,Ishft(ix,13))
    ix=Ieor(ix,Ishft(ix,-17))
    ix=Ieor(ix,Ishft(ix,5))

    !** The Park-Miller sequence by schrages's method. Period=2**31-2
    k=iy/CONST_IQ
    iy=CONST_IA*(iy-k*CONST_IQ)-CONST_IR*k
    If (iy < 0) iy=iy+CONST_IM

    !** Combine both sequences (to improve randomness, I guess)
    !** Note that the combined sequence has period of 2**64
    rranf_NR=am*Ior(Iand(CONST_IM,Ieor(ix,iy)),1)

    Return

  End Function rranf_NR

  !----------------------------------------------------------
  ! Chooses a trial according to its cumulative probability
  !----------------------------------------------------------
  Integer Function random_gettrial(orig_low,orig_high,cum_prob)
    Integer, Intent(In)         :: orig_low, orig_high
    Real(kind=RDbl), Dimension(:), Intent(In) :: cum_prob

    Integer                     :: low,high,mid,i 
    Real(kind=RDbl)             :: key

    low = orig_low
    high = orig_high
    If (low >  high) Then
      Write(0,'(1x,2a,i4,a, i4)') __FILE__," : ",__LINE__, &
          "Could not find the index "
      Stop
    Else If (low == high) Then
      random_gettrial = low
      Return
    End If

    key  = rranf()
    mid = (low + high)/2
    If(cum_prob(low) > key) Then
      random_gettrial = low
      Return
    Else
       Do i=low,high-1
        If ( key > cum_prob(i) .AND. key < cum_prob(i+1)) Then
           random_gettrial = i + 1
            Return
        End If
       End Do
    End IF
 
  End Function random_gettrial

  !----------------------------------------------------------------------------
  ! Generates a random number according the Gaussian distribution
  !* This function uses the Box-Muller method.
  !* Reference: "Monte Carlo Methods. Volume 1:Basics" 
  !*             by M. H. Kalos, P. A. Whitlock; pp. (45-48).
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function random_gaussian(ave,sigma)
    Real(Kind=RDbl), Intent(In) :: ave
    Real(Kind=RDbl), Intent(In) :: sigma
    Real(kind=RDbl), Parameter :: twopi = 2.0_RDbl*pi
    Real(Kind=RDbl), save :: rnum1, rnum2, zeta1, zeta2
    Integer, Save :: step = 1

    If (step == 1) Then
      step = 0
      rnum1 = rranf()
      rnum2 = rranf()
      If (rnum1 < 1.0E-10_RDbl) rnum1 = 1.0E-10_RDbl
      If (rnum2 < 1.0E-10_RDbl) rnum2 = 1.0E-10_RDbl
      zeta1 = Sqrt(-2.0_RDbl*Log(rnum1))*Cos(twopi*rnum2)
      zeta2 = Sqrt(-2.0_RDbl*Log(rnum1))*Sin(twopi*rnum2)
      random_gaussian = ave+sigma*zeta1
    Else If (step == 0) Then
      random_gaussian = ave+sigma*zeta2
      step = 1
    End If
  End Function random_gaussian

  !------------------------------------------------
  ! Returns the iseed
  !------------------------------------------------
  Integer Function random_getiseed()
    random_getiseed = iseed
  End Function random_getiseed

  !------------------------------------------------
  ! Returns the changing seed for the new generator
  !------------------------------------------------
  Function random_getnewseed()
    Character(len=lstrlen)   :: random_getnewseed

    Write(random_getnewseed,'(2i12)') ix,iy

  End Function random_getnewseed

  !------------------------------------------------
  ! Print out the current seed
  !------------------------------------------------
  Subroutine random_display(optunit)
    Integer, Optional, Intent(in)  :: optunit

    Integer   :: unitno

    If (Present(optunit)) Then
      unitno = optunit
    Else
      unitno = 6
    End If

    Write(unitno, '(7x,a,i10)') "The random seed is: ", iseed
  End Subroutine random_display

End Module random
