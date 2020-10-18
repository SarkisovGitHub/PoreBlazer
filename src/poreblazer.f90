
!--------------------------------------------------------------------------------------
! This is the main program for characterization of porous solids
!
! Poreblazer 4.0
!
! Lev Sarkisov 2020
!--------------------------------------------------------------------------------------

module parameters
        use fundcell, only: Fundamental_Cell, fundcell_getvolume
        type(Fundamental_Cell) :: fcell                  ! contains information about the strucutre of the simulation cell
        real*8  :: hicut, hicut2                         ! cutoff potential distance (A) and its squared value
        real*8  :: coeff_surface, coeff_surface2         ! additional coefficient to switch between hard sphere (1) or potential minima (1.122) accessible surface
        real*8  :: temp                                  ! temperature (K) required for Helium porosimetry calculations
        integer :: nsample                               ! number of tests for surface area
        integer :: iseed                                 ! random seed number
        integer :: vis_option                            ! visualization options
        character*20 :: property                         ! differentiates between total and accessible properties
end module parameters

module atoms
        integer*4, allocatable :: atype(:)               ! this array is natoms long, for each atom i containing the number corresponding to the type of atom
        character*10, allocatable :: aname(:)            ! this array is ntypes long, for each atom i containing the name corresponding to the type of atom
        real*8, allocatable :: asigma(:), asigma2(:)     ! these arrays are ntypes long, for each atom type i containing the collision diameter for that type and squared value of it
        real*8, allocatable :: aeps(:)                   ! this array is ntypes long, for each atom type i containing the LJ epsilon, K
        real*8, allocatable :: amass(:)                  ! this array is ntypes long, for each atom type i containing molar mass
        real*8, allocatable :: asigma2_he(:), aeps_he(:) ! these arrays contain mixed values for sigma and epsilon for atom type i interacting with atom of helium
        real*8, allocatable :: asigma_n(:), asigma2_n(:) ! these arrays contain mixed values for sigma and squared value of it for atom type i interacting with atom of nitrogen
end module atoms

module adsorbent
        use  vector, only:       vectype
        integer*4 :: natypes, natoms                                       ! number of atom types, and atoms in adsorbent structure
        real*8, allocatable :: coords(:, :), coords_temp(:,:)              ! coordinates of the adsorbent (natoms long array)
        type(vectype), allocatable :: matvec(:)                            ! 
        character*10, allocatable :: adsname(:)                            ! names of all atoms in the input coordinates file
end module adsorbent

module lattice
        real*8                   :: cube_size                                         ! size of the lattice cubelet
        integer                  :: ncubesx, ncubesy, ncubesz                         ! number of lattice cubes in each direction and other parameters
        integer                  :: ntot, spanning                                    ! total number of lattice cubes,  logical variable for percolation
        integer                  :: nhe_cubes, nn_cubes, ng_cubes                     ! number of lattice cubes accessible to helium, nitrogen and a point, respectively
        integer*2, allocatable   :: lattice_space(:,:,:)                              ! lattice cubes: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible)
        integer*2, allocatable   :: lattice_temp(:,:,:)                               ! temp array for various uses: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible)
        integer, allocatable     :: cl_summary(:,:)                                   ! cluster summary: stores information about percolating clusters and the number of percolated dimensions
        integer, allocatable     :: lattice_index(:,:)                                ! links the cubelet number and its i, j, k side indicies
        real*8, allocatable      :: lattice_rdist2(:,:,:)                             ! squared distances between the center of the cubelet and nearest atom of the structure
        integer*2, allocatable   :: lattice_space_he(:,:,:), lattice_space_n(:,:,:)   ! lattice cubes: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible) for helim and nitrogen
        integer, allocatable     :: he_cubes(:), n_cubes(:), g_cubes(:)               ! list of cubes available to helium, nitrogen and a point
        real*8, allocatable      :: lattice_lj_he(:)                                  ! for each cube i, summary of the LJ interaction between the He atom in the center of the cube and the structure
        real*8,  allocatable     :: PA1(:)                                            ! array used for percolation analysis, based on site occupation
        integer, allocatable     :: PA2(:), PA3(:), PA4(:)                            ! arrays used for percolation analysis, based on site occupation
end module lattice


module distributions
        real*4, allocatable :: psd_cumul(:), psd(:)                     ! cumulative pore volume Vp function and negative derivative of that function, i.e. pore size distribution
        real*4              :: binsize, maxpore                                       ! size of the bin for PSD analysis (A), maximum anticipated pore diameter (A)
        integer             :: nbins                                                  ! number of bins for PSD analysis
end module distributions

module results
        real*8 :: volume, smass, density, rho, pore_v_he, pore_geom, stotal, pore_lim, pore_max, free_volume, ffv
        integer :: sys_perc
end module results

program poreblazer
use parameters
use results
implicit none
interface initialize
    subroutine initialize(arg)
        character(len=*), optional :: arg
    end subroutine initialize
end interface initialize


character(len=256) :: arg

write(*,*) "!-------------------------------------------------------!"
write(*,*) "! Welcome to PB4.0                                      !"
write(*,*) "! Developed by: Lev Sarkisov, 2012/20                   !"
write(*,*) "! Cite: Sarkisov, Harrison, Molecular Simulation, 2011  !"
write(*,*) "!-------------------------------------------------------!"
write(*,*)

write(*,*) "!-------------------------------------------------------!"
write(*,*) "! The program first goes through several preliminary    !"
write(*,*) "! initializations, associated with reading parameters.  !"
write(*,*) "!-------------------------------------------------------!"
write(*,*) 

write(*,*) "Step 1/10"

if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    call initialize(arg)
else
    call initialize
end if

    volume = fundcell_getvolume(fcell)
    write(*,'(a,f12.3)') " System volume in A^3:    ",  volume
    write(*,'(a,f12.3)') " System mass, g/mol:      ",  smass
    write(*,'(a,f12.3)') " System density, g/cm^3:  ",  rho
    write(*,*)
    write(*,*) "-------------------------------------------------------"
    write(*,*) "If calculated density is  different from the  reported "
    write(*,*) "crystallographic one,  make sure you remove unresolved "
    write(*,*) "atoms from the adsorbent structure file"
    write(*,*) "------------------------------------------------------"
    write(*,*)
    write(100,'(a,f12.3)') "V_A^3 ",  volume
    write(100,'(a,f12.3)') "M _g/mol ",  smass
    write(100,'(a,f12.3)') "RHO_g/cm^3 ",  rho
 
write(*,*) "Step 2/10"
call lattice_calculations
call limiting_diameter           ! calculation of the pore limiting diamater

!--------------------------------------------------------------------------------------
!      Now we sequentially call several routines where the actual
!      analysis of the structure is performed
!--------------------------------------------------------------------------------------

write(*,*)
write(*,*) "!=======================================================!"
write(*,*)
write(*,*) "!--------------------------------------------------------!"
write(*,*) "! Phase 1: Total properties                              !"
write(*,*) "!--------------------------------------------------------!"
write(*,*)


property = "Total "
write(*,*) "Step 3/10"
call surface_area                ! surface area
write(*,*) "Step 4/10"
call pore_distribution           ! geometric pore size distribution
write(*,*) "Step 5/10"
call volumes                     ! various volumes

write(*,*)
write(*,*) "!=======================================================!"
write(*,*)
write(*,*) "!--------------------------------------------------------!"
write(*,*) "! Phase 2: Network-accessible properties                 !"
write(*,*) "!--------------------------------------------------------!"
write(*,*)

property = "Network-accessible "
write(*,*) "Step 6/10"
call geom_lattice                ! returns lattice accessible to point probe
write(*,*) "Step 7/10"
call helium_lattice              ! returns lattice accessible to helium
write(*,*) "Step 8/10"
call nitrogen_lattice            ! returns lattice accessible to nitrogen


write(*,*) "Step 9/10"
call surface_area                ! accessible surface area
write(*,*) "Step 10/10"
call pore_distribution           ! accessible geometric pore size distribution
write(*,*) "Step 11/10"
call volumes                     ! accessible volumes


!call finalize

write(*,*) "!-------------------------------------------------------!"
write(*,*) "! Simulation complete                                   !"
write(*,*) "!-------------------------------------------------------!"


end program

!--------------------------------------------------------------------------------------
! THE END OF THE PROGRAM
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
! SUBROUTINES
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
! Initialization: the main purpose of this subroutine is to read in various types of
! simulation parameters and initialize arrays
!--------------------------------------------------------------------------------------

subroutine initialize(filename)

    use parameters
    use atoms
    use adsorbent
    use lattice
    use distributions

    use defaults, only: rdbl
    use fundcell, only: fundcell_init, fundamental_cell,  fundcell_slant, fundcell_unslant
    use  vector, only:  vectype
    use random, only:   random_init
    use results

    implicit none
    character(len=*), optional :: filename                                              ! input filename, defaults stdin
    integer*4                             :: i, j, k, l                                 ! cycle indicies
    logical                               :: check                                      ! logical variable, required for input analysis
    character(256)                        :: filename1, filename2, filename3            ! data file names
    real*8                                :: sigma_he, eps_he, sigma_n                  ! sigma (A) and epsilon (K) of helium; sigma (A) of nitrogen atom
    type(vectype)                         ::  atvec1

    ! Initialization cycle

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Proceed with initialization and data input            !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    if (present(filename)) then
        open(4, file=trim(adjustl(filename)), status='old')
        filename3 = "defaults.dat"
        read(4,'(a)') filename2 !  name of the xyz coordinates file
        read(4,*) fcell%ell(1), fcell%ell(2), fcell%ell(3)    ! non-orthorhombic dimensions of the simulation cell
        read(4,*) fcell%anglebc, fcell%angleac, fcell%angleab ! cell angles  (RDbl)  (alpha,beta,gamma)
        close(4)
    else
        filename3 = "defaults.dat"
        read(*,*) filename2  !  name of the xyz coordinates file
        read(*,*) fcell%ell(1), fcell%ell(2), fcell%ell(3)    ! non-orthorhombic dimensions of the simulation cell
        read(*,*) fcell%anglebc, fcell%angleac, fcell%angleab ! cell angles  (RDbl)  (alpha,beta,gamma)
    end if

    open(2, file=Trim(adjustl(filename2)), status='old')
    open(3, file=Trim(adjustl(filename3)), status='old')

    read(3,*) filename1  !  name of the atom types file
    read(3,*) sigma_he, eps_he, temp, hicut 
    read(3,*) sigma_n
    read(3,*) nsample

    open(1, file=Trim(adjustl(filename1)), status='old')

    hicut2 = hicut*hicut

    coeff_surface = 1.122 ! hard-coded now
    coeff_surface2 = coeff_surface*coeff_surface

    ! Here we read in the number of atom types and their collision diameters

    read(1,*) natypes !  number of atom types
    allocate (aname(natypes),asigma(natypes),aeps(natypes), &
    amass(natypes), asigma2(natypes), asigma2_he(natypes), aeps_he(natypes), asigma_n(natypes), asigma2_n(natypes))

    do i=1, natypes  ! reading in characteristics of atoms
        read(1,*) aname(i), asigma(i), aeps(i), amass(i)
        asigma2(i)    = asigma(i)*asigma(i)
        asigma2_he(i) = (0.5*(asigma(i) + sigma_he))**2.0
        aeps_he(i)    = sqrt(aeps(i)*eps_he)
        asigma_n(i)   = (0.5*(asigma(i) + sigma_n))
        asigma2_n(i)  = (0.5*(asigma(i) + sigma_n))**2.0
    end do

    ! Now we initialize the adsrobent structure

    read(2,*) natoms !  number of atoms in the strucutre
    allocate (adsname(natoms), coords(3, natoms), coords_temp(3, natoms), atype(natoms), matvec(natoms))

    ! Here we read in x, y, z, coordinates of the atoms from the input file
    ! The format of the file is xyz  file

    read(2,*)

    smass = 0.0 ! mass of the unit cell g/mol

    do i=1, natoms
        read(2,*) adsname(i), coords(1, i), coords(2, i), coords(3, i)

        check=.False.

        do j=1, natypes
            if(aname(j)==trim(adjustl(adsname(i)))) then
                atype(i)=j
                smass = smass + amass(j)
                check=.True.
                exit
            end if
        end do
        if(.Not.check) then
            print*, 'No such atom type:   ', trim(adjustl(adsname(i))), aname(i)
            stop
        end if
    end do

    ! Parameters of the simulations



    read(3,*) cube_size            ! lattice cube_size size in A and number of rotations per molecule

    ! Parameters of the simulation box

    fcell%origin%comp=0.0                                 ! coordinates of the origin of the cell (0,0,0)

    call fundcell_init(fcell)

    do i=1, natoms
    matvec(i)%comp = coords(:,i)
    end do
    rho = (smass/fcell%volume)*1.660539

    fcell%eff(1) = fcell%ell(1)
    fcell%eff(2) = fcell%ell(2)
    fcell%eff(3) = fcell%ell(3)
  
     do i=1, natoms
     matvec(i) = fundcell_slant(fcell, matvec(i))
     end do

    do i=1, natoms
    coords_temp(:,i) =  matvec(i)%comp 
    end do

    fcell%mx =  minval(coords_temp(1, :))
    fcell%my =  minval(coords_temp(2, :))
    fcell%mz =  minval(coords_temp(3, :))


    coords_temp(1, :) = coords_temp(1, :) - fcell%mx
    coords_temp(2, :) = coords_temp(2, :) - fcell%my
    coords_temp(3, :) = coords_temp(3, :) - fcell%mz

   do i=1, natoms
   matvec(i)%comp = coords_temp(:,i)
   end do

   atvec1%comp(1) = fcell%mx
   atvec1%comp(2) = fcell%my
   atvec1%comp(3) = fcell%mz

   atvec1 = fundcell_unslant(fcell, atvec1)

    ! Adjusting the position of the system in case it is seriously outside the simulation cell

    coords(1, :) = coords(1, :) - atvec1%comp(1)
    coords(2, :) = coords(2, :) - atvec1%comp(2)
    coords(3, :) = coords(3, :) - atvec1%comp(3)

    open(40, file = 'new.xyz')
    write(40,*) natoms     !  number of atoms in the strucutre

    write(40,*)

    do i=1, natoms
        write(40,'(a, 3f10.5)') adsname(i), coords(1, i), coords(2, i), coords(3, i)
    end do

    close(40)

    ncubesx = int(fcell%eff(1)/cube_size)                ! number of cubelets on x side
    cube_size = fcell%eff(1)/dble(ncubesx)               ! corrected size of the cubelet
    ncubesx = int(fcell%eff(1)/cube_size)                ! number of cubelets on x side
    ncubesy = int(fcell%eff(2)/cube_size)                ! number of cubelets on y side
    ncubesz = int(fcell%eff(3)/cube_size)                ! number of cubelets on z side

    ntot = ncubesx*ncubesy*ncubesz                       ! total number of lattice cubelets in the system
    allocate(lattice_space(ncubesx, ncubesy, ncubesz), lattice_temp(ncubesx, ncubesy, ncubesz), g_cubes(ntot))
    allocate(lattice_index(3, ntot))
    allocate(lattice_rdist2(ncubesx, ncubesy, ncubesz))
    allocate(lattice_space_he(ncubesx, ncubesy, ncubesz), he_cubes(ntot))
    allocate(lattice_space_n(ncubesx, ncubesy, ncubesz), n_cubes(ntot))
    allocate(lattice_lj_he(ntot))
    allocate(cl_summary(20,2))

    lattice_space = 0; lattice_temp = 0; g_cubes = 0
    lattice_index = 0
    lattice_rdist2 = 0.0
    lattice_space_he = 0;  he_cubes = 0
    lattice_space_n = 0;   n_cubes = 0
    nhe_cubes = 0;  nn_cubes = 0; ng_cubes = 0
    lattice_lj_he = 0.0d0
    cl_summary = 0

    read(3,*) maxpore, binsize    ! information related to pore size distrubution
    nbins = int(maxpore/binsize)
    binsize = maxpore/real(nbins)

    allocate(psd_cumul(nbins+100), psd(nbins+100))

    read(3,*) iseed               ! random seed

    if(iseed>0) then
        iseed = -iseed
    end if

    call random_init(iseed)


    read(3,*) vis_option

    close(1)
    close(2)
    close(3)

    open(100, file = "summary.dat", status = 'unknown')
    write(100, *) Trim(adjustl(filename2))
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Initialization complete                               !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine initialize

!--------------------------------------------------------------------------------------
! Finalize: formally deallocate all arrays
!--------------------------------------------------------------------------------------                

subroutine finalize
    use atoms
    use adsorbent
    use lattice
    use distributions

    implicit none

    ! These arrays are allocated in initialize
    deallocate(aname,asigma,aeps, amass, asigma2, asigma2_he, aeps_he, asigma_n, asigma2_n)
    deallocate(adsname, coords, atype, matvec)
    deallocate(lattice_space, lattice_temp, g_cubes)
    deallocate(lattice_index)
    deallocate(lattice_rdist2)
    deallocate(lattice_space_he, he_cubes)
    deallocate(lattice_space_n, n_cubes)
    deallocate(lattice_lj_he)
    deallocate(cl_summary)
    deallocate(psd_cumul, psd)

    ! These arrays are allocated in lattice_calculations
    deallocate(PA1, PA2, PA3, PA4)
end subroutine finalize

!--------------------------------------------------------------------------------------
!       Now we will perform a number of preliminary calculations:
!       - distance from each lattice cube center to all atoms, detecting the shortest distance
!       - lattice accessible to a helium atom and cluster accessible to a helium atom
!       - preliminary information for helium and geometric pore volume
!       - lattice accessible to a nitrogen atom
!--------------------------------------------------------------------------------------

subroutine lattice_calculations

    use parameters
    use atoms
    use adsorbent
    use lattice

    use defaults, only:     rdbl
    use fundcell, only:     fundamental_cell, fundcell_snglminimage
    use vector, only:       vectype
    use percolation, only:  percolation_calc, percolation_calc_simple

    implicit none

    character(20)                      :: filename1, filename2, filename3
    type(vectype)                      :: atvec1, atvec2, sepvec
    real*8                             :: rdist, rdist2, rdist6, rdist12, rdist2_ref, rdist_surface, rdist_surface_ref
    real*8                             :: sigma, sigma6, sigma12,  sig2_rdist2, lj_energy
    integer*4                          :: i, j, k, l, icount, nstep,amin
    logical                            :: overlap


    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting preliminary lattice calculations             !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    icount = 0

    do l=1, ncubesz                                    ! we go cubelet by cubelet
        do k=1, ncubesy
            do j=1, ncubesx

                icount = ((l-1) * ncubesx * ncubesy) + ((k-1) * ncubesx) + j                                ! count cubelets
                ! create a look-up table connecting the cubelet number with its indicies
                lattice_index(:, icount) = (/ j, k, l /)

                ! calculate the coordinates of the center of each cubelet
                atvec1%comp = (/dble(j-1), dble(k-1), dble(l-1)/) * cube_size + (0.5 * cube_size)

                overlap = .False.

                amin = 1
                rdist2_ref = huge(0.0d0)
                rdist_surface_ref = huge(0.0d0)

                do i=1, natoms                                     ! for each cubelet go through the whole set of atoms of the adsorbent structure

                    call fundcell_snglMinImage(fcell,atvec1,matvec(i),sepvec,rdist2)     ! calculate the distance between each atom and center of cubelet icount

                    ! First check if the point is "inside" a sphere
                    if(rdist2<0.25*asigma2(atype(i))) then             ! ignore the cubelet icount if it is "inside" an atom
                        overlap=.True.
                        exit
                    end if

                    if(rdist2<rdist2_ref) then                         ! in the next few lines we identify if the returned distance squared rdist2 is the smallest so far between
                        rdist2_ref = rdist2                            ! cubelet icount and an atom of the structure, without an overlap
                        amin = i
                    end if

                    rdist_surface = sqrt(rdist2)- 0.5*asigma(atype(i))

                    if(rdist_surface<rdist_surface_ref) then
                        rdist_surface_ref = rdist_surface
                    end if

                    if(rdist2<hicut2) then                            ! ignore the atom i if it is beyond the cutoff radius
                        sig2_rdist2 = asigma2_he(atype(i))/rdist2         ! if an atom is within the cut-off, this is a convinient place to calculate its Lennard-Jones interaction
                        rdist6 = sig2_rdist2*sig2_rdist2*sig2_rdist2                         ! with a helium atom placed in the center of the cubelet icount for later use in the helium volume calculation
                        rdist12 = rdist6*rdist6                           ! based in the second virial approach
                        lj_energy  = aeps_he(atype(i))*(rdist12-rdist6)
                        lattice_lj_he(icount) = lattice_lj_he(icount) + lj_energy
                    end if
                end do

                if(overlap.eqv..True.) then                          ! if in the previous cycle an overlap was detected, this whole cubelet is ignored
                    cycle
                end if

                lattice_space(j,k,l) = 1                          ! otherwise we add it to the list of geometrically accessible cubelets lattice_space(j,k,l) = 1
                ng_cubes = ng_cubes + 1
                g_cubes(ng_cubes) = icount

                lattice_rdist2(j,k,l) = rdist_surface_ref*rdist_surface_ref ! lattice_rdist2(j,k,l) stores the shortest squared distance between cubelet j, k, l and nearest atom (without overlap)

                if(rdist2_ref>0.25*asigma2_he(atype(amin))) then  ! next few lines detect if the cubelet is accessible to helium atom and update the list of
                    lattice_space_he(j,k,l) = 1                   ! helium accessible cubelets lattice_space_He(j,k,l)
                    !$omp atomic
                    nhe_cubes = nhe_cubes + 1
                    he_cubes(nhe_cubes) = icount
                end if

                if(rdist2_ref>asigma2_n(atype(amin))) then        ! next few lines detect if the cubelet is accessible to nitrogen atom and update the list of
                    lattice_space_n(j,k,l) = 1                    ! nitrogen accessible cubelets lattice_space_N(j,k,l)
                    nn_cubes = nn_cubes + 1
                    n_cubes(nn_cubes) = icount
                end if

            end do
        end do
    end do

    allocate(PA1(ng_cubes), PA2(ng_cubes), PA3(ng_cubes), PA4(ng_cubes))
    PA1=0.0; PA2=0; PA3=0; PA4=0

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Preliminary lattice complete                          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine lattice_calculations
!--------------------------------------------------------------------------------------
!	Here we use the information obtained in lattice_calculations
!	to generate percolated lattice accessible to point probe
!--------------------------------------------------------------------------------------

subroutine geom_lattice
    use lattice
    use percolation, only: percolation_calc

    implicit none

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of  point-accessible lattice               !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    call percolation_calc(lattice_space, g_cubes, ng_cubes, cl_summary)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of point-accessible lattice complete	 !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine geom_lattice


!--------------------------------------------------------------------------------------
!       Here we use the information obtained in lattice_calculations
!       to generate percolated lattice accessible to helium
!--------------------------------------------------------------------------------------

subroutine helium_lattice
    use lattice
    use percolation, only: percolation_calc

    implicit none

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of helium-accessible lattice               !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    call percolation_calc(lattice_space_he, he_cubes, nhe_cubes, cl_summary)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of helium-accessible lattice complete      !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine helium_lattice

!--------------------------------------------------------------------------------------
!       Here we use the information obtained in lattice_calculations
!       to generate percolated lattice accessible to nitrogen atom/probe
!--------------------------------------------------------------------------------------

subroutine nitrogen_lattice
    use parameters
    use lattice
    use percolation, only:  percolation_calc
    use fundcell, only:     fundcell_unslant
    use vector, only:       vectype
    use defaults, only:     radTodeg

    implicit none

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of nitrogen-accessible lattice             !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    call percolation_calc(lattice_space_n, n_cubes, nn_cubes, cl_summary)

    ! vusalization options: 1 -xyz, 2 - grd, 3 - both; 0 - none
    if(vis_option == 1) then 
    call nitrogen_lattice_vis(1)
    elseif(vis_option == 2) then
    call nitrogen_lattice_vis(2)
    elseif(vis_option == 3) then
    call nitrogen_lattice_vis(1)
    call nitrogen_lattice_vis(2)
    endif

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of nitrogen-accessible lattice complete    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine nitrogen_lattice

!--------------------------------------------------------------------------------------
!       In this subroutine we calculate helium-accessible volume,
!       using the second virial approach, geometric volume
!       accessible to a point probe, and free volume enclosed bu Connolly surface
!--------------------------------------------------------------------------------------

subroutine volumes

    use parameters
    use lattice
    use adsorbent
    use fundcell, only: fundamental_cell, fundcell_getvolume
    use results

    implicit none
    integer ::                            i, j, cube_number
    real*8  ::                            lj_energy, bf !, pore_v_he, pore_geom, volume  ! Moved these variables to the results module.

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting pore volume calculations                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    pore_v_he = 0.0

    do i=1, nhe_cubes

        cube_number = he_cubes(i)

        lj_energy = lattice_lj_he(cube_number)      ! total LJ energy of cubelet cube_number and its environment within cut-off distance
	lj_energy = 4.0 * lj_energy
        bf = exp(-lj_energy/temp)                   ! Boltzmann factor for cubelet cube_number
        pore_v_he = pore_v_he + bf                  ! second virial as a sum of Boltzmann factors over all cubelets
    end do

    pore_v_he = pore_v_he/dble(ntot)            ! averaging over the whole sample

    pore_geom = dble(ng_cubes)/dble(ntot)       ! geometric volume is simply proportional to all cubelets accessible to the point probe

    volume = fundcell_getvolume(fcell)          ! volume of the simulation cell

    pore_v_he = volume*pore_v_he
    pore_geom = volume*pore_geom


    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " helium volume in A^3:                    ",  pore_v_he
    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " helium volume in cm^3/g:                 "&
    ,  pore_v_he*0.602214129/smass
    write(*,*)
    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " geometric volume in A^3:                 ",  pore_geom
    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " geometric volume in cm^3/g:              "&
    ,  pore_geom*0.602214129/smass
    write(*,*)
    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " probe-occupiable volume in A^3:           ",  free_volume
    write(*,'(2a,f12.3)') " "//Trim(adjustl(property)), " probe-occupiable volume in cm^3/g:        "&
    , free_volume*0.602214129/smass
    write(*,'(2a,f12.5)') " "//Trim(adjustl(property)), " fraction of free volume:                 ",  ffv
    write(*,*)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Pore volume calculations complete                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    write(100,*) property
    write(100,'(a,f12.3)') "V_He_A^3 ",  pore_v_he
    write(100,'(a,f12.3)') "V_He_cm^3/g ", pore_v_he*0.602214129/smass
    write(100,'(a,f12.3)') "V_G_A^3 ", pore_geom
    write(100,'(a,f12.3)') "V_G cm^3/g: ",  pore_geom*0.602214129/smass
    write(100,'(a,f12.3)') "V_PO A^3: ",  free_volume
    write(100,'(a,f12.3)') "V_PO cm^3/g: ", free_volume*0.602214129/smass
    write(100,'(a,f12.5)') "FV_PO: ",  ffv
    
end subroutine volumes

!--------------------------------------------------------------------------------------
!       Subroutine for calculation of the accessible surface area
!--------------------------------------------------------------------------------------

subroutine surface_area

    use parameters
    use atoms
    use adsorbent
    use lattice

    use defaults, only:   rdbl, pi
    use fundcell, only:   fundamental_cell, fundcell_snglminimage, fundcell_getvolume, fundcell_maptocell, fundcell_slant
    use vector, only:     vectype
    Use random, only:     rranf
    use results

    implicit none
    real*8                                :: phi, costheta, theta
    !real*8                                :: stotal, volume  ! Moved these variables to the results module.
    real*8                                :: sjreal, sfrac, stotalreduced, rdist2
    integer                               :: nx, ny, nz, ncount, i, j, k, nx_temp, ny_temp, nz_temp
    type(VecType)                         :: atvec1, atvec2, sepvec, atvec_temp
    logical                               :: deny

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting surface area calculations                    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    if (nn_cubes == 0) then
        ! If there are no nitrogen-accessible points, the surface area is zero.
        stotal = 0.0
        volume = fundcell_getvolume(fcell)
        stotalreduced=stotal/(volume)*1.E4
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area in A^2:                 ', stotal
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area per volume in m^2/cm^3: ', stotalreduced
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area per mass in m^2/g:      ', &
        stotalreduced / ( smass/(0.6022141*volume))
        write(*,*) "No nitrogen accessible surface area "
        write(*,*)

        write(*,*) "!-------------------------------------------------------!"
        write(*,*) "! Surface area calculations complete                    !"
        write(*,*) "!-------------------------------------------------------!"
        write(*,*)

        write(100,*) property
        write(100,'(a,f12.2)') "S_AC_A^2 ", stotal
        write(100,'(a,f12.2)') "S_AC_m^2/cm^3 ", stotalreduced
        write(100,'(a,f12.2)') "S_AC_m^2/g ", stotalreduced / ( smass/(0.6022141*volume))

        return
    end if

    stotal = 0.0      ! initialize cumulative accessible surface area

    do i=1, natoms    ! number of atoms in the structure

        ncount = 0

        do j=1, nsample   ! number of sample points for each atom

            ! generate random vector of length 1
            ! first generate phi -pi pi

            phi=pi - rranf()*2.0*pi

            ! generate theta -pi:pi
            costheta = 1 - rranf() * 2.0
            theta = acos(costheta)
            atvec1%comp(1) = sin(theta) * cos(phi)
            atvec1%comp(2) = sin(theta) * sin(phi)
            atvec1%comp(3) = costheta

            ! make this vector of (sigma+probe_diameter)/2.0 length

            atvec1%comp = atvec1%comp * (coeff_surface * asigma_n(atype(i))) +  coords(:,i)

            ! translate the center of the coordinate to the particle i center and apply PBC

            ! apply PBCs to ensure that the selected point is within the simulation cell

            atvec_temp = atvec1
            
            if(fcell%orthoflag.eqv..True.) then
            atvec1 = fundcell_maptocell(fcell,atvec1)
            else
            atvec1 = fundcell_slant(fcell, atvec1)
            end if

            ! locate the cubelet in which the point sits

            nx = int(atvec1%comp(1)/cube_size) + 1
            ny = int(atvec1%comp(2)/cube_size) + 1
            nz = int(atvec1%comp(3)/cube_size) + 1
            
            if(nx>ncubesx) nx = nx - (nx/ncubesx)*ncubesx
            if(nx<1) nx = nx + (1-(nx/ncubesx))*ncubesx
            if(ny>ncubesy) ny = ny - (ny/ncubesy)*ncubesy
            if(ny<1) ny = ny +  (1-(ny/ncubesy))*ncubesy
            if(nz>ncubesz) nz = nz - (nz/ncubesz)*ncubesz
            if(nz<1) nz = nz +  (1-(nz/ncubesz))*ncubesz

            if(lattice_space_n(nx,ny,nz)<1) cycle ! reject the point if it is within an atom

            !----------------------
            ! now we perform the overlap test, i.e. we ensure that the selected point is not inside
            ! any other surf_coeff*(atom_diamter+probe_diameter)/2 distance in the system

            deny=.False.

            do k=1, natoms
                if(k==i) cycle

                !atvec2%comp = coords(:, k)

                Call fundcell_snglMinImage(fcell,atvec1,matvec(k),sepvec,rdist2)
                if(rdist2<coeff_surface2*asigma2_n(atype(k))) then
                    deny=.True.
                    exit
                end if
            end do
            !----------------------


            if(deny) cycle
            ncount=ncount+1
        end do

        ! fraction of the accessible surface area for sphere i
        sfrac=dble(ncount)/dble(nsample)


        ! surface area for sphere i in real units (A^2)
        sjreal=4.0*pi*coeff_surface2*asigma2_n(atype(i))*sfrac
        stotal=stotal+sjreal

    end do

    ! converting stotal on Surface per Volume

    volume = fundcell_getvolume(fcell)
    stotalreduced=stotal/(volume)*1.E4

    ! report results

    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area in A^2:                 ', stotal
    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area per volume in m^2/cm^3: ', stotalreduced
    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), ' surface area per mass in m^2/g:      ', &
    stotalreduced / ( smass/(0.6022141*volume))
    write(*,*)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Surface area calculations complete                    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    write(100,*) property
    write(100,'(a,f12.2)') "S_AC_A^2 ", stotal
    write(100,'(a,f12.2)') "S_AC_m^2/cm^3 ", stotalreduced
    write(100,'(a,f12.2)') "S_AC_m^2/g ", stotalreduced / ( smass/(0.6022141*volume))


    return

end subroutine surface_area


!--------------------------------------------------------------------------------------
!       Subroutine for calculation of the PSD
!--------------------------------------------------------------------------------------

subroutine pore_distribution

    use parameters
    use atoms
    use adsorbent
    use lattice
    use distributions
    use results


    use fundcell, only: fundcell_init, fundamental_cell, fundcell_getvolume, fundcell_snglminimage, fundcell_slant, fundcell_unslant
    Use vector, only: vectype
    Use random, Only: rranf


    implicit none
    type(vectype)                         :: atvec1, atvec2, sepvec
    integer                               :: i, j, k, l, nx, ny, nz, nx1, ny1, nz1, icount, m, bin, isite, ncycles, ivis
    real*8                                :: sigma_ref, sigma2_ref, sigma_ref2, rdist2
    real*8                                :: deldis1, deldis2, deldis
    real*8                                :: x, y, z
    integer,allocatable                   :: connolly_volume(:)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting pore size distribution calculations          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    psd_cumul=0.0
    psd=0.0
    icount = 0
    ivis = 0
    ncycles = 10000
    ffv = 0.0
    free_volume = 0.0

    allocate(connolly_volume(ncubesz*ncubesy*ncubesx))
    connolly_volume = 0

     volume = fundcell_getvolume(fcell)

    ! here we go through all cubelets accessible to nitrogen, store distances between the centers of cubelets
    ! and the surface of the nearest neighbour atom in PA1 array, and sort PA1 in an ascending order

    do l=1, ncubesz
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space_n(j,k,l)<1) cycle
                icount = icount + 1
                PA1(icount) = lattice_rdist2(j,k,l)
                PA2(icount) = j
                PA3(icount) = k
                PA4(icount) = l
            end do
        end do
    end do

    call sort(nn_cubes,PA1,PA2,PA3,PA4)

    if (nn_cubes == 0) then
        ! If there are no nitrogen-accessible cubes, the pore size is zero.
        ! psd() and psd_cumul() have already been initialised to zero.
        write(*,*)  " No nitrogen accessible porosity "
        write(*,*)  
        write(*,*) 

        write(*,*) "!-------------------------------------------------------!"
        write(*,*) "! Pore size distribution calculations complete          !"
        write(*,*) "!-------------------------------------------------------!"
        write(*,*)
        return
    end if

    sigma2_ref=0.0

    do i = 1, ncycles
    sigma2_ref=0.0

        isite = int(rranf()*dble(ng_cubes)) + 1 ! randomly select an available cubelet
        if(isite > ng_cubes) isite = ng_cubes

        nx = lattice_index(1, g_cubes(isite))
        ny = lattice_index(2, g_cubes(isite))
        nz = lattice_index(3, g_cubes(isite))

        atvec1%comp(1) = dble(nx-1)*cube_size+0.5*cube_size ! this is the center of the selected cubelet
        atvec1%comp(2) = dble(ny-1)*cube_size+0.5*cube_size
        atvec1%comp(3) = dble(nz-1)*cube_size+0.5*cube_size


        do j=nn_cubes, 1, -1  ! now we go through all  cubelets and see if point atvec1 is within the
            ! distance between the center of a cubelet and the surface of the nearest neighbour atom

            nx1 = PA2(j)
            ny1 = PA3(j)
            nz1 = PA4(j)

            atvec2%comp(1) = dble(nx1-1)*cube_size+0.5*cube_size
            atvec2%comp(2) = dble(ny1-1)*cube_size+0.5*cube_size
            atvec2%comp(3) = dble(nz1-1)*cube_size+0.5*cube_size

!!LS 09 01 2018            atvec2 = fundcell_slant(fcell, atvec2)
            call fundcell_snglminimage(fcell,atvec1,atvec2,sepvec,rdist2)

            if(rdist2>lattice_rdist2(nx1,ny1,nz1)) then    ! if not, cycle
                cycle
            else
                sigma2_ref = lattice_rdist2(nx1,ny1,nz1)       ! if yes, this will be the largest sphere

                sigma_ref =  sqrt(sigma2_ref)
                sigma_ref2 = 2.0*sigma_ref

                ffv = ffv + 1.0                                ! part of the free volume enclosed by Connolly surface
                ivis = ivis + 1
                connolly_volume(ivis) = isite
                exit                                           ! within which point  atvec1 can sit
            end if

        end do

        if(sigma2_ref==0.0) cycle

        sigma_ref = sqrt(sigma2_ref)                   ! sigma here is distance, not diameter
        bin=int(2.0*sigma_ref/binsize)+1               ! distribution bin (2 is needed to make sigma proper diameter)

        if (bin > ubound(psd_cumul, 1)) bin = ubound(psd_cumul, 1)
        do m=1, bin
            psd_cumul(m)=psd_cumul(m)+1                    ! update cumulative distribution
        end do

    end do

    ffv = ffv/real(ncycles)
    ffv = ffv*dble(ng_cubes)/dble(ntot)
    free_volume = ffv*volume


    open(13, file=Trim(adjustl(property))//'_psd_cumulative.txt', status='unknown') ! The file for cumulative Vp function
    open(14, file=Trim(adjustl(property))//'_psd.txt', status='unknown')            ! File containing pore size distribution

    write(13, *) '# Cumulative accessible volume distribution as a function of probe diameter'
    write(13, *) '# '
    write(13, *) '# d(probe)                Volume Fraction'

    psd_cumul=psd_cumul/psd_cumul(1)

    do i=1, nbins
        write(13,*) binsize*real(i-1)-binsize/2.0, psd_cumul(i)
    end do


    psd(1)=0.0
    psd(2)=0.0
    psd(nbins)=0.0

    do i=2, nbins-1               ! numerical differentiation where for calculation in point i, we use points i-1 and i+1
        deldis1=psd_cumul(i+1)
        deldis2=psd_cumul(i-1)
        deldis=deldis1-deldis2
        psd(i)=-1.0*(deldis)/(binsize*2.0) ! note that this way of taking derivatives is very crude
    end do

    write(14,*) '# Derivative distribution function -dV(r)/dr (or -dV(d)/dd) vs d'

    do i=2, nbins-1
        write(14,*) binsize*real(i-1)-binsize/2.0, psd(i)
    end do

    close(13)
    close(14)

    open(20, file = "probe_occupiable_volume.xyz", status = "unknown")

    write(20, *) ivis
    write(20, *)

    do i=1, ivis    
    isite = connolly_volume(i)

        nx = lattice_index(1, g_cubes(isite))
        ny = lattice_index(2, g_cubes(isite))
        nz = lattice_index(3, g_cubes(isite))

        x = dble(nx-1)*cube_size+0.5*cube_size ! this is the center of the selected cubelet
        y = dble(ny-1)*cube_size+0.5*cube_size
        z = dble(nz-1)*cube_size+0.5*cube_size

        x = x +  fcell%mx
        y = y +  fcell%my
        z = z +  fcell%mz

        atvec1%comp(1) = x
        atvec1%comp(2) = y
        atvec1%comp(3) = z

        atvec1 = fundcell_unslant(fcell, atvec1)

        write(20,*) "He ", atvec1%comp(1), atvec1%comp(2), atvec1%comp(3)
    end do
    close(20)

    
    write(*,*)  Trim(adjustl(property))," cumulative PSD and differential"
    write(*,*)  "PSD have been stored in files ", Trim(adjustl(property))//"_psd_cumulative.txt"
    write(*,*)  "and ", Trim(adjustl(property))//"_ psd.txt"
    write(*,*)
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Pore size distribution calculations complete          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine pore_distribution

!--------------------------------------------------------------------------------------
!       Subroutine for for limiting pore diameter and for maximum pore size
!--------------------------------------------------------------------------------------

subroutine limiting_diameter

    use parameters
    use atoms
    use lattice

    use percolation, only: percolation_calc_simple
    use results

    implicit none
    real*8  :: rhigh, rlow, rmiddle, rmax, rdiff, rdiff_old
    integer :: i, j, k, l, percolation_type
    integer :: icount

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting limiting pore diameter and maximum pore size !"
    write(*,*) "! analysis                                              !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    icount = 0

    do l=1, ncubesz
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space(j,k,l)<1) cycle
                icount = icount + 1
                PA1(icount) = lattice_rdist2(j,k,l)
                PA2(icount) = j
                PA3(icount) = k
                PA4(icount) = l
            end do
        end do
    end do

    call sort(ng_cubes,PA1,PA2,PA3,PA4)


    rmax  = maxval(PA1)
    rhigh = maxval(PA1)
    rlow  = minval(PA1)
    rdiff = 0.0

    do i=1, 100000 ! here we employ bisection method + percolation analysis to find the largest probe
        ! with respect to which the system remains lattice-percolated

        rdiff_old = rdiff
        rdiff = rhigh-rlow
        if(abs(rdiff)<0.25.or.(rdiff==rdiff_old)) exit ! convergence is achieved
        rmiddle = 0.5*(rhigh+rlow)

        lattice_temp = 0

        do j=ng_cubes, 1, -1
            if(PA1(j) < rmiddle) exit
            lattice_temp(PA2(j), PA3(j), PA4(j)) = 1
        end do

        call percolation_calc_simple(lattice_temp, cl_summary, spanning)

        if(spanning>0) then
            rlow  = PA1(j-1)
            percolation_type = spanning
            cycle
        else
            rhigh = PA1(j+1)
            cycle
        end if

    end do

    write(*,'(a,f12.2)') ' Pore limiting diameter in A: ', 2.0*sqrt(rhigh)
    write(*,'(a,f12.2)') ' Maximum pore diameter in A:  ', 2.0*sqrt(rmax)
    if(percolation_type==1) then
        write(*, *) 'The system is percolated in ', percolation_type, ' dimension (channels)'
    elseif(percolation_type==2) then
        write(*, *) 'The system is percolated in ', percolation_type, ' dimensions (slits)'
    elseif(percolation_type==3) then
        write(*, *) 'The system is percolated in ', percolation_type, ' dimensions (3D pores)'
    else
    end if

    ! Store the output in the results module.
    pore_lim = 2.0*sqrt(rhigh)
    pore_max = 2.0*sqrt(rmax)
    sys_perc = percolation_type

    write(*,*)
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Limiting pore diameter and maximum pore size          !"
    write(*,*) "! analysis complete                                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    write(100,'(a,f12.2)') "PLD_A ", 2.0*sqrt(rhigh)
    write(100,'(a,f12.2)') "LCD_A ", 2.0*sqrt(rmax)
    write(100, *) "D ", percolation_type
    return

end subroutine limiting_diameter

!--------------------------------------------------------------------------------------
!       Here we visualize the nitrogen accessible network
!       using two options: xyz format or grd format
!--------------------------------------------------------------------------------------

subroutine nitrogen_lattice_vis(option)
    use parameters
    use lattice
    use adsorbent
    use percolation, only:  percolation_calc
    use fundcell, only:     fundcell_unslant
    use vector, only:       vectype
    use defaults, only:     radTodeg

    implicit none
    integer                            :: option
    integer                            :: nx, ny, nz, nxlow, nxup, nylow, nyup, nzlow, nzup, i, j, k
    integer                            :: ncx, ncy, ncz
    real                               :: x, y, z, field
    type(vectype)                      :: atvec1

    
    if(option==1) then
    ! Creating a grid of lattice sites accessible to nitrogen probe in xyz fornat
    open(110, file='nitrogen_network.xyz', status='unknown')

    write(110,*) nn_cubes
    write(110,*)

    do i=1, nn_cubes
        nx = lattice_index(1, n_cubes(i))
        ny = lattice_index(2, n_cubes(i))
        nz = lattice_index(3, n_cubes(i))

        x = dble(nx-1)*cube_size+0.5*cube_size ! this is the center of the selected cubelet
        y = dble(ny-1)*cube_size+0.5*cube_size
        z = dble(nz-1)*cube_size+0.5*cube_size

        x = x +  fcell%mx
        y = y +  fcell%my
        z = z +  fcell%mz

        atvec1%comp(1) = x
        atvec1%comp(2) = y
        atvec1%comp(3) = z

        atvec1 = fundcell_unslant(fcell, atvec1)

        write(110,*) "N ", atvec1%comp(1), atvec1%comp(2), atvec1%comp(3)
    end do
    close(110)
    return
    endif

    if(option==2) then
    open(111, file='nitrogen_network.grd', status='unknown')

    write(111,*) "Poreblazer nitrogen accessible grid: field = r"
    write(111,*) "(1p,e12.5)"
    write(111,'(6f12.3)') fcell%ell(1), fcell%ell(2), fcell%ell(3), radTodeg*fcell%anglebc, radTodeg*fcell%angleac, &
    radTodeg*fcell%angleab

    ncx = ncubesx-1
    ncy = ncubesy-1
    ncz = ncubesz-1  

    write(111,*) ncx, ncy, ncz   

    write(111,*) "1 ",  " 0 ", ncx, " 0 ", ncy, " 0 ", ncz

   do k=1,  ncubesz
    do j=1, ncubesy
     do i=1, ncubesx

       if(lattice_rdist2(i,j,k) == 0.0) then
       field = 0.0
       write(111,'(e12.5)') field
       else
       field = sqrt(lattice_rdist2(i,j,k))
       write(111,'(e12.5)') field
       endif
    end do
    end do
    end do
    close(111)
    return
    end if

    return

end subroutine nitrogen_lattice_vis

!==============================================================

subroutine sort(n,arr,brr,crr,drr)
    integer n,m,nstack
    real*8 arr(n)
    integer brr(n),crr(n), drr(n)
    parameter (M=7,NSTACK=50)
    integer i,ir,j,jstack,k,l,istack(NSTACK)
    real*8 a,temp
    integer b, c, d, tempi
    jstack=0
    l=1
    ir=n
    1     if(ir-l.lt.M)then
    do 12 j=l+1,ir
    a=arr(j)
    b=brr(j)
    c=crr(j)
    d=drr(j)
    do 11 i=j-1,1,-1
    if(arr(i).le.a)goto 2
    arr(i+1)=arr(i)
    brr(i+1)=brr(i)
    crr(i+1)=crr(i)
    drr(i+1)=drr(i)
    11        continue
    i=0
    2         arr(i+1)=a
    brr(i+1)=b
    crr(i+1)=c
    drr(i+1)=d
    12      continue
    if(jstack.eq.0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
    else
    k=(l+ir)/2
    temp=arr(k)
    arr(k)=arr(l+1)
    arr(l+1)=temp
    tempi=brr(k)
    brr(k)=brr(l+1)
    brr(l+1)=tempi
    tempi=crr(k)
    crr(k)=crr(l+1)
    crr(l+1)=tempi
    tempi=drr(k)
    drr(k)=drr(l+1)
    drr(l+1)=tempi
    if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1)
    arr(l+1)=arr(ir)
    arr(ir)=temp
    tempi=brr(l+1)
    brr(l+1)=brr(ir)
    brr(ir)=tempi
    tempi=crr(l+1)
    crr(l+1)=crr(ir)
    crr(ir)=tempi
    tempi=drr(l+1)
    drr(l+1)=drr(ir)
    drr(ir)=tempi
    endif
    if(arr(l).gt.arr(ir))then
    temp=arr(l)
    arr(l)=arr(ir)
    arr(ir)=temp
    tempi=brr(l)
    brr(l)=brr(ir)
    brr(ir)=tempi
    tempi=crr(l)
    crr(l)=crr(ir)
    crr(ir)=tempi
    tempi=drr(l)
    drr(l)=drr(ir)
    drr(ir)=tempi
    endif
    if(arr(l+1).gt.arr(l))then
    temp=arr(l+1)
    arr(l+1)=arr(l)
    arr(l)=temp
    tempi=brr(l+1)
    brr(l+1)=brr(l)
    brr(l)=tempi
    tempi=crr(l+1)
    crr(l+1)=crr(l)
    crr(l)=tempi
    tempi=drr(l+1)
    drr(l+1)=drr(l)
    drr(l)=tempi
    endif
    i=l+1
    j=ir
    a=arr(l)
    b=brr(l)
    c=crr(l)
    d=drr(l)
    3       continue
    i=i+1
    if(arr(i).lt.a)goto 3
    4       continue
    j=j-1
    if(arr(j).gt.a)goto 4
    if(j.lt.i)goto 5
    temp=arr(i)
    arr(i)=arr(j)
    arr(j)=temp
    tempi=brr(i)
    brr(i)=brr(j)
    brr(j)=tempi
    tempi=crr(i)
    crr(i)=crr(j)
    crr(j)=tempi
    tempi=drr(i)
    drr(i)=drr(j)
    drr(j)=tempi
    goto 3
    5       arr(l)=arr(j)
    arr(j)=a
    brr(l)=brr(j)
    brr(j)=b
    crr(l)=crr(j)
    crr(j)=c
    drr(l)=drr(j)
    drr(j)=d
    jstack=jstack+2
    if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
    if(ir-i+1.ge.j-l)then
    istack(jstack)=ir
    istack(jstack-1)=i
    ir=j-1
    else
    istack(jstack)=j-1
    istack(jstack-1)=l
    l=i
    endif
    endif
    goto 1
end subroutine sort
!  (C) Copr. 1986-92 Numerical Recipes Software .2.
