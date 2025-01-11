!! @author Diana Powell
!! @version Oct-2019

program carma_day
  implicit none

  call test_day()

  write(*,*) "Done"
end program

subroutine test_day()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carmagas_mod
  use carmasolute_mod
  use carma_mod

  implicit none

  integer, parameter    :: NZ           = 89
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 18
  integer, parameter    :: NBIN         = 80
  integer, parameter    :: NGROUP       = 11
  integer, parameter    :: NSOLUTE      = 1
  integer, parameter    :: NGAS         = 10
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: ilongitude   = 64

  real(kind=f), parameter   :: dtime  = 100._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 200._f
  real(kind=f), parameter   :: zmin   = 0._f
  real(kind=f), parameter   :: rplanet = 6.991e9_f

 ! integer, parameter    :: irestart     = 1  ! =1 to restart
  integer, parameter    :: irestart     = 0  ! =1 to restart
 integer, parameter    :: idiag     = 1  ! =1 to output diagnostic
  ! integer, parameter    :: idiag     = 0  ! =1 to output diagnostic
  integer, parameter    :: iskip        = 1000 ! Output every iskip steps; no steps skipped if = 1
  integer, parameter    :: nstep        = 1000000

  integer, parameter        :: I_KCL       = 1
  integer, parameter        :: I_ZNS       = 2
  integer, parameter        :: I_NA2S      = 3
  integer, parameter        :: I_MNS       = 4
  integer, parameter        :: I_CR        = 5
  integer, parameter        :: I_MG2SIO4   = 6
  integer, parameter        :: I_FE        = 7
  integer, parameter        :: I_TIO2      = 8
  integer, parameter        :: I_AL2O3     = 9

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0

  real(kind=f), allocatable   :: xc(:)
  real(kind=f), allocatable   :: dx(:)
  real(kind=f), allocatable   :: yc(:)
  real(kind=f), allocatable   :: dy(:)
  real(kind=f), allocatable   :: zc(:)
  real(kind=f), allocatable   :: zl(:)
  real(kind=f), allocatable   :: dz_test(:)
  real(kind=f), allocatable   :: p(:)
  real(kind=f), allocatable   :: pl(:)
  real(kind=f), allocatable   :: t(:)
  real(kind=f), allocatable   :: rho_atm_cgs(:)

  real(kind=f), allocatable   :: mmr(:,:,:)
  real(kind=f), allocatable   :: numden(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: mmr_gas_old(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  real(kind=f), allocatable   :: satliq_old(:,:)
  real(kind=f), allocatable   :: satice_old(:,:)
  real(kind=f), allocatable   :: svpliq(:,:)
  real(kind=f), allocatable   :: wtpct(:,:)
  real(kind=f), allocatable   :: gflux(:,:)
  real(kind=f), allocatable   :: pflux(:,:,:)
  real(kind=f), allocatable   :: winds(:)
!  real(kind=f), allocatable   :: ekz(:)
  real(kind=f), allocatable   :: prodrate(:,:,:)
  real(kind=f), allocatable   :: prodrate_mass(:,:,:)
  real(kind=f), allocatable   :: prodrate_gas(:,:)
  real(kind=f), allocatable   :: totmass(:)

  real(kind=f), allocatable   :: ftopp(:,:)
  real(kind=f), allocatable   :: fbotp(:,:)
  real(kind=f), allocatable   :: pctop(:,:)
  real(kind=f), allocatable   :: pcbot(:,:)
  real(kind=f), allocatable   :: gctop(:)
  real(kind=f), allocatable   :: gcbot(:)
  real(kind=f), allocatable   :: ftopg(:)
  real(kind=f), allocatable   :: fbotg(:)

  real(kind=f), allocatable   :: gasprod(:,:)
  real(kind=f), allocatable   :: rhompe(:,:,:)
  real(kind=f), allocatable   :: growpe(:,:,:)
  real(kind=f), allocatable   :: evappe(:,:,:)
  real(kind=f), allocatable   :: growlg(:,:,:)
  real(kind=f), allocatable   :: evaplg(:,:,:)

  real(kind=f), allocatable   :: zsubsteps(:)

  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rlow(:)
  real(kind=f), allocatable   :: rup(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:,:)

  real(kind=f)          :: lat
  real(kind=f)          :: lon

  integer               :: i
  integer               :: j
  integer               :: istep, istep_old
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: ibinm
  integer               :: iz
  integer, parameter    :: lunerr = 48
  integer, parameter    :: lun = 42
  integer, parameter    :: lunp = 43
  integer, parameter    :: lunf = 44
  integer, parameter    :: lunfp = 45
  integer, parameter    :: lunres = 46
  integer, parameter    :: lundiagn = 47
  integer               :: nsubsteps
  integer               :: binmultiple
  integer               :: lastsub = 0

  integer, parameter    :: lun_atm = 49
  integer, parameter    :: lun_kzz = 50

  integer, parameter    :: lunrates = 55
  integer, parameter    :: lunratesp = 56

  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f
  real(kind=f)          :: t_orig

  real(kind=f)          :: t_in
  real(kind=f)          :: p_in
  real(kind=f)          :: pl_in
  real(kind=f)          :: z_in
  real(kind=f)          :: zl_in
  real(kind=f)          :: mu_in
  real(kind=f)          :: g_in
  real(kind=f)          :: kzz_in
  real(kind=f)          :: tio2_in
  real(kind=f)          :: t0_in

  real(kind=f)          :: met	! Metallicity = [Fe/H]

  real(kind=f)          :: startcd
  real(kind=f)          :: endcd
  real(kind=f)          :: inputrate
  real(kind=f)          :: vertpartflux
  real(kind=f)          :: vertgasflux

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, wtmol

  character(30)      	:: name
  character(30)      	:: gname
  character(len=3)	:: fileprefix = 'hj_'
  character(len=5)	:: temp = 'temp_'
  character(len=5)	:: flux = 'flux_'
  character(len=5)	:: diag = 'diag_'
  character(len=6)	:: rates = 'rates_'
  character(len=4)	:: filesuffix = '.txt'
  character(len=4)	:: filesuffix_restart = '.dat'
  character(len=100)	:: filename_restart = 'diamondback_test'
  character(len=100)	:: filename = 'diamondback_test'

  ! KCl     = 1
  ! ZnS     = 2
  ! Na2S    = 3
  ! MnS     = 4
  ! Cr      = 5
  ! Mg2SiO3 = 6
  ! Fe      = 7
  ! TiO2    = 8
  ! Al2O3   = 9

  ! Naming scheme:
  ! Heterogeneous Nucleation: (numbers of mantle material)nuc(number of core)
  ! Homogeneous Nucleation: (numbers of condensing material)hom

  real(kind=f)          :: tempr(NZ), pre(NZ), prel(NZP1), alt(NZ), altl(NZP1), wtmol_air(NZ), grav(NZ), ekz(NZP1), ekzl(NZP1)
  real(kind=f)          :: temp_equator(NZ, ilongitude), p_equator_center(NZ), p_equator_level(NZP1), velocity(ilongitude), longitudes(ilongitude)

  real(kind=f)          :: distance_btwn_elements, circumference, rotation_counter, slope, intercept
  real(kind=f)          :: current_distance, closeto_temp_profile, num_steps_btwn, current_step, RPLANET_DAT, velocity_avg

  write(*,*) ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12, file = 'diamondbacktest_centers.txt')

 do i = 1, NZ
 	read(12,*) alt(i), pre(i), tempr(i)
 end do

 close(12)

 open(13, file = 'diamondbacktest_centers.txt')

 do i=1, NZP1
   read(13,*) altl(i), prel(i), ekz(i)
 end do

 close(13)

!  open(14, file = 'w39b_kzz_big.txt')

!  do i=1, NZP1
!    read(14,*) ekzl(i)
!  end do

!  close(14)

! Taken from 2D Carma 1100 K
! do i = 1, NZP1
!  	if (pl(i) < 1.e5_f) then
! 		ekz(i) = 5.e7_f*(pl(i)*1.e-5_f)**(-0.45)
!  	else
!  		ekz(i) = 5.e7_f
!  	end if
! end do

  wtmol_air(:) = 2.3202_f
  grav(:) = 426._f !cm/s^2
  met = 1._f
  tio2_in = 5.913558445076336e-09_f

  t0_in = MAXVAL(tempr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Open the output text file
  open(unit=lun,file = fileprefix // filename(1:len_trim(filename)) // filesuffix, status="unknown")
  open(unit=lunf,file = fileprefix // flux // filename(1:len_trim(filename)) // filesuffix, status="unknown")
  open(unit=lunrates,file = fileprefix // rates // filename(1:len_trim(filename)) // filesuffix, status="unknown")

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), dz_test(NZ), p(NZ), pl(NZP1), &
           t(NZ), rho_atm_cgs(NZ))
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(numden(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(mmr_gas_old(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(satliq_old(NZ,NGAS))
  allocate(satice_old(NZ,NGAS))
  allocate(svpliq(NZ,NGAS))
  allocate(wtpct(NZ,NGAS))
  allocate(gflux(NZP1,NGAS))
  allocate(pflux(NZP1,NBIN,NELEM))
  allocate(winds(NZ))
!  allocate(ekz(NZP1))
  allocate(prodrate(NZ,NBIN,NELEM))
  allocate(prodrate_mass(NZ,NBIN,NELEM))
  allocate(prodrate_gas(NZ,NGAS))
  allocate(totmass(NZ))
  allocate(zsubsteps(NZ))
  allocate(r(NBIN))
  allocate(rlow(NBIN))
  allocate(rup(NBIN))
  allocate(dr(NBIN))
  allocate(rmass(NBIN,NGROUP))
  allocate(ftopp(NBIN,NELEM))
  allocate(fbotp(NBIN,NELEM))
  allocate(pctop(NBIN,NELEM))
  allocate(pcbot(NBIN,NELEM))
  allocate(gctop(NGAS))
  allocate(gcbot(NGAS))
  allocate(ftopg(NGAS))
  allocate(fbotg(NGAS))
  allocate(gasprod(NZ,NGAS))
  allocate(rhompe(NZ,NBIN,NELEM))
  allocate(growpe(NZ,NBIN,NELEM))
  allocate(evappe(NZ,NBIN,NELEM))
  allocate(growlg(NZ,NBIN,NGROUP))
  allocate(evaplg(NZ,NBIN,NGROUP))

  ! Define the particle-grid extent of the CARMA test
  write(*,*) "Create CARMA Object ..."

  if (idiag .eq. 1) then
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6, lundiag=lundiagn)
  else
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  end if
  if (rc < 0) stop "    *** FAILED in CARMA_Create ***"

	carma_ptr => carma

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the groups
  write(*,*) "  Add Group(s) ..."
  write(*,*) " "

  write(*,*) "Add TiO2 ..."
  call CARMAGROUP_Create(carma, 1, "Pure TiO2", 1e-8_f, 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Al2O3 ..."
  call CARMAGROUP_Create(carma, 2, "Al2O3 on TiO2", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Fe ..."
  call CARMAGROUP_Create(carma, 3, "Fe on TiO2", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Mg2SiO4 ..."
  call CARMAGROUP_Create(carma, 4, "Mg2SiO4 on TiO2", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Cr ..."
  call CARMAGROUP_Create(carma, 5, "Cr on TiO2", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add MnS ..."
  call CARMAGROUP_Create(carma, 6, "MnS on TiO2", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Na2S ..."
  call CARMAGROUP_Create(carma, 7, "Na2S on TiO2", 1e-8_f * 2._f**(2._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Fe ..."
  call CARMAGROUP_Create(carma, 8, "Pure Fe", 1e-8_f, 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Cr ..."
  call CARMAGROUP_Create(carma, 9, "Pure Cr", 1e-8_f, 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add KCl ..."
  call CARMAGROUP_Create(carma, 10, "Pure KCl", 1e-8_f, 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add ZnS ..."
  call CARMAGROUP_Create(carma, 11, "ZnS on KCl", 1e-8_f * 2._f**(1._f/3._f), 2._f, I_SPHERE, 1._f, &
                         .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the element
  write(*,*) "  Add Element(s) ..."
  write(*,*) " "

  write(*,*) "Add TiO2 ..."
  call CARMAELEMENT_Create(carma, 1, 1, "Pure TiO2", RHO_TIO2, I_VOLATILE, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Al2O3 Mantle ..."
  call CARMAELEMENT_Create(carma, 2, 2, "Al2O3 Mantle", RHO_AL2O3, I_VOLATILE, I_AL2O3, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (Al2O3 Cloud) ..."
  call CARMAELEMENT_Create(carma, 3, 2, "TiO2 Core (Al2O3)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Fe Mantle ..."
  call CARMAELEMENT_Create(carma, 4, 3, "Fe Mantle", RHO_FE, I_VOLATILE, I_FE, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (Fe Cloud) ..."
  call CARMAELEMENT_Create(carma, 5, 3, "TiO2 Core (Fe)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Mg2SiO4 Mantle ..."
  call CARMAELEMENT_Create(carma, 6, 4, "Mg2SiO4 Mantle", RHO_MG2SIO4, I_VOLATILE, I_MG2SIO4, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (Mg2SiO4 Cloud) ..."
  call CARMAELEMENT_Create(carma, 7, 4, "TiO2 Core (Mg2SiO4)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Cr Mantle ..."
  call CARMAELEMENT_Create(carma, 8, 5, "Cr Mantle", RHO_CR, I_VOLATILE, I_CR, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (Cr Cloud) ..."
  call CARMAELEMENT_Create(carma, 9, 5, "TiO2 Core (Cr)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add MnS Mantle ..."
  call CARMAELEMENT_Create(carma, 10, 6, "MnS Mantle", RHO_MNS, I_VOLATILE, I_MNS, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (MnS Cloud) ..."
  call CARMAELEMENT_Create(carma, 11, 6, "TiO2 Core (MnS)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Na2S Mantle ..."
  call CARMAELEMENT_Create(carma, 12, 7, "Na2S Mantle", RHO_NA2S, I_VOLATILE, I_NA2S, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Core (Na2S Cloud) ..."
  call CARMAELEMENT_Create(carma, 13, 7, "TiO2 Core (Na2S)", RHO_TIO2, I_COREMASS, I_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Fe ..."
  call CARMAELEMENT_Create(carma, 14, 8, "Pure Fe", RHO_FE, I_VOLATILE, I_FE, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Cr ..."
  call CARMAELEMENT_Create(carma, 15, 9, "Pure Cr", RHO_CR, I_VOLATILE, I_CR, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add KCl ..."
  call CARMAELEMENT_Create(carma, 16, 10, "Pure KCl", RHO_KCL, I_VOLATILE, I_KCL, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add ZnS Mantle ..."
  call CARMAELEMENT_Create(carma, 17, 11, "ZnS Mantle", RHO_ZNS, I_VOLATILE, I_ZNS, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add KCl Core (ZnS Cloud) ..."
  call CARMAELEMENT_Create(carma, 18, 11, "KCl Core (ZnS)", RHO_KCL, I_COREMASS, I_KCL, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the solutes
  write(*,*) "  Add Solute(s) ..."
  write(*,*) " "

  write(*,*) "Add H2SO4 Solute ..."
  call CARMASOLUTE_Create(carma, 1, "Sulfuric Acid", 2, WTMOL_H2SO4,RHO_H2SO4, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the gases
  write(*,*) "  Add Gase(s) ..."
  write(*,*) " "

  write(*,*) "Add H2O Vapour ..."
  call CARMAGAS_Create(carma, 1, "Water Vapour", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add TiO2 Vapour ..."
  call CARMAGAS_Create(carma, 2, "TiO2 Vapour", WTMOL_TIO2, I_VAPRTN_TIO2_HELLING2001, I_GCOMP_TIO2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Fe Vapour ..."
  call CARMAGAS_Create(carma, 3, "Fe Vapour", WTMOL_FE, I_VAPRTN_FE_VISSCHER2010, I_GCOMP_FE, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Mg2SiO4 Vapour ..."
  call CARMAGAS_Create(carma, 4, "Mg2SiO4 Vapour", WTMOL_MG2SIO4, I_VAPRTN_MG2SIO4_VISSCHER2010, I_GCOMP_MG2SIO4, rc, wtmol_dif=WTMOL_MG)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Cr Vapour ..."
  call CARMAGAS_Create(carma, 5, "Cr Vapour", WTMOL_CR, I_VAPRTN_CR_MORLEY2012, I_GCOMP_CR, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add MnS Vapour ..."
  call CARMAGAS_Create(carma, 6, "MnS Vapour", WTMOL_MNS, I_VAPRTN_MNS_MORLEY2012, I_GCOMP_MNS, rc, wtmol_dif=WTMOL_MN)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Na2S Vapour ..."
  call CARMAGAS_Create(carma, 7, "Na2S Vapour", WTMOL_NA2S, I_VAPRTN_NA2S_MORLEY2012, I_GCOMP_NA2S, rc, wtmol_dif=WTMOL_NA)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add ZnS Vapour ..."
  call CARMAGAS_Create(carma, 8, "ZnS Vapour", WTMOL_ZNS, I_VAPRTN_ZNS_MORLEY2012, I_GCOMP_ZNS, rc, wtmol_dif=WTMOL_ZN)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add KCl Vapour ..."
  call CARMAGAS_Create(carma, 9, "KCl Vapour", WTMOL_KCL, I_VAPRTN_KCL_MORLEY2012, I_GCOMP_KCL, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "Add Al2O3 Vapour ..."
  call CARMAGAS_Create(carma, 10, "Al2O3 Vapour", WTMOL_Al2O3, I_VAPRTN_AL2O3_WAKEFORD2017, I_GCOMP_AL2O3, rc, wtmol_dif=WTMOL_AL)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup the CARMA processes to exercise growth, nucleation, and coagulation.

  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 1, I_HOMGEN, 0._f, rc, igas=2)
  if (rc < 0) stop "    *** FAILED ***"

  !Al2O3 on TiO2 cos(43.6) = 0.724172 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 3, I_HETGEN, 0._f, rc, igas=10, ievp2elem=1, mucos = 0.724172_f )
  if (rc < 0) stop "    *** FAILED ***"

 !Fe on TiO2 cos(77.2) = 0.221548 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 5, I_HETGEN, 0._f, rc, igas=3, ievp2elem=1, mucos = 0.221548_f)
  if (rc < 0) stop "    *** FAILED ***"

  !Mg2SiO4 on TiO2 cos(0.1) = 0.995 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 7, I_HETGEN, 0._f, rc, igas=4, ievp2elem=1, mucos = 0.995_f)
  if (rc < 0) stop "    *** FAILED ***"

  !Cr on TiO2 cos(74.8) = 0.262189 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 9, I_HETGEN, 0._f, rc, igas=5, ievp2elem=1, mucos = 0.262189_f)
  if (rc < 0) stop "    *** FAILED ***"

  !MnS on TiO2 cos(77.6) = 0.214735 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 11, I_HETGEN, 0._f, rc, igas=6, ievp2elem=1, mucos = 0.214735_f)
  if (rc < 0) stop "    *** FAILED ***"

  !Na2S on TiO2 cos(61) = 0.48481 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 1, 13, I_HETGEN, 0._f, rc, igas=7, ievp2elem=1, mucos = 0.48481_f)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 14, 14, I_HOMGEN, 0._f, rc, igas=3)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 15, 15, I_HOMGEN, 0._f, rc, igas=5)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 16, 16, I_HOMGEN, 0._f, rc, igas=9)
  if (rc < 0) stop "    *** FAILED ***"

  !ZnS on KCl cos(81.7) = 0.144356 from Gao+2020
  write(*,*) "  Add Nucleation ..."
  call CARMA_AddNucleation(carma, 16, 18, I_HETGEN, 0._f, rc, igas=8, ievp2elem=16, mucos = 0.144356_f)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 1, 2, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 2, 10, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 4, 3, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 6, 4, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 8, 5, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 10, 6, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 12, 7, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 14, 3, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 15, 5, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 16, 9, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) "  Add Growth ..."
  call CARMA_AddGrowth(carma, 17, 8, rc)
  if (rc < 0) stop "    *** FAILED ***"

 ! write(*,*) "  Add Coagulation ..."
 ! call CARMA_AddCoagulation(carma, 2, 2, 2, I_COLLEC_FUCHS, rc)
 ! if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup the CARMA processes to exercise
  write(*,*) "Initialize CARMA ..."

  call CARMA_Initialize(carma, rc, do_cnst_rlh =.FALSE., do_coag=.FALSE., do_fixedinit=.TRUE., do_grow=.TRUE., &
                        do_explised=.FALSE., do_substep=.TRUE., do_print_init=.TRUE., &
                        do_vdiff=.TRUE., do_vtran=.TRUE., maxsubsteps=10, maxretries=20, &
                        itbnd_pc=I_FLUX_SPEC, ibbnd_pc=I_FIXED_CONC, itbnd_gc=I_FLUX_SPEC, ibbnd_gc=I_FIXED_CONC)
  if (rc < 0) stop "    *** FAILED CARMA_Initialize ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.

  ! Setup model atmosphere

  lat = 45.0_f
  lon = -105.0_f

  ! Horizonal centers
  dx(:) = deltax
  xc(:) = dx(:) / 2._f
  dy(:) = deltay
  yc(:) = dy(:) / 2._f

  ! Vertical center
  do i = 1, NZ
    zc(i) = alt(i)
    t(i) = tempr(i)
    p(i) = pre(i)
  end do

  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
  do i = 1, NZ
    write(*,'(i6,3f12.3)') i, zc(i), p(i), t(i)
  end do

  ! Vertical edge
  do i = 1, NZP1
    zl(i) = altl(i)
    pl(i) = prel(i)
  end do

  write(*,*) ""
  write(*,'(a6, 2a12)') "level", "zl", "pl"
  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
  do i = 1, NZP1
    write(*,'(i6,2f12.3)') i, zl(i), pl(i)
  end do

  ! ekz(:) = ekzl(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup up a mass mixing ratio of H20 and H2SO4 vapor

  mmr_gas(:,:) = 1e-50_f

  mmr_gas_old = mmr_gas

  ! Setup an intial mmr for condensation nuclei
  mmr(:,:,:) = 0._f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Write output for the test

  if (idiag .eq. 1) then
    open(unit=lundiagn,file = fileprefix // diag // filename(1:len_trim(filename)) // filesuffix, status="unknown")
    write(lundiagn,*) 'PART 0: HEADER'
    write(lundiagn,*) 'SIMULATION TYPE:'
    write(lundiagn,*) 'VENUS'
    write(lundiagn,*) 'DIMENSIONS:'
    write(lundiagn,'(6A9)') 'NZ', 'NGROUP', 'NELEM', 'NBIN', 'NGAS', 'NSTEP'
    write(lundiagn,'(6I9)') NZ, NGROUP, NELEM, NBIN, NGAS, nstep + 1
    write(lundiagn,*) 'GROUP INFORMATION:'
    write(lundiagn,'(A13,2A10,A17,A14,A19,A15)') 'IGROUP', 'IBIN', 'R', 'MASS', 'dR', 'R_LOWBOUND', 'R_UPBOUND'
    write(lundiagn,'(A10,A12,A15,A11,A19,2A15)') '','','(microns)', '(g)','(microns)','(microns)','(microns)'
  end if

  write(lun,'(7i10)') NZ, NGROUP, NELEM, NBIN, NGAS, nstep + 1, iskip

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rlow=rlow, rup=rup, dr=dr, rmass=rmass(:,igroup))
    if (rc < 0) stop "    *** FAILED CARMAGROUP_Get ***"

    do ibin = 1, NBIN
      write(lun,'(2i4,5e15.5)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin,igroup), dr(ibin) * 1e4_f, rlow(ibin) * 1e4_f, rup(ibin) * 1e4_f
      if (idiag .eq. 1) then
        write(lundiagn,'(i10,i12,5e15.5)') igroup, ibin, r(ibin) * 1e4_f, &
	rmass(ibin,igroup), dr(ibin) * 1e4_f, rlow(ibin) * 1e4_f, rup(ibin) * 1e4_f
      end if
    end do
  end do

  if (idiag .eq. 1) then
    write(lundiagn,*) 'ELEMENT INFORMATION:'
    write(lundiagn,'(2A10,A18)') 'IELEM', 'IGROUP', 'NAME'

    do ielem = 1, NELEM
      call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, name=name)
      if (rc < 0) stop "    *** FAILED CARMA_ELEMENT_Get ***"
      write(lundiagn,'(i8,i10,A35)') ielem, igroup, name
    end do

    write(lundiagn,*) 'GAS INFORMATION:'
    write(lundiagn,'(A6,A15,A35)') 'IGAS', 'NAME', 'WTMOL (g/mol)'

    do igas = 1, NGAS
      call CARMAGAS_Get(carma, igas, rc, name=gname, wtmol=wtmol)
      if (rc < 0) stop "    *** FAILED CARMAGAS_Get ***"
      write(lundiagn,'(i4,A40,e10.3)') igas, gname, wtmol
    end do

    write(lundiagn,*) 'ATMOSPHERE INFORMATION:'
    write(lundiagn,'(A5,4A15)') 'Z', 'ALTITUDE', 'dALT', 'PRESSURE', 'TEMPERATURE'
    write(lundiagn,'(A5,4A15)') '', '(km)', '(km)', '(mbars)', '(K)'
  end if

  do i = 1, NZ
    write(lun,'(i3,5e15.5)') i, zc(i), zl(i+1)-zl(i), p(i) * 10._f, t(i), ekz(i)
    if (idiag .eq. 1) then
      write(lundiagn,'(i5,4e15.5)') i, zc(i) / 1000._f, (zl(i+1)-zl(i)) / 1000._f, p(i) / 100._f, t(i)
    end if
  end do

  if (idiag .eq. 1) then
    write(lundiagn,*) ' '
  end if

  write(*,*) ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define production rates of particles and gases

  prodrate(:,:,:) = 0._f
  prodrate_gas(:,:) = 0._f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define constant wind speeds (temporally varying winds defined within step loop below)

  dz_test(:) = abs(zl(2:NZP1) - zl(1:NZ)) * 100._f
  rho_atm_cgs = p(:) / (RGAS/wtmol_air(:) * t(:))

  !met = 0._f

  !do i = 1, NZ

    !winds(i) = 0._f  				   ! No wind

    !winds(i) = 8.0e-5_f / rho_atm_cgs(i)  	   ! Upward Hadley circulation (no cut off)

    !if (i .lt. 105) then                         ! Upward Hadley circulation (gaussian cut off at 70 km)
    !  winds(i) = 8.0e-5_f / rho_atm_cgs(i)
    !else
    !  winds(i) = (8.0e-5_f / rho_atm_cgs(i)) * exp( - (((zc(i) - zc(105)) / 5000._f) ** 2) / 2 )
    !end if

    !winds(i) = -8.0e-5_f / rho_atm_cgs(i) 	   ! Downward Hadley circulation (no cut off)

    !if (i .lt. NZ/2) then                          ! Upward gust (linear cut off at 75 km)
    !  winds(i) = 8.0e-3_f / rho_atm_cgs(i)
    !else
    !  winds(i) = (8.0e-3_f - (i - NZ/2)*3.2e-4_f) / rho_atm_cgs(i)
    !  if (winds(i) .lt. 0._f) winds(i) = 0._f
    !end if

    !if (i .lt. NZ/2) then                          ! Upward gust (gaussian cut off at 75 km)
    !  winds(i) = 8.0e-3_f / rho_atm_cgs(i)
    !else
    !  winds(i) = (8.0e-3_f / rho_atm_cgs(i)) *  exp( - (((zc - zc(NZ/2)) / 4000._f) ** 2) / 2 )
    !end if

  !end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Boundary Conditions

  ftopp(:,:) = 0._f
  fbotp(:,:) = 0._f
  pctop(:,:) = 0._f
  pcbot(:,:) = 0._f
  gctop(:) = 1e-50_f
  gcbot(:) = 1e-50_f
  ftopg(:) = 0._f
  fbotg(:) = 0._f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(lun,*) 0
  write(lunf,*) 0
  do j = 1, NBIN
   do i = 1, NZ
    write(lun,'(i3,i4,29e11.3,f8.0)') &
     j, i, &
     real(mmr(i,1,j) * rho_atm_cgs(i) / rmass(j,1)), &
     real(mmr(i,2,j) * rho_atm_cgs(i) / rmass(j,2)), &
     real(mmr(i,4,j) * rho_atm_cgs(i) / rmass(j,3)), &
     real(mmr(i,6,j) * rho_atm_cgs(i) / rmass(j,4)), &
     real(mmr(i,8,j) * rho_atm_cgs(i) / rmass(j,5)), &
     real(mmr(i,10,j) * rho_atm_cgs(i) / rmass(j,6)), &
     real(mmr(i,12,j) * rho_atm_cgs(i) / rmass(j,7)), &
     real(mmr(i,14,j) * rho_atm_cgs(i) / rmass(j,8)), &
     real(mmr(i,15,j) * rho_atm_cgs(i) / rmass(j,9)), &
     real(mmr(i,16,j) * rho_atm_cgs(i) / rmass(j,10)), &
     real(mmr(i,17,j) * rho_atm_cgs(i) / rmass(j,11)), &
     real(mmr_gas(i,2) * 1.0e6_f / (WTMOL_TIO2 / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,3) * 1.0e6_f / (WTMOL_FE / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,4) * 1.0e6_f / (WTMOL_MG / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,5) * 1.0e6_f / (WTMOL_CR / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,6) * 1.0e6_f / (WTMOL_MN / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,7) * 1.0e6_f / (WTMOL_NA / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,8) * 1.0e6_f / (WTMOL_ZN / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,9) * 1.0e6_f / (WTMOL_KCL / wtmol_air(i))), 0._f, &
     real(mmr_gas(i,10) * 1.0e6_f / (WTMOL_AL / wtmol_air(i))), 0._f, 0
   end do
  end do

  binmultiple = int(NBIN / 10._f)

  if (irestart .eq. 1) then
    !open(unit=lunres,file="venus_atmosphere_300z_45bins_noinit_lf_h2ofixed_10nmCNs_t0030et_2.dat",form='unformatted',status="unknown")
    open(unit=lunres,file=fileprefix // filename_restart(1:len_trim(filename_restart)) // filesuffix_restart,form='unformatted',status="unknown")
    read(lunres) istep_old, xc, dx, yc, dy, &
         zc, zl, p, pl, t, rho_atm_cgs, &
         mmr, mmr_gas, mmr_gas_old, satliq, &
         satice, satliq_old, satice_old, svpliq, wtpct, &
         zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
	 prodrate, prodrate_mass, prodrate_gas, totmass, &
	 winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
	 gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
         evappe, growlg, evaplg
     write(*,*)'read restart file'
     rewind(lunres)
!     istep = istep + 1
    close(unit=lunres)
  endif

  call CARMASTATE_CreateFromReference(cstate, carma_ptr, time, dtime, NZ, &
                         	      I_CART, I_CART, lat, lon, &
                         	      xc(:), dx(:), &
                         	      yc(:), dy(:), &
                         	      zc(:), zl(:), p(:), &
                         	      pl(:), t(:), wtmol_air(:), grav(:), rplanet, rc, winds=winds(:), ekz=ekz(:), met=met, t0 = t0_in)
  if (rc < 0) stop "    *** FAILED CARMASTATE_CreateFromReference ***"

  ! Iterate the model over a few time steps.
  do istep = 1, nstep

    !open(unit=lunres,file="venus_atmosphere_300z_45bins_noinit_lf_h2ofixed_10nmCNs_t0030et_2.dat",form='unformatted',status="unknown")
    open(unit=lunres,file=fileprefix // filename(1:len_trim(filename)) // filesuffix_restart,form='unformatted',status="unknown")
    open(unit=lunp,file = fileprefix // temp // filename(1:len_trim(filename)) // filesuffix, status="unknown")
    open(unit=lunfp,file = fileprefix // temp // flux // filename(1:len_trim(filename)) // filesuffix, status="unknown")
    open(unit=lunratesp,file = fileprefix // temp // rates // filename(1:len_trim(filename)) // filesuffix, status="unknown")

    ! Calculate the model time.
    time = (istep - 1) * dtime

    write(*,*) 'istep, time', istep, time, filename(1:len_trim(filename))

    if (idiag .eq. 1) then
      do iz = 1, NZ
        totmass(iz) = sum((mmr(iz,1,:)+mmr(iz,2,:)) * abs((pl(iz+1) - pl(iz))) / 88.7_f) / deltaz / 100._f
      end do
      write(lundiagn,'(A6,I10,A12,f15.2,A8)') 'STEP:', istep, 'TIME:', istep*dtime, 'SECONDS'
      write(lundiagn,*) ' '
      write(lundiagn,*) 'PART 1: TOTAL MASS AT START OF TIME STEP'
      write(lundiagn,*) '****************************************'
      write(lundiagn,*) 'Z: ALTITUDE LEVEL INDEX'
      write(lundiagn,*) 'PARTMASS: TOTAL MASS DENSITY OF GAS AT Z [g/cm3]'
      write(lundiagn,*) 'GASMASS: TOTAL MASS DENSITY OF PARTICLES AT Z [g/cm3]'
      write(lundiagn,*) '****************************************'
      write(lundiagn,'(A6,2A15)') 'Z', 'PARTMASS', 'GASMASS'
      do iz = 1, NZ
	write(lundiagn,'(i6,2e15.3)') iz, totmass(iz), (mmr_gas(iz,1)+mmr_gas(iz,2)) * abs((pl(iz+1) - pl(iz))) / 88.7_f / deltaz / 100._f
      end do
      write(lundiagn,*) ' '
      write(lundiagn,'(A6,2e15.3)') 'TOT:',sum(totmass), &
                                    sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f / deltaz / 100._f )
      write(lundiagn,*) ' '
      write(lundiagn,*) 'COLUMN TOTAL MASS = TOT * NZ * deltaz (in meters) * (100 cm / m)'
      startcd = sum(totmass) * NZ * deltaz * 100._f + sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f ) * NZ
      write(lundiagn,'(A28,e28.15)') 'COLUMN TOTAL MASS [g/cm2]: ', startcd
      write(lundiagn,*) ' '
    end if

 !   mmr(:,1,1) = mmr(:,1,1) + prodrate(:,1) * dtime / rho_atm_cgs
 !   mmr(:,1,10) = mmr(:,1,10) + prodrate(:,10) * dtime / rho_atm_cgs
 !   mmr_gas(:,3) = mmr_gas(:,3) + prodrate_gas(:,3) * dtime / rho_atm_cgs(:)
 !   mmr_gas(:,2) = mmr_gas(:,2) + prodrate_gas(:,2) * dtime / rho_atm_cgs(:)


    ! To do: change gas input rate; add gaussian distribution to size of CNs being added; change nucleation rate with
    ! substepping; look closer at eddy diffusion machinery, share with Bardeen; read Bardeen's emails more!

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                           I_CART, I_CART, lat, lon, &
                           xc(:), dx(:), &
                           yc(:), dy(:), &
                           zc(:), zl(:), p(:), &
                           pl(:), t(:), wtmol_air(:), grav(:), rplanet, rc, told=t(:), winds=winds(:), ekz=ekz(:), &
			   ftopp=ftopp,fbotp=fbotp,pctop=pctop,pcbot=pcbot, &
			   gctop=gctop,gcbot=gcbot,ftopg=ftopg,fbotg=fbotg,met=met)
    if (rc < 0) stop "    *** FAILED CARMASTATE_Create ***"


    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc < 0) stop "    *** FAILED CARMASTATE_SetBin ***"
      end do
    end do

    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc, mmr_old=mmr_gas_old(:,igas), &
                             satice_old=satice(:,igas), satliq_old=satliq(:,igas))
      if (rc < 0) stop "    *** FAILED CARMASTATE_SetGas ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc < 0) stop "    *** FAILED CARMASTATE_Step ***"

    ! Get the retry stats and the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries, zsubsteps=zsubsteps)
    if (rc < 0) stop "    ***  FAILED CARMASTATE_Get ***"

    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc, numberDensity=numden(:,ielem,ibin), pflux=pflux(:,ibin,ielem))
        if (rc < 0) stop "    *** FAILED CARMASTATE_GetBin ***"
      end do
    end do

    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, mmr_gas(:,igas), rc, satliq=satliq(:,igas), &
                             satice=satice(:,igas), eqliq = svpliq(:,igas), wtpct=wtpct(:,igas), gflux=gflux(:,igas))
      if (rc < 0) stop "    *** FAILED CARMASTATE_GetGas ***"
    end do

    call CARMASTATE_GetDiag(cstate,rc,rhompe_tot=rhompe,growpe_tot=growpe,evappe_tot=evappe, &
                            growlg_tot=growlg,evaplg_tot=evaplg,gasprod_tot=gasprod)
    if (rc < 0) stop "    *** FAILED CARMASTATE_GetDiag ***"

    mmr_gas_old = mmr_gas

    ! Get the updated temperature.
!    call CARMASTATE_GetState(cstate, rc, t=t(:))
!    if (rc < 0) stop "    *** FAILED ***"

    ! Write output

 ! write(*,*) rmass


     if (istep .eq. 1) then
       !mmr_gas(1,2) = min(1.69e-7_f * (WTMOL_TIO2 / wtmol_air(1)) * 10._f**met,svpliq(1,2) * (WTMOL_TIO2 / wtmol_air(1)))
       mmr_gas(1,2) = min(tio2_in * (WTMOL_TIO2 / wtmol_air(1)),svpliq(1,2) * (WTMOL_TIO2 / wtmol_air(1)))
       mmr_gas(1,3) = min(5.78e-5_f * (WTMOL_FE / wtmol_air(1)) * 10._f**met,svpliq(1,3) * (WTMOL_FE / wtmol_air(1)))
       mmr_gas(1,4) = min(5.936e-5_f * (WTMOL_MG / wtmol_air(1)) * 10._f**met,svpliq(1,4) * (WTMOL_MG / wtmol_air(1)))
       mmr_gas(1,5) = min(8.87e-7_f * (WTMOL_CR / wtmol_air(1)) * 10._f**met,svpliq(1,5) * (WTMOL_CR / wtmol_air(1)))
       mmr_gas(1,6) = min(5.41e-7 * (WTMOL_MN / wtmol_air(1)) * 10._f**met,svpliq(1,6) * (WTMOL_MN / wtmol_air(1)))
       mmr_gas(1,7) = min(3.34e-6_f * (WTMOL_NA / wtmol_air(1)) * 10._f**met,svpliq(1,7) * (WTMOL_NA / wtmol_air(1)))
       mmr_gas(1,8) = min(7.65e-8_f * (WTMOL_ZN / wtmol_air(1)) * 10._f**met,svpliq(1,8) * (WTMOL_ZN / wtmol_air(1)))
       mmr_gas(1,9) = min(2.2e-7_f * (WTMOL_KCL / wtmol_air(1)) * 10._f**met,svpliq(1,9) * (WTMOL_KCL / wtmol_air(1)))
       mmr_gas(1,10) = min(4.937e-6_f * (WTMOL_AL / wtmol_air(1)) * 10._f**met,svpliq(1,10) * (WTMOL_AL / wtmol_air(1)))
       !gcbot(2) = min(1.69e-7_f * (WTMOL_TIO2 / wtmol_air(1)) * 10._f**met,svpliq(1,2) * (WTMOL_TIO2 / wtmol_air(1)))
       gcbot(2) = min(tio2_in * (WTMOL_TIO2 / wtmol_air(1)),svpliq(1,2) * (WTMOL_TIO2 / wtmol_air(1)))
       gcbot(3) = min(5.78e-5_f * (WTMOL_FE / wtmol_air(1)) * 10._f**met,svpliq(1,3) * (WTMOL_FE / wtmol_air(1)))
       gcbot(4) = min(5.936e-5_f * (WTMOL_MG / wtmol_air(1)) * 10._f**met,svpliq(1,4) * (WTMOL_MG / wtmol_air(1)))
       gcbot(5) = min(8.87e-7_f * (WTMOL_CR / wtmol_air(1)) * 10._f**met,svpliq(1,5) * (WTMOL_CR / wtmol_air(1)))
       gcbot(6) = min(5.41e-7 * (WTMOL_MN / wtmol_air(1)) * 10._f**met,svpliq(1,6) * (WTMOL_MN / wtmol_air(1)))
       gcbot(7) = min(3.34e-6_f * (WTMOL_NA / wtmol_air(1)) * 10._f**met,svpliq(1,7) * (WTMOL_NA / wtmol_air(1)))
       gcbot(8) = min(7.65e-8_f * (WTMOL_ZN / wtmol_air(1)) * 10._f**met,svpliq(1,8) * (WTMOL_ZN / wtmol_air(1)))
       gcbot(9) = min(2.2e-7_f * (WTMOL_KCL / wtmol_air(1)) * 10._f**met,svpliq(1,9) * (WTMOL_KCL / wtmol_air(1)))
       gcbot(10) = min(4.937e-6_f * (WTMOL_AL / wtmol_air(1)) * 10._f**met,svpliq(1,10) * (WTMOL_AL / wtmol_air(1)))
     endif


    if (MOD (istep, iskip) .eq. 0) then

      write(*,*) 'Recorded'
      write(lun,*) (istep)*dtime
      write(lunp,*) (istep)*dtime
      write(lunf,*) (istep)*dtime
      write(lunfp,*) (istep)*dtime
      write(lunrates,*) (istep)*dtime
      write(lunratesp,*) (istep)*dtime

      do j = 1, NBIN
        do i = 1, NZ
          write(lun,'(i3,i4,29e11.3,f8.0)') &
          j, i, &
          real(numden(i,1,j)), &
          real(numden(i,2,j)), &
          real(numden(i,4,j)), &
          real(numden(i,6,j)), &
          real(numden(i,8,j)), &
          real(numden(i,10,j)), &
          real(numden(i,12,j)), &
          real(numden(i,14,j)), &
          real(numden(i,15,j)), &
          real(numden(i,16,j)), &
          real(numden(i,17,j)), &
          real(mmr_gas(i,2) * 1.0e6_f / (WTMOL_TIO2 / wtmol_air(i))), &
          real(svpliq(i,2) * 1.0e6_f), &
          real(mmr_gas(i,3) * 1.0e6_f / (WTMOL_FE / wtmol_air(i))), &
          real(svpliq(i,3) * 1.0e6_f), &
          real(mmr_gas(i,4) * 1.0e6_f / (WTMOL_MG / wtmol_air(i))), &
          real(svpliq(i,4) * 1.0e6_f), &
          real(mmr_gas(i,5) * 1.0e6_f / (WTMOL_CR / wtmol_air(i))), &
          real(svpliq(i,5) * 1.0e6_f), &
          real(mmr_gas(i,6) * 1.0e6_f / (WTMOL_MN / wtmol_air(i))), &
          real(svpliq(i,6) * 1.0e6_f), &
          real(mmr_gas(i,7) * 1.0e6_f / (WTMOL_NA / wtmol_air(i))), &
          real(svpliq(i,7) * 1.0e6_f), &
          real(mmr_gas(i,8) * 1.0e6_f / (WTMOL_ZN / wtmol_air(i))), &
          real(svpliq(i,8) * 1.0e6_f), &
          real(mmr_gas(i,9) * 1.0e6_f / (WTMOL_KCL / wtmol_air(i))), &
          real(svpliq(i,9) * 1.0e6_f), &
          real(mmr_gas(i,10) * 1.0e6_f / (WTMOL_AL / wtmol_air(i))), &
          real(svpliq(i,10) * 1.0e6_f), &
          zsubsteps(i)
          write(lunp,'(i3,i4,29e11.3,f8.0)') &
          j, i, &
          real(numden(i,1,j)), &
          real(numden(i,2,j)), &
          real(numden(i,4,j)), &
          real(numden(i,6,j)), &
          real(numden(i,8,j)), &
          real(numden(i,10,j)), &
          real(numden(i,12,j)), &
          real(numden(i,14,j)), &
          real(numden(i,15,j)), &
          real(numden(i,16,j)), &
          real(numden(i,17,j)), &
          real(mmr_gas(i,2) * 1.0e6_f / (WTMOL_TIO2 / wtmol_air(i))), &
          real(svpliq(i,2) * 1.0e6_f), &
          real(mmr_gas(i,3) * 1.0e6_f / (WTMOL_FE / wtmol_air(i))), &
          real(svpliq(i,3) * 1.0e6_f), &
          real(mmr_gas(i,4) * 1.0e6_f / (WTMOL_MG / wtmol_air(i))), &
          real(svpliq(i,4) * 1.0e6_f), &
          real(mmr_gas(i,5) * 1.0e6_f / (WTMOL_CR / wtmol_air(i))), &
          real(svpliq(i,5) * 1.0e6_f), &
          real(mmr_gas(i,6) * 1.0e6_f / (WTMOL_MN / wtmol_air(i))), &
          real(svpliq(i,6) * 1.0e6_f), &
          real(mmr_gas(i,7) * 1.0e6_f / (WTMOL_NA / wtmol_air(i))), &
          real(svpliq(i,7) * 1.0e6_f), &
          real(mmr_gas(i,8) * 1.0e6_f / (WTMOL_ZN / wtmol_air(i))), &
          real(svpliq(i,8) * 1.0e6_f), &
          real(mmr_gas(i,9) * 1.0e6_f / (WTMOL_KCL / wtmol_air(i))), &
          real(svpliq(i,9) * 1.0e6_f), &
          real(mmr_gas(i,10) * 1.0e6_f / (WTMOL_AL / wtmol_air(i))), &
          real(svpliq(i,10) * 1.0e6_f), &
          zsubsteps(i)
!          write(lunrates,'(i3,i4,6e11.3)') &
!          j, i, real(rhompe(i,j,1) * rmass(j,1)), &
!          real(growpe(i,j,1) * rmass(j,1)), &
!          real(evappe(i,j,1) * rmass(j,1)), &
!          real(growlg(i,j,1) * rmass(j,1)), &
!          real(evaplg(i,j,1) * rmass(j,1)), &
!          real(gasprod(i,2))
!          write(lunratesp,'(i3,i4,6e11.3)') &
!          j, i, real(rhompe(i,j,1) * rmass(j,1)), &
!          real(growpe(i,j,1) * rmass(j,1)), &
!          real(evappe(i,j,1) * rmass(j,1)), &
!          real(growlg(i,j,1) * rmass(j,1)), &
!          real(evaplg(i,j,1) * rmass(j,1)), &
!          real(gasprod(i,2))
        end do
!        do i = 1, NZP1
!          write(lunf,'(i3,i4,3e11.3)') &
!          j, i, real(gflux(i,1)), real(gflux(i,2)), &
!          real(pflux(i,j,1)*rmass(j,1))
!          write(lunfp,'(i3,i4,3e11.3)') &
!          j, i, real(gflux(i,1)), real(gflux(i,2)), &
!          real(pflux(i,j,1)*rmass(j,1))
!        end do
      end do

    !endif


    if (idiag .eq. 1) then
      do iz = 1, NZ
        totmass(iz) = sum((mmr(iz,1,:)+mmr(iz,2,:)) * abs((pl(iz+1) - pl(iz))) / 88.7_f) / deltaz / 100._f
      end do
      write(lundiagn,*) 'PART 6: TOTAL MASS AT END OF TIME STEP'
      write(lundiagn,*) '****************************************'
      write(lundiagn,*) 'Z: ALTITUDE LEVEL INDEX'
      write(lundiagn,*) 'PARTMASS: TOTAL MASS DENSITY OF GAS AT Z [g/cm3]'
      write(lundiagn,*) 'GASMASS: TOTAL MASS DENSITY OF PARTICLES AT Z [g/cm3]'
      write(lundiagn,*) '****************************************'
      write(lundiagn,'(A6,2A15)') 'Z', 'PARTMASS', 'GASMASS'
      do iz = 1, NZ
	write(lundiagn,'(i6,2e15.3)') iz, totmass(iz), (mmr_gas(iz,1)+mmr_gas(iz,2)) * abs((pl(iz+1) - pl(iz))) / 88.7_f / deltaz / 100._f
      end do
      write(lundiagn,*) ' '
      write(lundiagn,'(A6,2e15.3)') 'TOT:',sum(totmass), &
                                    sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f / deltaz / 100._f )
      write(lundiagn,*) ' '
      write(lundiagn,*) 'COLUMN TOTAL MASS = TOT * NZ * deltaz (in meters) * (100 cm / m)'
      endcd = sum(totmass) * NZ * deltaz * 100._f + sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f ) * NZ
      write(lundiagn,'(A28,e28.15)') 'COLUMN TOTAL MASS [g/cm2]: ', endcd
      write(lundiagn,*) ' '
      write(lundiagn,*) 'PART 7: MASS CONSERVATION SUMMARY'
      write(lundiagn,*) '******************************************************************************'
      write(lundiagn,*) ''
      write(lundiagn,'(A62)') 'MASS CONSERVATION: Ab - Aa = (B + C + D) * E'
      write(lundiagn,*) ''
      write(lundiagn,'(A45,e28.15)') 'Aa. TOTAL COLUMN DENSTY AT START (g/cm2):', startcd
      write(lundiagn,'(A45,e28.15)') 'Ab. TOTAL COLUMN DENSTY AT END (g/cm2):', endcd
      write(lundiagn,'(A45,e28.15)') 'B. TOTAL INPUT RATE (g/cm2/s):', inputrate
      write(lundiagn,'(A45,e28.15)') 'C. TOTAL VERTICAL PARTICLE FLUX (g/cm2/s):', vertpartflux
      write(lundiagn,'(A45,e28.15)') 'D. TOTAL VERTICAL GAS FLUX (g/cm2/s):', vertgasflux
      write(lundiagn,'(A45,e28.15)') 'E. TIME STEP (s):', dtime
      write(lundiagn,*) ''
      write(lundiagn,'(A45,e28.15)') 'Ab - Aa = ', endcd - startcd
      write(lundiagn,'(A45,e28.15)') '(B + C + D) * E = ', (inputrate + vertgasflux + vertpartflux) * dtime
      write(lundiagn,'(A45,e28.15)') 'Ab - Aa -[(B + C + D) * E] = ', endcd - startcd - &
	[(inputrate + vertpartflux + vertgasflux) * dtime]
      write(lundiagn,*) ''
      write(lundiagn,*) '************************************************************10000000******************'
      write(lundiagn,*) ''
      write(lundiagn,*) ''
    end if


    !if (istep .ge. 89000) then
       write(lunres) istep+1, xc, dx, yc, dy, &
         zc, zl, p, pl, t, rho_atm_cgs, &
         mmr, mmr_gas, mmr_gas_old, satliq, &
         satice, satliq_old, satice_old, svpliq, wtpct, &
         zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
	 prodrate, prodrate_mass, prodrate_gas, totmass, &
	 winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
	 gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
         evappe, growlg, evaplg
       write(*,*)'write restart file'
       rewind(lunres)
    !endif
    endif

    close(unit=lunp)
    close(unit=lunfp)
    close(unit=lunratesp)
    close(unit=lunres)

  end do   ! time loop

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc < 0) stop "    *** FAILED CARMASTATE_Destroy ***"

  ! Close the output file
  close(unit=lun)
  close(unit=lunf)
  close(unit=lunrates)
  if (idiag .eq. 1) then
    close(unit=lundiagn)
  end if

  if (rc < 0) stop "    *** FAILED ***"

  write(*,*)  ""
  write(*,*) "CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc < 0) stop "    *** FAILED CARMA_Destroy ***"

end subroutine
