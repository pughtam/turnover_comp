program pft_phenology_v2
!------------------------------------------------------------------------------------------
! This program classifies simulated plant functional types (PFTs) into forest types using
! simulated LAI (1985-2014 30-year mean) and PFT phenological traits. See Pugh et al.
! (in review, Biogeosciences Discuss.) for additional details.
!
! Author: Sarah L. Shafer, U.S. Geological Survey (sshafer@usgs.gov)
!
! Version: 2.0
! Last update: 2020-02-06
!
! Note: This code was designed for use on workstations with 192 GB of
! memory and may need to be modified to run on other systems. NetCDF
! software is required to run this code
! (https://www.unidata.ucar.edu/software/netcdf/)
!
! USGS Disclaimer: This software is being provided to meet the need for timely best science.
! No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the
! functionality of the software and related material nor shall the fact of release
! constitute any such warranty. The software is provided on the condition that neither the
! USGS nor the U.S. Government shall be held liable for any damages resulting from the
! authorized or unauthorized use of the software. The USGS or the U.S. Government shall not
! be held liable for improper or incorrect use of the data described and/or contained
! herein. Any use of trade, firm, or product names is for descriptive purposes only and
! does not imply endorsement by the U.S. Government.
!------------------------------------------------------------------------------------------

use netCDF_subs ! Subroutines for reading and writing netCDF files.

implicit none

integer, parameter      :: nd = 4, nphen = 10, nphen_trees = 7
real(8), allocatable    :: clon(:),clat(:)
real(8)                 :: fillvalue,fillvalue_input,missingvalue,scalefactor,addoffset
real(8)                 :: lai_pass(5),lai_max,fillvalue_dble,pft_lai_mean_count(7)

real(8), allocatable    :: lai_month_model(:,:,:,:),lai_month(:,:,:,:)
real(8), allocatable    :: lai_annual(:,:,:,:),total_lai_annual(:,:,:)
real(8), allocatable    :: lai(:,:,:,:),fpc(:,:,:,:)
real(8), allocatable    :: total_lai(:,:,:),total_tree_lai(:,:,:),tree_max_lai(:,:,:)
real(8), allocatable    :: boreal_tree_lai(:,:,:),temperate_tree_lai(:,:,:)
real(8), allocatable    :: tropical_tree_lai(:,:,:),grass_lai(:,:,:),time_months(:),time(:)
real(8), allocatable    :: phen_total_lai(:,:,:,:)
real(8), allocatable    :: pft_max_lai(:,:,:),phen_max_lai(:,:,:),phen_lai_annual(:,:,:,:)
real(8), allocatable    :: pft_lai_to_total_lai_ratio(:,:,:,:),pft_max_lai_to_total_lai_ratio(:,:,:)
real(8), allocatable    :: pft_lai_to_total_tree_lai_ratio(:,:,:,:),pft_max_lai_to_total_tree_lai_ratio(:,:,:)
real(8), allocatable    :: phen_lai_to_total_lai_ratio(:,:,:,:),phen_max_lai_to_total_lai_ratio(:,:,:)
real(8), allocatable    :: phen_lai_to_total_tree_lai_ratio(:,:,:,:),phen_max_lai_to_total_tree_lai_ratio(:,:,:)
real(8), allocatable    :: pft_lai_mean(:,:,:,:),phen_lai_mean(:,:,:,:),total_lai_mean(:,:,:),total_tree_lai_mean(:,:,:)

real(4)                 :: fillvalue_input_real,missingvalue_real,scalefactor_real,addoffset_real,filltemp
real(4), allocatable    :: lai_month_float(:,:,:,:),height(:,:),fpc_float(:,:,:,:),lai_month_temp(:,:,:,:)
real(4), allocatable    :: lai_annual_float(:,:,:,:),landcover_frac(:,:,:,:)

integer, allocatable    :: biome(:,:,:),forest(:,:,:),mask(:,:),veg(:,:)
integer, allocatable    :: tree_max_lai_pft(:,:,:)
integer, allocatable    :: pft_max_lai_pftnum(:,:,:),phen_max_lai_phennum(:,:,:)
integer, allocatable    :: pft_phenclass_crosswalk(:),pft(:)

integer                 :: modelnum,nveg,ntime,nyrs,nmonths,dimid,pftdimid,year,fill_check,phen(nphen)
integer                 :: t_mean,phen_trees(nphen_trees),phentreesdimid,debug,nan_count(7)

integer                 :: pft_id,phen_id,biome_id,forest_id,total_lai_id,total_tree_lai_id
integer                 :: tree_max_lai_id,tree_max_lai_pft_id,boreal_tree_lai_id,temperate_tree_lai_id
integer                 :: tropical_tree_lai_id,grass_lai_id,lai_id,tree_lai_id,phentrees_id,nclon,nclat
integer                 :: forest_count,forest20_count,grass10,grass20,grass30
integer                 :: forest_count30lai,forest_count35lai,forest_count40lai
integer                 :: skip,totlai_id
integer                 :: phendimid
integer                 :: pft_max_lai_id,pft_max_lai_pftnum_id,pft_lai_ratio_id,pft_max_lai_ratio_id
integer                 :: pft_lai_id,phen_lai_id
integer                 :: phen_max_lai_id,phen_max_lai_phennum_id,phen_lai_ratio_id,phen_max_lai_ratio_id
integer                 :: pft_max_tree_lai_ratio_id,pft_tree_lai_ratio_id
integer                 :: phen_max_tree_lai_ratio_id,phen_tree_lai_ratio_id
 
integer                 :: i,m,p,t,x,y,z,fillvalue_int,vegnum,ntree,begin_veg

! netCDF compression variables
integer                 :: dflevel
logical                 :: shuffle_filter
integer                 :: chunk3d(3),chunk4d(4)

character(3)            :: var
character(20)           :: year_txt
character(19)           :: time_stamp
character(64)           :: veg_name,lai_name,latname,lonname

character(1024)         :: ncpath,ncoutpath,infile,maskfile,outfile,outfile_mask,fpcfile,vegfile,landcoverfile
character(1024)         :: pft_names,phen_names,tree_phen_names

phen = (/1,2,3,4,5,6,7,8,9,10/)
phen_trees = (/1,2,3,4,5,6,7/)

! netCDF compression parameters
dflevel = 3
shuffle_filter = "True"

ncpath = "\projects\mortality\data\"
ncoutpath = "\projects\mortality\"

latname = "lat"
lonname = "lon"
fillvalue = -1.0E32
fillvalue_dble = -1.0D32
fillvalue_int = -9999
pft_phenclass_crosswalk = -9999

! Phenology class array
phen_names = "1=needleleaved evergreen, &
2=needleleaved deciduous, &
3=boreal broadleaved deciduous, &
4=temperate broadleaved evergreen, &
5=temperate broadleaved deciduous, &
6=tropical broadleaved evergreen, &
7=tropical broadleaved raingreen, &
8=grass, &
9=shrub, &
10=tundra"

tree_phen_names = "1=needleleaved evergreen, &
2=needleleaved deciduous, &
3=boreal broadleaved deciduous, &
4=temperate broadleaved evergreen, &
5=temperate broadleaved deciduous, &
6=tropical broadleaved evergreen, &
7=tropical broadleaved raingreen"

! Models
! modelnum = 0  ! SEIB
! modelnum = 1  ! JULES
! modelnum = 2  ! LPJmL
! modelnum = 3  ! LPJ-GUESS
! modelnum = 4  ! LPJ-wsl
! modelnum = 5  ! CABLE-POP
! modelnum = 6  ! ORCHIDEE

do modelnum = 0,6
    begin_veg = 1
    debug = 0
    if (modelnum.eq.6) begin_veg = 2

! SEIB
    if (modelnum.eq.0) then
        infile = trim(ncpath)//"seib_cruncep_lai_month_1901_2014.nc4"
        maskfile = trim(ncoutpath)//"seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"seib_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"seib_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        veg_name = "vegtype"
        lai_name = "lai"
        latname = "latitude"
        lonname = "longitude"
        source = "seib_cruncep_lai_month_1901_2014.nc4"
        nclat = 360
        nclon = 720
        ntree = 12
        pft_names = "1=Tropical broad-leaved evergreen tree (Canopy, Shade tolerant), &
        2=Tropical broad-leaved evergreen tree (sub-Canopy, Shade tolerant), &
        3=Tropical broad-leaved evergreen tree (Pioneer, Light demanding), &
        4=Tropical broad-leaved evergreen tree (Short, Intermediate shade tolerant), &
        5=Tropical broad-leaved evergreen tree (Only distributes on Africa), &
        6=Tropical broad-leaved raingreen tree (Only distributes on Africa), &
        7=Temperate needle-leaved evergreen tree, &
        8=Temperate broad-leaved evergreen tree, &
        9=Temperate broad-leaved summergreen tree, &
        10=Boreal needle-leaved evergreen tree, &
        11=Boreal needle-leaved summergreen tree, &
        12=Boreal broad-leaved summergreen tree, &
        13=Temperate herbaceous (C3), &
        14=Temperate herbaceous (C3)"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(14))
        pft_phenclass_crosswalk(1) = 6
        pft_phenclass_crosswalk(2) = 6
        pft_phenclass_crosswalk(3) = 6
        pft_phenclass_crosswalk(4) = 6
        pft_phenclass_crosswalk(5) = 6
        pft_phenclass_crosswalk(6) = 7
        pft_phenclass_crosswalk(7) = 1
        pft_phenclass_crosswalk(8) = 4
        pft_phenclass_crosswalk(9) = 5
        pft_phenclass_crosswalk(10) = 1
        pft_phenclass_crosswalk(11) = 2
        pft_phenclass_crosswalk(12) = 3
        pft_phenclass_crosswalk(13) = 8
        pft_phenclass_crosswalk(14) = 8

! JULES
    else if (modelnum.eq.1) then
        infile = trim(ncpath)//"JULESC2_cruncep_lai_Monthly_1901_2014.nc"
        maskfile = trim(ncoutpath)//"JULESC2_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"JULESC2_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"JULESC2_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        veg_name = "pft"
        lai_name = "lai"
        source = "JULESC2_cruncep_lai_Monthly_1901_2014.nc"
        latname = "latitude"
        lonname = "longitude"
        nclat = 112
        nclon = 192
        ntree = 5
        pft_names = "1=Tropical broadleaved evergreen tree, &
        2=Temperate broadleaved evergreen tree, &
        3=Broadleaf deciduous tree, &
        4=Needleleaf evergreen tree, &
        5=Needleleaf deciduous tree, &
        6=C3 grass, &
        7=C4 grass, &
        8=evergreen shrub, &
        9=deciduous shrub"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(9))
        pft_phenclass_crosswalk(1) = 6
        pft_phenclass_crosswalk(2) = 4
        pft_phenclass_crosswalk(3) = -5 ! Broadleaf summergreen; Divide into boreal,temperate,tropical zones
        pft_phenclass_crosswalk(4) = 1
        pft_phenclass_crosswalk(5) = 2
        pft_phenclass_crosswalk(6) = 8
        pft_phenclass_crosswalk(7) = 8
        pft_phenclass_crosswalk(8) = 9
        pft_phenclass_crosswalk(9) = 9

! LPJmL
    else if (modelnum.eq.2) then
        infile = trim(ncpath)//"lpjml_cruncep_lai_annual_1901_2014.nc"
        maskfile = trim(ncoutpath)//"lpjml_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"lpjml_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"lpjml_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        veg_name = "npft"
        lai_name = "lai"
        source = "lpjml_cruncep_lai_annual_1901_2014.nc"
        latname = "latitude"
        lonname = "longitude"
        nclat = 279
        nclon = 720
        ntree = 7
        pft_names = "1=Tropical broad-leaved evergreen tree, &
        2=Tropical broad-leaved raingreen tree, &
        3=Temperate needle-leaved evergreen tree, &
        4=Temperate broad-leaved evergreen tree, &
        5=Temperate broad-leaved summergreen tree, &
        6=Boreal needle-leaved evergreen tree, &
        7=Boreal broad-leaved summergreen tree, &
        8=Temperate herbaceous, &
        9=Tropical herbaceous"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(9))
        pft_phenclass_crosswalk(1) = 6
        pft_phenclass_crosswalk(2) = 7
        pft_phenclass_crosswalk(3) = 1
        pft_phenclass_crosswalk(4) = 4
        pft_phenclass_crosswalk(5) = 5
        pft_phenclass_crosswalk(6) = 1
        pft_phenclass_crosswalk(7) = 3
        pft_phenclass_crosswalk(8) = 8
        pft_phenclass_crosswalk(9) = 8

! LPJ-GUESS
    else if (modelnum.eq.3) then
        infile = trim(ncpath)//"lpj-guess_cruncep_lai_monthly_1901_2014.nc4"
        maskfile = trim(ncoutpath)//"lpj-guess_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"lpj-guess_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"lpj-guess_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        veg_name = "vegtype"
        lai_name = "lai"
        latname = "latitude"
        lonname = "longitude"
        source = "lpj-guess_cruncep_lai_monthly_1901_2014.nc4"
        nclat = 360
        nclon = 720
        ntree = 9
        pft_names = "1=Boreal needleleaf evergreen, &
        2=Boreal shade-intolerant needleleaf evergreen, &
        3=Boreal needleleaf summergreen, &
        4=Temperate broadleaf summergreen, &
        5=Temperate shade-intolerant broadleaf summergreen, &
        6=Temperate broadleaf evergreen, &
        7=Tropical broadleaf evergreen, &
        8=Tropical shade-intolerant broadleaf evergreen, &
        9=Tropical broadleaf raingreen, &
        10=C3 grasses, &
        11=C4 grasses"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(11))
        pft_phenclass_crosswalk(1) = 1
        pft_phenclass_crosswalk(2) = 1
        pft_phenclass_crosswalk(3) = 2
        pft_phenclass_crosswalk(4) = 5
        pft_phenclass_crosswalk(5) = 5
        pft_phenclass_crosswalk(6) = 4
        pft_phenclass_crosswalk(7) = 6
        pft_phenclass_crosswalk(8) = 6
        pft_phenclass_crosswalk(9) = 7
        pft_phenclass_crosswalk(10) = 8
        pft_phenclass_crosswalk(11) = 8

! LPJ-wsl
    else if (modelnum.eq.4) then
        infile = trim(ncpath)//"LPJ_pft_mlai_BONE.nc"
        landcoverfile = trim(ncpath)//"LPJ_landCoverFrac.nc"
        maskfile = trim(ncoutpath)//"lpj-wsl_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"lpj-wsl_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"lpj-wsl_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        lai_name = "lai"
        latname = "lat"
        lonname = "lon"
        source = "LPJ_pft_mlai_[pft name].nc, LPJ_landCoverFrac.nc"
        nclat = 360
        nclon = 720
        ntree = 7
        pft_names = "1=Tropical evergreen (TREV), &
        2=Tropical raingreen (TRRG), &
        3=Temperate needleleaf evergreen (TENE), &
        4=Temperate broadleaf evergreen (TEBE), &
        5=Temperate broadleaf summergreen (TEBS), &
        6=Boreal needleleaf evergreen (BONE), &
        7=Boreal broadleaf summergreen (BOBS), &
        8=C3, &
        9=C4"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(9))
        pft_phenclass_crosswalk(1) = 6
        pft_phenclass_crosswalk(2) = 7
        pft_phenclass_crosswalk(3) = 1
        pft_phenclass_crosswalk(4) = 4
        pft_phenclass_crosswalk(5) = 5
        pft_phenclass_crosswalk(6) = 1
        pft_phenclass_crosswalk(7) = 3
        pft_phenclass_crosswalk(8) = 8
        pft_phenclass_crosswalk(9) = 8

! CABLE-POP
    else if (modelnum.eq.5) then
        infile = trim(ncpath)//"CABLE-POP_cruncep_lai_month_1901_2015.nc4"
        maskfile = trim(ncoutpath)//"CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"CABLE-POP_cruncep_lai_annual_1901_2015_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"CABLE-POP_cruncep_lai_annual_1901_2015_phenology_mask.nc"

        lai_name = "lai"
        latname = "y"
        lonname = "x"
        source = "CABLE-POP_cruncep_lai_month_1901_2015.nc"
        nclat = 360
        nclon = 720
        ntree = 4
        pft_names = "1=Evergreen Needleleaf Forest, &
        2=Evergreen Broadleaf Forest, &
        3=Deciduous Needleleaf Forest, &
        4=Deciduous Broadleaf Forest, &
        5=shrub, &
        6=C3 grass, &
        7=C4 grass, &
        8=tundra"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(8))
        pft_phenclass_crosswalk(1) = 1
        pft_phenclass_crosswalk(2) = -7 ! Broadleaf evergreen, divide into temperate and tropical
        pft_phenclass_crosswalk(3) = 2
        pft_phenclass_crosswalk(4) = -5 ! Deciduous broadleaf forest, divide into boreal, temperate, and tropical
        pft_phenclass_crosswalk(5) = 9
        pft_phenclass_crosswalk(6) = 8
        pft_phenclass_crosswalk(7) = 8
        pft_phenclass_crosswalk(8) = 10

! ORCHIDEE
    else if (modelnum.eq.6) then
        infile = trim(ncpath)//"ORCHIDEEr3085_cruncep_LAI_13pft_year_1901_2014.nc"
        fpcfile = trim(ncpath)//"ORCHIDEEr3085_cruncep_veget_max_13pft_year_1901_2014.nc"
        maskfile = trim(ncoutpath)//"ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
        outfile = trim(ncoutpath)//"ORCHIDEE_cruncep_lai_annual_1901_2014_phenology_draft3.nc"
        outfile_mask = trim(ncoutpath)//"ORCHIDEE_cruncep_lai_annual_1901_2014_phenology_mask.nc"

        lai_name = "lai"
        latname = "lat"
        lonname = "lon"
        source = "ORCHIDEEr3085_cruncep_LAI_13pft_year_1901_2014.nc, ORCHIDEEr3085_cruncep_veget_max_13pft_year_1901_2014.nc"
        nclat = 360
        nclon = 720
        ntree = 9
        pft_names = "1=Bare ground, &
        2=tropical broad-leaved evergreen, &
        3=tropical broad-leaved raingreen, &
        4=temperate needleleaf evergreen, &
        5=temperate broad-leaved evergreen, &
        6=temperate broad-leaved summergreen, &
        7=boreal needleleaf evergreeen, &
        8=boreal broad-leaved summergreen, &
        9=boreal needleleaf summergreen, &
        10=C3 grass, &
        11=C4 grass"

        ! Model PFT to phenology class crosswalk
        allocate (pft_phenclass_crosswalk(11))
        pft_phenclass_crosswalk(1) = fillvalue_int
        pft_phenclass_crosswalk(2) = 6
        pft_phenclass_crosswalk(3) = 7
        pft_phenclass_crosswalk(4) = 1
        pft_phenclass_crosswalk(5) = 4
        pft_phenclass_crosswalk(6) = 5
        pft_phenclass_crosswalk(7) = 1
        pft_phenclass_crosswalk(8) = 3
        pft_phenclass_crosswalk(9) = 2
        pft_phenclass_crosswalk(10) = 8
        pft_phenclass_crosswalk(11) = 8
    end if

    allocate (clat(nclat),clon(nclon))

    write (*,*) "Reading input data..."

    call netcdf_open(infile,0,ncid)
    call netcdf_get_1d_double(ncid,nclon,trim(lonname),clon)
    call netcdf_get_1d_double(ncid,nclat,trim(latname),clat)

    if (modelnum.lt.6) then
        call netcdf_get_dimlen(ncid,"time",nmonths)
        allocate (time_months(nmonths))
        call netcdf_get_1d_double(ncid,nmonths,"time",time_months)
    else if (modelnum.eq.6) then
        call netcdf_get_dimlen(ncid,"time_counter",ntime) ! Annual data
        allocate (time(ntime))
        call netcdf_get_1d_double(ncid,ntime,"time_counter",time)
    endif

    if (modelnum.lt.4) call netcdf_get_dimlen(ncid,veg_name,nveg)
    if (modelnum.eq.4) nveg = 9  ! LPJ-wsl
    if (modelnum.eq.5) nveg = 8  ! CABLE-POP
    if (modelnum.eq.6) nveg = 11 ! ORCHIDEE

    allocate (pft(nveg))
    do i = 1,nveg
        pft(i) = i
    end do

    if (modelnum.eq.0) then ! SEIB
        allocate (lai_month_model(nclon,nclat,nmonths,nveg))
        allocate (lai_month(nclon,nclat,nveg,nmonths))
        call netcdf_get_4d_double(ncid,nclon,nclat,nmonths,nveg,trim(lai_name),lai_month_model)
        call netcdf_close(infile,ncid)

        ! Restructure array
        write (*,*) "Restructuring LAI array..."
        do x = 1,nclon
            write (*,*) "Restructuring longitude = ",x
            do y = 1,nclat
                do z = 1,nveg
                    do t = 1,nmonths
                        lai_month(x,y,z,t) = lai_month_model(x,y,t,z)
                    end do
                end do
            end do
        end do
        deallocate (lai_month_model)

    else if (modelnum.eq.1) then ! JULES
        allocate (lai_month_float(nclon,nclat,nveg,nmonths))
        allocate (lai_month(nclon,nclat,nveg,nmonths))
        call netcdf_get_4d_real(ncid,nclon,nclat,nveg,nmonths,"lai",lai_month_float)
        call netcdf_close(infile,ncid)

        do x = 1,nclon
            do y = 1,nclat
                do z = 1,nveg
                    do t = 1,nmonths
                        if (lai_month_float(x,y,z,t).ge.0.00000000) then
                            lai_month(x,y,z,t) = dble(lai_month_float(x,y,z,t))
                        else
                            lai_month(x,y,z,t) = fillvalue_dble
                        endif
                    end do
                end do
           end do
       end do
       deallocate (lai_month_float)

    else if (modelnum.eq.2) then ! LPJmL
        allocate (lai_month_float(nclon,nclat,nveg,nmonths))
        allocate (lai_annual(nclon,nclat,nveg,nmonths))
        lai_annual = 0.0d0
        call netcdf_get_4d_real(ncid,nclon,nclat,nveg,nmonths,"lai",lai_month_float)
        call netcdf_close(infile,ncid)

        do x = 1,nclon
            do y = 1,nclat
                do z = 1,nveg
                    do t = 1,nmonths
                        lai_annual(x,y,z,t) = dble(lai_month_float(x,y,z,t))
                    end do
                end do
            end do
        end do
        deallocate (lai_month_float)

    else if (modelnum.eq.3) then ! LPJ-GUESS
        allocate (lai_month(nclon,nclat,nveg,nmonths))
        call netcdf_get_4d_double(ncid,nclon,nclat,nveg,nmonths,trim(lai_name),lai_month)
        call netcdf_close(infile,ncid)

    else if (modelnum.eq.4) then ! LPJ-wsl
        allocate (lai_month(nclon,nclat,nveg,nmonths))
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,6,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_BOBS.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,7,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_TENE.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,3,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_TEBS.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,5,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_TEBE.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,4,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_TREV.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,1,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_TRRG.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,2,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_C3.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,8,:))
        call netcdf_close(infile,ncid)

        infile = trim(ncpath)//"LPJ_pft_mlai_C4.nc"
        call netcdf_open(infile,0,ncid)
            call netcdf_get_3d_double(ncid,nclon,nclat,nmonths,"lai",lai_month(:,:,9,:))
        call netcdf_close(infile,ncid)

    else if (modelnum.eq.5) then ! CABLE-POP
        allocate (lai_month(nclon,nclat,nveg,nmonths),lai_month_temp(nclon,nclat,nveg,nmonths))
        call netcdf_get_4d_real(ncid,nclon,nclat,nveg,nmonths,trim(lai_name),lai_month_temp)
        call netcdf_get_att_real(ncid,trim(lai_name),"_FillValue",filltemp)
        call netcdf_close(infile,ncid)

        write (*,*) filltemp

        do y = 1,nclat
            do x = 1,nclon
                do z = 1,nveg
                     do t = 1,nmonths
                        if (lai_month_temp(x,y,z,t).ne.filltemp) then
                            lai_month(x,y,z,t) = dble(lai_month_temp(x,y,z,t))
                        else
                            lai_month(x,y,z,t) = fillvalue_dble
                        end if
                    end do
                end do
            end do
        end do
        deallocate (lai_month_temp)

    else if (modelnum.eq.6) then ! ORCHIDEE
        allocate (lai_annual_float(nclon,nclat,nveg,ntime))
        allocate (lai_annual(nclon,nclat,nveg,ntime))
        allocate (fpc_float(nclon,nclat,nveg,ntime))
        lai_annual = 0.0d0

        call netcdf_get_4d_real(ncid,nclon,nclat,nveg,ntime,trim(lai_name),lai_annual_float)
        call netcdf_get_att_real(ncid,trim(lai_name),"_FillValue",filltemp)
        call netcdf_close(infile,ncid)

        call netcdf_open(fpcfile,0,ncid) !veget_max = Maximum vegetation fraction (LAI -> infinity)
            call netcdf_get_4d_real(ncid,nclon,nclat,nveg,ntime,"veget_max",fpc_float)
        call netcdf_close(fpcfile,ncid)

        do x = 1,nclon
            do y = 1,nclat
                do z = 1,nveg
                    do t = 1,ntime
                        ! A fill value of -99999.0f is assigned to the file, but only used for year 1. After year 1,
                        ! ocean and other non-modeled points have a value >1.0e30.
                        if (lai_annual_float(x,y,z,t).ne.filltemp.and.lai_annual_float(x,y,z,t).le.999999.0) then
                            lai_annual(x,y,z,t) = dble(lai_annual_float(x,y,z,t) * fpc_float(x,y,z,t))
                        else if (lai_annual_float(x,y,z,t).eq.filltemp.or.lai_annual_float(x,y,z,t).gt.999999.0) then
                            lai_annual(x,y,z,t) = fillvalue_dble
                        endif
                    end do
                end do
            end do
        end do
        deallocate (lai_annual_float,fpc_float)
    endif

    ! SEIB
    if (modelnum.eq.0) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! JULES
    else if (modelnum.eq.1) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! LPJmL
    else if (modelnum.eq.2) then
        ntime = nmonths
        allocate (time(ntime))

    ! LPJ-GUESS
    else if (modelnum.eq.3) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! LPJ-wsl
    else if (modelnum.eq.4) then
        ntime = nmonths / 12
        allocate (time(ntime))

        allocate (landcover_frac(nclon,nclat,ntime,nveg))

        call netcdf_open(landcoverfile,0,ncid)
            call netcdf_get_4d_real(ncid,nclon,nclat,ntime,nveg,"landCoverFrac",landcover_frac)
        call netcdf_close(landcoverfile,ncid)

    ! CABLE-POP
    else if (modelnum.eq.5) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! ORCHIDEE
    else if (modelnum.eq.6) then
        ! Intentionally blank
    end if

    do t = 1,ntime
        time(t) = (dble(t) * 1.0d0) + 1900.0d0
    end do

    if (modelnum.ne.2.and.modelnum.ne.6) then
        allocate (lai_annual(nclon,nclat,nveg,ntime))
        lai_annual = 0.0d0
    end if

    if (modelnum.ne.2.and.modelnum.ne.6) then
        do y = 1,nclat
            if (mod(y,100).eq.0) write (*,*) "row = ",y
            do x = 1,nclon
                year = 0
                do t = 1,nmonths,12
                    year = year + 1
                    do i = 1,nveg
                        fill_check = 0
                        lai_max = -99999.0d0
                        do m = 0,11
                            if (lai_month(x,y,i,t+m).ge.0.0d0) then
                                if (lai_month(x,y,i,t+m).ge.0.00000000d0) fill_check = 1
                                if (lai_month(x,y,i,t+m).gt.lai_max) lai_max = lai_month(x,y,i,t+m)
                            end if
                        end do
                        lai_annual(x,y,i,year) = lai_max
                        if (fill_check.eq.0) lai_annual(x,y,i,year) = fillvalue_dble
                    end do
                end do
            end do
        end do
        deallocate (lai_month)
    else if (modelnum.eq.2) then ! LPJmL
        do y = 1,nclat
            if (mod(y,100).eq.0) write (*,*) "row = ",y
            do x = 1,nclon
                year = 0
                do t = 1,nmonths
                    do i = 1,nveg
                        ! LPJmL fill value = -9999.
                        if (lai_annual(x,y,i,t).lt.-9990.0d0) lai_annual(x,y,i,t) = fillvalue_dble
                    end do
                end do
            end do
        end do
    else if (modelnum.eq.6) then ! ORCHIDEE
        ! For ORCHIDEE data no additional modification needed.
    endif

    allocate (mask(nclon,nclat))
    call netcdf_open(maskfile,0,ncid)
        call netcdf_get_2d_int(ncid,nclon,nclat,"forest_30yr_any_10_years",mask)
    call netcdf_close(maskfile,ncid)

    allocate (total_tree_lai(nclon,nclat,ntime),tree_max_lai(nclon,nclat,ntime))
    allocate (tree_max_lai_pft(nclon,nclat,ntime))
    allocate (phen_lai_annual(nclon,nclat,nphen,ntime),total_lai_annual(nclon,nclat,ntime))
    allocate (pft_lai_mean(nclon,nclat,nveg,ntime),phen_lai_mean(nclon,nclat,nphen_trees,ntime))
    allocate (total_lai_mean(nclon,nclat,ntime),total_tree_lai_mean(nclon,nclat,ntime))
    allocate (pft_max_lai(nclon,nclat,ntime),pft_max_lai_pftnum(nclon,nclat,ntime))
    allocate (phen_max_lai(nclon,nclat,ntime),phen_max_lai_phennum(nclon,nclat,ntime))
    allocate (pft_lai_to_total_lai_ratio(nclon,nclat,nveg,ntime))
    allocate (pft_max_lai_to_total_lai_ratio(nclon,nclat,ntime))
    allocate (phen_lai_to_total_lai_ratio(nclon,nclat,nphen_trees,ntime))
    allocate (phen_max_lai_to_total_lai_ratio(nclon,nclat,ntime))
    allocate (pft_lai_to_total_tree_lai_ratio(nclon,nclat,nveg,ntime))
    allocate (pft_max_lai_to_total_tree_lai_ratio(nclon,nclat,ntime))
    allocate (phen_lai_to_total_tree_lai_ratio(nclon,nclat,nphen_trees,ntime))
    allocate (phen_max_lai_to_total_tree_lai_ratio(nclon,nclat,ntime))

    total_tree_lai = 0.0d0
    total_lai_annual = 0.0d0
    phen_lai_annual = 0.0d0
    total_lai_mean = 0.0d0
    total_tree_lai_mean = 0.0d0
    pft_lai_mean = 0.0d0
    phen_lai_mean = 0.0d0

    tree_max_lai = fillvalue_dble
    tree_max_lai_pft = fillvalue_int

    pft_max_lai = -9999.0d0
    pft_max_lai_pftnum = fillvalue_int

    p = fillvalue_int
    phen_max_lai = -9999.0d0
    phen_max_lai_phennum = fillvalue_int

    pft_lai_to_total_lai_ratio = fillvalue_dble
    pft_max_lai_to_total_lai_ratio = fillvalue_dble

    phen_lai_to_total_lai_ratio = fillvalue_dble
    phen_max_lai_to_total_lai_ratio = fillvalue_dble

    pft_lai_to_total_tree_lai_ratio = fillvalue_dble
    pft_max_lai_to_total_tree_lai_ratio = fillvalue_dble

    phen_lai_to_total_tree_lai_ratio = fillvalue_dble
    phen_max_lai_to_total_tree_lai_ratio = fillvalue_dble

    write (*,*) "Calculating total lai..."
    do y = 1,nclat
        write (*,*) "          ",y
        do x = 1,nclon
            do t = 1,ntime

                if (debug.eq.1) write (*,*) x,y,t

                do z = begin_veg,nveg
                    p = pft_phenclass_crosswalk(z)
                    if (lai_annual(x,y,z,t).ge.0..and.mask(x,y).eq.1) then

                        if (modelnum.eq.0.or.modelnum.eq.1.or.modelnum.eq.2.or.modelnum.eq.3.or.modelnum.eq.5.or.modelnum.eq.6) then
                            if (lai_annual(x,y,z,t).gt.0.) then
                                total_lai_annual(x,y,t) = total_lai_annual(x,y,t) + lai_annual(x,y,z,t)

                                if (z.lt.(ntree+1)) total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            end if
                        else if (modelnum.eq.4)then
                            if (lai_annual(x,y,z,t).gt.0..and.landcover_frac(x,y,t,z).gt.0.) then
                                total_lai_annual(x,y,t) = total_lai_annual(x,y,t) + &
                                    (lai_annual(x,y,z,t) * dble(landcover_frac(x,y,t,z)))

                                if (z.lt.(ntree+1)) total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + &
                                    (lai_annual(x,y,z,t) * dble(landcover_frac(x,y,t,z)))
                            end if
                        endif

                        if (p.gt.0.and.p.lt.(nphen_trees+1)) then
                            if (modelnum.eq.0.or.modelnum.eq.1.or.modelnum.eq.2.or.modelnum.eq.3.or.modelnum.eq.5.or.modelnum.eq.6) then
                                phen_lai_annual(x,y,p,t) = phen_lai_annual(x,y,p,t) + lai_annual(x,y,z,t)
 
                            else if (modelnum.eq.4) then
                                phen_lai_annual(x,y,p,t) = phen_lai_annual(x,y,p,t) + &
                                    (lai_annual(x,y,z,t) * dble(landcover_frac(x,y,t,z)))
                            end if

                        else if (p.eq.-5.and.modelnum.eq.1) then ! JULES
                            if (clat(y).gt.55.0d0) then ! Boreal
                                phen_lai_annual(x,y,3,t) = phen_lai_annual(x,y,3,t) + lai_annual(x,y,3,t)
                            else if ((clat(y).le.55.0d0.and.clat(y).gt.23.0d0).or.(clat(y).lt.-23.0d0)) then ! Temperate
                                phen_lai_annual(x,y,5,t) = phen_lai_annual(x,y,5,t) + lai_annual(x,y,3,t)
                            else if (clat(y).le.23.0d0.and.clat(y).ge.-23.0d0) then ! Tropical
                                phen_lai_annual(x,y,7,t) = phen_lai_annual(x,y,7,t) + lai_annual(x,y,3,t)
                            end if

                        else if (p.eq.-5.and.modelnum.eq.5) then ! CABLE-POP deciduous broadleaf
                            if (clat(y).gt.55.0d0) then !Boreal
                                phen_lai_annual(x,y,3,t) = phen_lai_annual(x,y,3,t) + lai_annual(x,y,4,t)
                            else if ((clat(y).le.55.0d0.and.clat(y).gt.23.0d0).or.(clat(y).lt.-23.0d0)) then ! Temperate
                                phen_lai_annual(x,y,5,t) = phen_lai_annual(x,y,5,t) + lai_annual(x,y,4,t)
                            else if (clat(y).le.23.0d0.and.clat(y).ge.-23.0d0) then ! Tropical
                                phen_lai_annual(x,y,7,t) = phen_lai_annual(x,y,7,t) + lai_annual(x,y,4,t)
                            end if

                        else if (p.eq.-7.and.modelnum.eq.5) then ! CABLE-POP evergreen broadleaf PFT
                            if (clat(y).gt.55.0d0) then !Boreal
                                stop "CABLE-POP boreal evergreen broadleaf PFT"
                                phen_lai_annual(x,y,4,t) = phen_lai_annual(x,y,4,t) + lai_annual(x,y,2,t)
                            else if ((clat(y).le.55.0d0.and.clat(y).gt.23.0d0).or.(clat(y).lt.-23.0d0)) then ! Temperate
                                phen_lai_annual(x,y,4,t) = phen_lai_annual(x,y,4,t) + lai_annual(x,y,2,t)
                            else if (clat(y).le.23.0d0.and.clat(y).ge.-23.0d0) then ! Tropical
                                phen_lai_annual(x,y,6,t) = phen_lai_annual(x,y,6,t) + lai_annual(x,y,2,t)
                            end if

                        ! Trap errors in PFT x phenology class crosswalk
                        else if (p.le.0.and.modelnum.ne.1.and.modelnum.ne.5) then
                            stop "p.le.0"
                        else if (p.lt.-5.and.modelnum.eq.1) then
                            stop "p.lt.-5"
                        else if (p.lt.-7.and.modelnum.eq.5) then
                            stop "p.lt.-7"
                        endif

                    else if (mask(x,y).ne.1 &
                        .or.lai_annual(x,y,z,t).lt.0.or.(landcover_frac(x,y,t,z).le.-9990.and.modelnum.eq.4)) then
                        total_lai_annual(x,y,t) = fillvalue_dble
                        total_tree_lai(x,y,t) = fillvalue_dble
                        tree_max_lai(x,y,t) = fillvalue_dble
                        tree_max_lai_pft(x,y,t) = fillvalue_int
                        pft_max_lai(x,y,t) = fillvalue_dble
                        phen_max_lai(x,y,t) = fillvalue_dble
                        if (mask(x,y).ne.1) pft_max_lai_pftnum(x,y,t) = fillvalue_int
                        if (mask(x,y).ne.1) phen_max_lai_phennum(x,y,t) = fillvalue_int
                        if (p.gt.0.and.p.lt.(nphen_trees+1)) phen_lai_annual(x,y,p,t) = fillvalue_dble
                        if (p.eq.-5.and.(modelnum.eq.1.or.modelnum.eq.5)) then ! JULES and CABLE-POP
                            phen_lai_annual(x,y,3,t) = fillvalue_dble
                            phen_lai_annual(x,y,5,t) = fillvalue_dble
                            phen_lai_annual(x,y,7,t) = fillvalue_dble
                        else if (p.eq.-7.and.modelnum.eq.5) then ! CABLE-POP
                            phen_lai_annual(x,y,4,t) = fillvalue_dble
                            phen_lai_annual(x,y,6,t) = fillvalue_dble
                        end if
                    endif
                end do ! PFTs

                if (modelnum.eq.2.or.modelnum.eq.4) phen_lai_annual(x,y,2,t) = fillvalue_dble ! Phenology class 2 not simulated by LPJ-wsl or LPJml
                if (modelnum.eq.3) phen_lai_annual(x,y,3,t) = fillvalue_dble ! Phenology class 3 not simulated by LPJ-GUESS

                if (mask(x,y).eq.1) then ! Calculate phen_lai_mean
                    if (t.gt.29) then
                        nan_count(:) = 0
                        do t_mean = t-29,t
                            do p = 1,nphen_trees ! Phenology classes
                                if (modelnum.eq.0.or.modelnum.eq.1.or.modelnum.eq.5.or.modelnum.eq.6) then !SEIB, JULES, CABLE-POP, ORCHIDEE
                                    if (phen_lai_annual(x,y,p,t_mean).gt.0.) then
                                        phen_lai_mean(x,y,p,t) = phen_lai_mean(x,y,p,t) + &
                                            (phen_lai_annual(x,y,p,t_mean) / 30.0d0)
                                    else if (phen_lai_annual(x,y,p,t_mean).lt.0) then
                                        phen_lai_mean(x,y,p,t) = fillvalue_dble
                                    endif

                                else if (modelnum.eq.2.or.modelnum.eq.4) then ! LPJmL and LPJ-wsl
                                    if (p.ne.2.and.phen_lai_annual(x,y,p,t_mean).gt.0.) then
                                        phen_lai_mean(x,y,p,t) = phen_lai_mean(x,y,p,t) + &
                                            (phen_lai_annual(x,y,p,t_mean) / 30.0d0)
                                    else if (phen_lai_annual(x,y,p,t_mean).lt.0) then
                                        phen_lai_mean(x,y,p,t) = fillvalue_dble
                                    endif

                                else if (modelnum.eq.3) then ! LPJ-GUESS
                                    if (p.ne.3.and.phen_lai_annual(x,y,p,t_mean).gt.0.) then
                                        phen_lai_mean(x,y,p,t) = phen_lai_mean(x,y,p,t) + &
                                            (phen_lai_annual(x,y,p,t_mean) / 30.0d0)
                                    else if (phen_lai_annual(x,y,p,t_mean).lt.0) then
                                        phen_lai_mean(x,y,p,t) = fillvalue_dble
                                    endif
                                endif
                            end do

                            do z = begin_veg,nveg ! Tree PFTs
                                if (modelnum.eq.0.or.modelnum.eq.1.or.modelnum.eq.2.or.modelnum.eq.3.or.modelnum.eq.5.or.modelnum.eq.6) then

                                    if (lai_annual(x,y,z,t_mean).gt.0.) then
                                        pft_lai_mean(x,y,z,t) = pft_lai_mean(x,y,z,t) + (lai_annual(x,y,z,t_mean) / 30.0d0)
                                    end if

                                else if (modelnum.eq.4) then
                                    if (lai_annual(x,y,z,t_mean).gt.0. &
                                        .and.landcover_frac(x,y,t_mean,z).gt.-99990) then
                                        pft_lai_mean(x,y,z,t) = pft_lai_mean(x,y,z,t) + ((lai_annual(x,y,z,t_mean)*dble(landcover_frac(x,y,t_mean,z))) / 30.0d0)
                                    end if
                                endif
                            end do

                            if (total_lai_annual(x,y,t_mean).gt.0.) &
                                total_lai_mean(x,y,t) = total_lai_mean(x,y,t) + (total_lai_annual(x,y,t_mean) / 30.0d0)
                            if (total_tree_lai(x,y,t_mean).gt.0.) &
                                total_tree_lai_mean(x,y,t) = total_tree_lai_mean(x,y,t) + (total_tree_lai(x,y,t_mean) / 30.0d0)
                        end do
                    endif
                else
                    phen_lai_mean(x,y,:,t) = fillvalue_dble
                    pft_lai_mean(x,y,:,t) = fillvalue_dble
                    total_lai_mean(x,y,t) = fillvalue_dble
                    total_tree_lai_mean(x,y,t) = fillvalue_dble
                end if

                do z = begin_veg,ntree ! Tree PFTs
                    if (pft_lai_mean(x,y,z,t).ge.0. &
                        .and.mask(x,y).eq.1) then
                        if (pft_lai_mean(x,y,z,t).gt.pft_max_lai(x,y,t)) then
                            pft_max_lai(x,y,t) = pft_lai_mean(x,y,z,t)
                            pft_max_lai_pftnum(x,y,t) = z
                        endif

                    else if (pft_lai_mean(x,y,z,t).le.fillvalue_dble.or.mask(x,y).ne.1) then
                        pft_max_lai_pftnum(x,y,t) = fillvalue_int
                    end if
                end do

                do p = 1,nphen_trees ! Tree phenology classes
                    if (phen_lai_mean(x,y,p,t).ge.0. &
                        .and.mask(x,y).eq.1) then
                        if (phen_lai_mean(x,y,p,t).gt.phen_max_lai(x,y,t)) then
                            phen_max_lai(x,y,t) = phen_lai_mean(x,y,p,t)
                            phen_max_lai_phennum(x,y,t) = p
                        endif
                    endif
                end do
                ! Mask is for 1985-2014. If forest mask = 2 in earlier years, phen_lai_mean = fill
                ! so phen_max_lai value never changes.
                if (phen_max_lai(x,y,t).lt.-9990.) phen_max_lai(x,y,t) = fillvalue_dble

                do z = begin_veg,nveg ! PFT LAI / Total LAI
                    if (pft_lai_mean(x,y,z,t).gt.fillvalue_dble.and.total_lai_mean(x,y,t).gt.fillvalue_dble) then
                        if (pft_lai_mean(x,y,z,t).gt.0.0d0.and.total_lai_mean(x,y,t).gt.0.0d0) then
                            pft_lai_to_total_lai_ratio(x,y,z,t) = pft_lai_mean(x,y,z,t) / total_lai_mean(x,y,t)
                        else if (pft_lai_mean(x,y,z,t).eq.0.0d0.or.total_lai_mean(x,y,t).eq.0.0d0) then
                            pft_lai_to_total_lai_ratio(x,y,z,t) = 0.0d0
                        else if (pft_lai_mean(x,y,z,t).lt.0.0d0.or.total_lai_mean(x,y,t).lt.0.0d0) then
                            stop "PFT LAI / Total LAI"
                        end if
                    else
                        pft_lai_to_total_lai_ratio(x,y,z,t) = fillvalue_dble
                    end if
                end do

                do z = begin_veg,nveg ! PFT LAI / Total tree LAI
                    if (pft_lai_mean(x,y,z,t).gt.fillvalue_dble.and.total_tree_lai_mean(x,y,t).gt.fillvalue_dble) then
                        if (pft_lai_mean(x,y,z,t).gt.0.0d0.and.total_tree_lai_mean(x,y,t).gt.0.0d0) then
                            pft_lai_to_total_tree_lai_ratio(x,y,z,t) = pft_lai_mean(x,y,z,t) / total_tree_lai_mean(x,y,t)
                        else if (pft_lai_mean(x,y,z,t).eq.0.0d0.or.total_tree_lai_mean(x,y,t).eq.0.0d0) then
                            pft_lai_to_total_tree_lai_ratio(x,y,z,t) = 0.0d0
                        else if (pft_lai_mean(x,y,z,t).lt.0.0d0.or.total_tree_lai_mean(x,y,t).lt.0.0d0) then
                            stop "PFT LAI / Total tree LAI"
                        end if
                    else
                        pft_lai_to_total_tree_lai_ratio(x,y,z,t)=fillvalue_dble
                    end if
                end do

                ! PFT max LAI / Total LAI
                if (pft_max_lai_pftnum(x,y,t).gt.0.and.pft_max_lai_pftnum(x,y,t).lt.(ntree+1) &
                    .and.pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).gt.fillvalue_dble.and.total_lai_mean(x,y,t).gt.fillvalue_dble) then
                    if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).gt.0.0d0.and.total_lai_mean(x,y,t).gt.0.0d0) then
                        pft_max_lai_to_total_lai_ratio(x,y,t) = pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t) / total_lai_mean(x,y,t)

                    else if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).eq.0.0d0.or.total_lai_mean(x,y,t).eq.0.0d0) then
                        pft_max_lai_to_total_lai_ratio(x,y,t) = 0.0d0

                    else if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).lt.0.0d0.or.total_lai_mean(x,y,t).lt.0.0d0) then
                        stop "PFT max LAI / Total LAI"
                    endif
                else
                    pft_max_lai_to_total_lai_ratio(x,y,t) = fillvalue_dble
                end if

                ! PFT max LAI / Total tree LAI
                if (pft_max_lai_pftnum(x,y,t).gt.0.and.pft_max_lai_pftnum(x,y,t).lt.(ntree+1) &
                    .and.pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).gt.fillvalue_dble.and.total_tree_lai_mean(x,y,t).gt.fillvalue_dble) then
                    if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).gt.0.0d0.and.total_tree_lai_mean(x,y,t).gt.0.0d0) then
                        pft_max_lai_to_total_tree_lai_ratio(x,y,t) = pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t) / total_tree_lai_mean(x,y,t)

                    else if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).eq.0.0d0.or.total_tree_lai_mean(x,y,t).eq.0.0d0) then
                        pft_max_lai_to_total_tree_lai_ratio(x,y,t) = 0.0d0

                    else if (pft_lai_mean(x,y,pft_max_lai_pftnum(x,y,t),t).lt.0.0d0.or.total_tree_lai_mean(x,y,t).lt.0.0d0) then
                        stop "PFT max LAI / Total tree LAI"
                    endif
                else
                    pft_max_lai_to_total_tree_lai_ratio(x,y,t) = fillvalue_dble
                end if

                do p = 1,nphen_trees ! Phenology class LAI / Total LAI
                    if (phen_lai_mean(x,y,p,t).gt.fillvalue_dble.and.total_lai_mean(x,y,t).gt.fillvalue_dble) then
                        if (phen_lai_mean(x,y,p,t).gt.0.0d0.and.total_lai_mean(x,y,t).gt.0.0d0) then
                            phen_lai_to_total_lai_ratio(x,y,p,t) = phen_lai_mean(x,y,p,t) / total_lai_mean(x,y,t)
                        else if (phen_lai_annual(x,y,p,t).eq.0.0d0.or.total_lai_mean(x,y,t).eq.0.0d0) then
                            phen_lai_to_total_lai_ratio(x,y,p,t) = 0.0d0
                        else if (phen_lai_annual(x,y,p,t).lt.0.0d0.or.total_lai_mean(x,y,t).lt.0.0d0) then
                            stop "Phenology class LAI / Total LAI"
                        endif
                    else
                        phen_lai_to_total_lai_ratio(x,y,p,t) = fillvalue_dble
                    end if
                end do

                do p = 1,nphen_trees ! Phenology class LAI / Total tree LAI
                    if (phen_lai_mean(x,y,p,t).gt.fillvalue_dble.and.total_tree_lai_mean(x,y,t).gt.fillvalue_dble) then
                        if (phen_lai_mean(x,y,p,t).gt.0.0d0.and.total_tree_lai_mean(x,y,t).gt.0.0d0) then
                            phen_lai_to_total_tree_lai_ratio(x,y,p,t) = phen_lai_mean(x,y,p,t) / total_tree_lai_mean(x,y,t)
                        else if (phen_lai_annual(x,y,p,t).eq.0.0d0.or.total_tree_lai_mean(x,y,t).eq.0.0d0) then
                            phen_lai_to_total_tree_lai_ratio(x,y,p,t) = 0.0d0
                        else if (phen_lai_annual(x,y,p,t).lt.0.0d0.or.total_tree_lai_mean(x,y,t).lt.0.0d0) then
                            stop "Phenology class LAI / Total tree LAI"
                        endif
                    else
                        phen_lai_to_total_tree_lai_ratio(x,y,p,t)=fillvalue_dble
                    end if
                end do

                ! Phenology class maximum LAI / Total LAI
                if (phen_max_lai_phennum(x,y,t).gt.0.and.phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).gt.fillvalue_dble &
                    .and.total_lai_mean(x,y,t).gt.fillvalue_dble) then

                    if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).gt.0.0d0.and.total_lai_mean(x,y,t).gt.0.0d0) then
                        phen_max_lai_to_total_lai_ratio(x,y,t) = phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t) / total_lai_mean(x,y,t)
                    else if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).eq.0.0d0.or.total_lai_mean(x,y,t).eq.0.0d0) then
                        phen_max_lai_to_total_lai_ratio(x,y,t) = 0.0d0
                    else if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).lt.0.0d0.or.total_lai_mean(x,y,t).lt.0.0d0) then
                        stop "Phenology class maximum LAI / Total LAI"
                    endif
                else
                    phen_max_lai_to_total_lai_ratio(x,y,t) = fillvalue_dble
                end if

                ! Phenology class maximum LAI / Total tree LAI
                if (phen_max_lai_phennum(x,y,t).gt.0.and.phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).gt.fillvalue_dble &
                    .and.total_tree_lai_mean(x,y,t).gt.fillvalue_dble) then
                    if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).gt.0.0d0.and.total_tree_lai_mean(x,y,t).gt.0.0d0) then
                        phen_max_lai_to_total_tree_lai_ratio(x,y,t) = phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t) / total_tree_lai_mean(x,y,t)
                    else if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).eq.0.0d0.or.total_tree_lai_mean(x,y,t).eq.0.0d0) then
                        phen_max_lai_to_total_tree_lai_ratio(x,y,t) = 0.0d0
                    else if (phen_lai_mean(x,y,phen_max_lai_phennum(x,y,t),t).lt.0.0d0.or.total_tree_lai_mean(x,y,t).lt.0.0d0) then
                        stop "Phenology class maximum LAI / Total tree LAI"
                    endif
                else
                    phen_max_lai_to_total_tree_lai_ratio(x,y,t) = fillvalue_dble
                end if

                ! If all PFTs have mean LAIs = 0 then assign a 0 max LAI PFT and phenology class number.
                if (pft_max_lai(x,y,t).eq.0.0d0) pft_max_lai_pftnum(x,y,t) = 0
                if (phen_max_lai(x,y,t).eq.0.0d0) phen_max_lai_phennum(x,y,t) = 0
            end do ! Time
        end do ! Longitude
    end do ! Latitude

    if (modelnum.eq.4) deallocate (landcover_frac)

    ! Create netCDF output file
    call netcdf_create(outfile,nf90_netcdf4,ncid)
    call netcdf_defdim(ncid,"lon",nclon,londimid)
    call netcdf_defdim(ncid,"lat",nclat,latdimid)
    call netcdf_defdim(ncid,"PFT",nveg,pftdimid)
    call netcdf_defdim(ncid,"phenology_class",nphen_trees,phentreesdimid)
    call netcdf_defdim(ncid,"time",ntime-29,timedimid)

    call netcdf_defcoord(ncid,"lon",nf90_double,"longitude","degrees_east","X",londimid,lonid)
    call netcdf_put_att_text(ncid,lonid,"standard_name","longitude")
    call netcdf_put_att_text(ncid,lonid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,lonid,"_CoordinateAxisType","Lon")

    call netcdf_defcoord(ncid,"lat",nf90_double,"latitude","degrees_north","Y",latdimid,latid)
    call netcdf_put_att_text(ncid,latid,"standard_name","latitude")
    call netcdf_put_att_text(ncid,latid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,latid,"_CoordinateAxisType","Lat")

    call netcdf_deftime(ncid,"time",nf90_double,"time","time","year","T","noleap",timedimid,timeid)

    call netcdf_defvar_1d(ncid,"PFT",nf90_int,"plant functional types (PFTs)","1",pftdimid,pft_id)
    call netcdf_put_att_text(ncid,pft_id,"comment",pft_names)

    call netcdf_defvar_1d(ncid,"phenology_class",nf90_int,"phenology class","1",phentreesdimid,phentrees_id)
    call netcdf_put_att_text(ncid,phentrees_id,"comment",tree_phen_names)

    ! Define netCDF output file variables
    if (modelnum.gt.-1) then
! Total LAI
        chunk3d = (/nclon,nclat,1/)
        call netcdf_defvar_3d_compress(ncid,"lai_total",nf90_double, &
        "total LAI (30-year mean)","m2 m-2",londimid,latdimid,timedimid,lai_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,lai_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,lai_id,"comment","sum of the LAI for tree, shrub, grass, and tundra PFTs (if simulated by the model)")

! Total tree PFT LAI
        chunk3d = (/nclon,nclat,1/)
        call netcdf_defvar_3d_compress(ncid,"tree_PFT_lai_total",nf90_double, &
        "total tree PFT LAI (30-year mean)","m2 m-2",londimid,latdimid,timedimid,tree_lai_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,tree_lai_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,tree_lai_id,"comment","sum of the LAI for the tree PFTs simulated by the model")

! PFT 30-yr mean LAI
        chunk4d = (/nclon,nclat,nveg,1/)
        call netcdf_defvar_4d_compress(ncid,"pft_lai_mean",nf90_double, &
        "PFT LAI (30-year mean)","m2 m-2",londimid,latdimid,pftdimid,timedimid,pft_lai_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_lai_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,pft_lai_id,"comment","means calculated for all PFTs (tree and non-tree)")

! Phenology 30-yr mean LAI
        chunk4d = (/nclon,nclat,nphen_trees,1/)
        call netcdf_defvar_4d_compress(ncid,"phen_lai_mean",nf90_double, &
        "phenology class LAI (30-year mean)","m2 m-2",londimid,latdimid,phentreesdimid,timedimid,phen_lai_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_lai_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_lai_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_lai_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_lai_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! PFT max LAI
        chunk3d = (/nclon,nclat,1/)
        call netcdf_defvar_3d_compress(ncid,"pft_max_lai",nf90_double, &
        "PFT maximum LAI (30-year mean)","m2 m-2",londimid,latdimid,timedimid,pft_max_lai_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_max_lai_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,pft_max_lai_id,"comment","calculated for tree PFTs")

! PFT max LAI PFT number
        call netcdf_defvar_3d_compress(ncid,"pft_max_lai_pft_number",nf90_int, &
        "PFT number of the tree PFT with the maximum LAI (30-year mean)","1",londimid,latdimid,timedimid,pft_max_lai_pftnum_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_int(ncid,pft_max_lai_pftnum_id,"_FillValue",fillvalue_int)
        call netcdf_put_att_text(ncid,pft_max_lai_pftnum_id,"comment","0 = the maximum LAI for all &
        tree PFTs is 0.0 m2 m-2")

! PFT max LAI / total LAI
        call netcdf_defvar_3d_compress(ncid,"pft_max_lai_divided_by_total_lai",nf90_double, &
        "PFT maximum LAI (30-year mean) divided by total LAI (30-year mean)","m2 m-2 / m2 m-2",londimid,latdimid,timedimid,pft_max_lai_ratio_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_max_lai_ratio_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,pft_max_lai_ratio_id,"comment","PFT maximum LAI calculated for tree PFTs")

! PFT max LAI / total tree LAI
        call netcdf_defvar_3d_compress(ncid,"pft_max_lai_divided_by_total_tree_lai",nf90_double, &
        "PFT maximum LAI (30-year mean) divided by total tree LAI (30-year mean)","m2 m-2 / m2 m-2",londimid,latdimid,timedimid,pft_max_tree_lai_ratio_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_max_tree_lai_ratio_id,"_FillValue",fillvalue_dble)
        call netcdf_put_att_text(ncid,pft_max_tree_lai_ratio_id,"comment","PFT maximum LAI calculated for tree PFTs")

! Phenology max LAI
        call netcdf_defvar_3d_compress(ncid,"phen_max_lai",nf90_double, &
        "phenology class maximum LAI (30-year mean)","m2 m-2",londimid,latdimid,timedimid,phen_max_lai_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_max_lai_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_max_lai_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_max_lai_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_max_lai_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! Phenology max LAI phenology class number
        call netcdf_defvar_3d_compress(ncid,"phen_max_lai_phen_number",nf90_int, &
        "phenology class number of the phenology class with the maximum LAI (30-year mean)","1",londimid,latdimid,timedimid,phen_max_lai_phennum_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_int(ncid,phen_max_lai_phennum_id,"_FillValue",fillvalue_int)
        call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment","0 = the maximum LAI for all &
        tree phenology classes is 0.0 m2 m-2")
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! Phenology max LAI / total LAI
        call netcdf_defvar_3d_compress(ncid,"phen_max_lai_divided_by_total_lai",nf90_double, &
        "phenology class maximum LAI (30-year mean) divided by total LAI (30-year mean)","m2 m-2 / m2 m-2",londimid,latdimid,timedimid,phen_max_lai_ratio_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_max_lai_ratio_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_max_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_max_lai_ratio_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_max_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! Phenology max LAI / total tree LAI
        call netcdf_defvar_3d_compress(ncid,"phen_max_lai_divided_by_total_tree_lai",nf90_double, &
        "phenology class maximum LAI (30-year mean) divided by total tree_LAI (30-year mean)","m2 m-2 / m2 m-2",londimid,latdimid,timedimid,phen_max_tree_lai_ratio_id &
        ,chunk3d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_max_tree_lai_ratio_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_max_tree_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_max_tree_lai_ratio_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_max_tree_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! PFT LAI / total LAI
        chunk4d = (/nclon,nclat,nveg,1/)
        call netcdf_defvar_4d_compress(ncid,"pft_lai_divided_by_total_lai",nf90_double, &
        "PFT LAI (30-year mean) divided by total LAI (30-year mean)","m2 m-2 / m2 m-2", &
        londimid,latdimid,pftdimid,timedimid,pft_lai_ratio_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_lai_ratio_id,"_FillValue",fillvalue_dble)

! Phenology LAI / total LAI
        chunk4d = (/nclon,nclat,nphen_trees,1/)
        call netcdf_defvar_4d_compress(ncid,"phen_lai_divided_by_total_lai",nf90_double, &
        "phenology class LAI (30-year mean) divided by total LAI (30-year mean)","m2 m-2 / m2 m-2", &
        londimid,latdimid,phentreesdimid,timedimid,phen_lai_ratio_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_lai_ratio_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_lai_ratio_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

! PFT LAI / total tree LAI
        chunk4d = (/nclon,nclat,nveg,1/)
        call netcdf_defvar_4d_compress(ncid,"pft_lai_divided_by_total_tree_lai",nf90_double, &
        "PFT LAI (30-year mean) divided by total tree LAI (30-year mean)","m2 m-2 / m2 m-2", &
        londimid,latdimid,pftdimid,timedimid,pft_tree_lai_ratio_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,pft_tree_lai_ratio_id,"_FillValue",fillvalue_dble)

! Phenology LAI / total tree LAI
        chunk4d = (/nclon,nclat,nphen_trees,1/)
        call netcdf_defvar_4d_compress(ncid,"phen_lai_divided_by_total_tree_lai",nf90_double, &
        "phenology class LAI (30-year mean) divided by total tree LAI (30-year mean)","m2 m-2 / m2 m-2", &
        londimid,latdimid,phentreesdimid,timedimid,phen_tree_lai_ratio_id &
        ,chunk4d,dflevel,shuffle_filter)
        call netcdf_put_att_double(ncid,phen_tree_lai_ratio_id,"_FillValue",fillvalue_dble)
        if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_tree_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
        if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_tree_lai_ratio_id,"comment","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
        if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_tree_lai_ratio_id,"comment","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

    endif

    ! Add global attributes
    if (modelnum.ne.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","30-minute")
    if (modelnum.eq.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","1.25 degrees latitude x 1.875 degrees longitude")
    call netcdf_put_att_text(ncid,nf90_global,"projection","geographic")
    call current_time(time_stamp)
    history = adjustl(trim(time_stamp))//" -- File created by S.L. Shafer"
    call netcdf_put_att_text(ncid,nf90_global,"history",history)
    call netcdf_put_att_text(ncid,nf90_global,"source_file",source)
    call netcdf_put_att_text(ncid,nf90_global,"references","LAI data: Pugh et al. (in prep.)")
    call netcdf_put_att_text(ncid,nf90_global,"contact","Sarah L. Shafer (sshafer@usgs.gov)")
    call netcdf_put_att_text(ncid,nf90_global,"institution","U.S. Geological Survey")
    call netcdf_put_att_text(ncid,nf90_global,"comment","means calculated using data for the previous 30 years")
    call netcdf_put_att_text(ncid,nf90_global,"comment2","analyses done for grid cells identified as forest for 1985-2014 by the forest_30yr_any_10_years variable in the forest mask file")

    call netcdf_put_att_text(ncid,nf90_global,"Conventions","CF-1.0")
    call netcdf_put_att_text(ncid,nf90_global,"disclaimer","Although these data have been processed successfully on a computer system at &
the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the &
display or utility of the data on any system, or for general or scientific purposes, nor &
shall the act of distribution constitute any such warranty. The USGS shall not be held liable &
for improper or incorrect use of the data described and/or contained herein. Any use of trade, &
firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.")

    ! Leave netCDF output file define mode
    call netcdf_enddef(outfile,ncid)

    ! Put data into the netCDF output file
    call netcdf_put_1d_double(ncid,lonid,nclon,clon)
    call netcdf_put_1d_double(ncid,latid,nclat,clat)
    call netcdf_put_1d_double(ncid,timeid,ntime-29,time(30:114))
    call netcdf_put_1d_int(ncid,pft_id,nveg,pft)
    call netcdf_put_1d_int(ncid,phentrees_id,nphen_trees,phen_trees)

    call netcdf_put_3d_double(ncid,lai_id,nclon,nclat,ntime-29,total_lai_mean(:,:,30:114))
    call netcdf_put_3d_double(ncid,tree_lai_id,nclon,nclat,ntime-29,total_tree_lai_mean(:,:,30:114))

    call netcdf_put_4d_double(ncid,pft_lai_id,nclon,nclat,nveg,ntime-29,pft_lai_mean(:,:,:,30:114))
    call netcdf_put_4d_double(ncid,phen_lai_id,nclon,nclat,nphen_trees,ntime-29,phen_lai_mean(:,:,:,30:114))

    call netcdf_put_3d_double(ncid,pft_max_lai_id,nclon,nclat,ntime-29,pft_max_lai(:,:,30:114))
    call netcdf_put_3d_int(ncid,pft_max_lai_pftnum_id,nclon,nclat,ntime-29,pft_max_lai_pftnum(:,:,30:114))
    call netcdf_put_3d_double(ncid,phen_max_lai_id,nclon,nclat,ntime-29,phen_max_lai(:,:,30:114))
    call netcdf_put_3d_int(ncid,phen_max_lai_phennum_id,nclon,nclat,ntime-29,phen_max_lai_phennum(:,:,30:114))

    call netcdf_put_3d_double(ncid,pft_max_lai_ratio_id,nclon,nclat,ntime-29,pft_max_lai_to_total_lai_ratio(:,:,30:114))
    call netcdf_put_3d_double(ncid,phen_max_lai_ratio_id,nclon,nclat,ntime-29,phen_max_lai_to_total_lai_ratio(:,:,30:114))

    call netcdf_put_4d_double(ncid,pft_lai_ratio_id,nclon,nclat,nveg,ntime-29,pft_lai_to_total_lai_ratio(:,:,:,30:114))
    call netcdf_put_4d_double(ncid,phen_lai_ratio_id,nclon,nclat,nphen_trees,ntime-29,phen_lai_to_total_lai_ratio(:,:,:,30:114))

    call netcdf_put_3d_double(ncid,pft_max_tree_lai_ratio_id,nclon,nclat,ntime-29,pft_max_lai_to_total_tree_lai_ratio(:,:,30:114))
    call netcdf_put_3d_double(ncid,phen_max_tree_lai_ratio_id,nclon,nclat,ntime-29,phen_max_lai_to_total_tree_lai_ratio(:,:,30:114))

    call netcdf_put_4d_double(ncid,pft_tree_lai_ratio_id,nclon,nclat,nveg,ntime-29,pft_lai_to_total_tree_lai_ratio(:,:,:,30:114))
    call netcdf_put_4d_double(ncid,phen_tree_lai_ratio_id,nclon,nclat,nphen_trees,ntime-29,phen_lai_to_total_tree_lai_ratio(:,:,:,30:114))

    ! Close the netCDF output file
    call netcdf_close(outfile,ncid)

    ! Create netCDF output file with single variable
    call netcdf_create(outfile_mask,nf90_netcdf4,ncid)
    call netcdf_defdim(ncid,"lon",nclon,londimid)
    call netcdf_defdim(ncid,"lat",nclat,latdimid)
    call netcdf_defdim(ncid,"phen_class",nphen_trees,phentreesdimid)

    call netcdf_defcoord(ncid,"lon",nf90_double,"longitude","degrees_east","X",londimid,lonid)
    call netcdf_put_att_text(ncid,lonid,"standard_name","longitude")
    call netcdf_put_att_text(ncid,lonid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,lonid,"_CoordinateAxisType","Lon")

    call netcdf_defcoord(ncid,"lat",nf90_double,"latitude","degrees_north","Y",latdimid,latid)
    call netcdf_put_att_text(ncid,latid,"standard_name","latitude")
    call netcdf_put_att_text(ncid,latid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,latid,"_CoordinateAxisType","Lat")

    call netcdf_defvar_1d(ncid,"phenology_class",nf90_int,"phenology class","1",phentreesdimid,phentrees_id)
    call netcdf_put_att_text(ncid,phentrees_id,"comment",tree_phen_names)

    ! Phenology max LAI phenology class number
    call netcdf_defvar_2d(ncid,"phen_max_lai_phen_number",nf90_int, &
        "phenology class number of the phenology class with the maximum LAI (1985-2014 30-year mean)","1",londimid,latdimid,phen_max_lai_phennum_id)
    call netcdf_put_att_int(ncid,phen_max_lai_phennum_id,"_FillValue",fillvalue_int)
    call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment","0 = the maximum LAI for all &
        tree phenology classes is 0.0 m2 m-2")
    if (modelnum.eq.2) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 2 (needleleaf deciduous) is not simulated by LPJmL")
    if (modelnum.eq.3) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 3 (boreal broadleaved deciduous) is not simulated by LPJ-GUESS")
    if (modelnum.eq.4) call netcdf_put_att_text(ncid,phen_max_lai_phennum_id,"comment2","phenology class 2 (needleleaf deciduous) is not simulated by LPJ-wsl")

    ! Add global attributes
    if (modelnum.ne.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","30-minute")
    if (modelnum.eq.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","1.25 degrees latitude x 1.875 degrees longitude")
    call netcdf_put_att_text(ncid,nf90_global,"projection","geographic")
    call current_time(time_stamp)
    history = adjustl(trim(time_stamp))//" -- File created by S.L. Shafer"
    call netcdf_put_att_text(ncid,nf90_global,"history",history)
    call netcdf_put_att_text(ncid,nf90_global,"source_file",source)
    call netcdf_put_att_text(ncid,nf90_global,"references","LAI data: Pugh et al. (in prep.)")
    call netcdf_put_att_text(ncid,nf90_global,"contact","Sarah L. Shafer (sshafer@usgs.gov)")
    call netcdf_put_att_text(ncid,nf90_global,"institution","U.S. Geological Survey")
    call netcdf_put_att_text(ncid,nf90_global,"comment","mean calculated using data for 1985-2014")
    call netcdf_put_att_text(ncid,nf90_global,"comment2","analyses done for grid cells identified as forest for 1985-2014 by the forest_30yr_any_10_years variable in the forest mask file")

    call netcdf_put_att_text(ncid,nf90_global,"Conventions","CF-1.0")
    call netcdf_put_att_text(ncid,nf90_global,"disclaimer","Although these data have been processed successfully on a computer system at &
the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the &
display or utility of the data on any system, or for general or scientific purposes, nor &
shall the act of distribution constitute any such warranty. The USGS shall not be held liable &
for improper or incorrect use of the data described and/or contained herein. Any use of trade, &
firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.")

    ! Leave netCDF output file define mode
    call netcdf_enddef(outfile_mask,ncid)

    ! Put data into the netCDF output file
    call netcdf_put_1d_double(ncid,lonid,nclon,clon)
    call netcdf_put_1d_double(ncid,latid,nclat,clat)
    call netcdf_put_1d_int(ncid,phentrees_id,nphen_trees,phen_trees)

    call netcdf_put_2d_int(ncid,phen_max_lai_phennum_id,nclon,nclat,phen_max_lai_phennum(:,:,114))

    ! Close the netCDF output file
    call netcdf_close(outfile_mask,ncid)

    ! Deallocate arrays
    deallocate (clat,clon,time)
    if (modelnum.lt.6) deallocate (time_months)
    deallocate (pft_phenclass_crosswalk,pft,mask)
    deallocate (lai_annual,total_lai_annual,total_tree_lai,phen_lai_annual)
    deallocate (tree_max_lai,tree_max_lai_pft)
    deallocate (pft_lai_mean,phen_lai_mean,total_lai_mean,total_tree_lai_mean)
    deallocate (pft_max_lai,pft_max_lai_pftnum,phen_max_lai,phen_max_lai_phennum)
    deallocate (pft_lai_to_total_lai_ratio,pft_max_lai_to_total_lai_ratio)
    deallocate (phen_lai_to_total_lai_ratio,phen_max_lai_to_total_lai_ratio)
    deallocate (pft_lai_to_total_tree_lai_ratio,pft_max_lai_to_total_tree_lai_ratio)
    deallocate (phen_lai_to_total_tree_lai_ratio,phen_max_lai_to_total_tree_lai_ratio)

end do ! Model loop
end
