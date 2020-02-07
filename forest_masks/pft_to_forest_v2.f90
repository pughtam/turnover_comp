program pft_to_forest_v2
!------------------------------------------------------------------------------------------
! This program creates a forest mask from simulated monthly or annual LAI. See Pugh et al.
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

integer, parameter      :: nd = 4, npft = 10

real(8)                 :: fillvalue
real(8)                 :: lai_max,fillvalue_dble

real(8), allocatable    :: clon(:),clat(:)
real(8), allocatable    :: lai_month_model(:,:,:,:),lai_month(:,:,:,:),lai_annual(:,:,:,:)
real(8), allocatable    :: fpc(:,:,:,:)
real(8), allocatable    :: total_tree_lai(:,:,:),tree_max_lai(:,:,:)
real(8), allocatable    :: time_months(:),time(:)

real(4)                 :: filltemp

real(4), allocatable    :: clon_float(:),clat_float(:),time_months_float(:)
real(4), allocatable    :: lai_month_float(:,:,:,:),fpc_float(:,:,:,:),lai_month_temp(:,:,:,:)
real(4), allocatable    :: lai_annual_float(:,:,:,:)

integer, allocatable    :: forest(:,:,:),mask(:,:),time_months_int(:)
integer, allocatable    :: forest30_10(:,:),tree_max_lai_pft(:,:,:)

integer                 :: modelnum,nveg,ntime,nmonths,year,fill_check,begin30,end30
integer                 :: forest_id
integer                 :: nclon,nclat
integer                 :: forest_count
integer                 :: grass30
integer                 :: forest30_10_id
 
integer                 :: i,m,x,y,t,z,fillvalue_int

character(4)            :: run
character(9)            :: text30
character(19)           :: time_stamp
character(64)           :: veg_name,lai_name,latname,lonname

character(1024)         :: ncpath,ncoutpath,infile,infile2,outfile,fpcfile

ncpath = "\projects\mortality\data\"
ncoutpath = "\projects\mortality\"

latname = "lat"
lonname = "lon"
fillvalue = -1.0E32
fillvalue_dble = -1.0D32
fillvalue_int = -9999

! Models
! modelnum = 0  ! SEIB
! modelnum = 1  ! JULES
! modelnum = 2  ! LPJmL
! modelnum = 3  ! LPJ-GUESS
! modelnum = 4  ! LPJ-wsl
! modelnum = 5  ! CABLE-POP
! modelnum = 6  ! ORCHIDEE

! Choose data to input
!run = "ncep"
run = "ipsl"

if (run.eq."ncep") begin30 = 85  ! 1985
if (run.eq."ncep") end30 = 114   ! 2014
if (run.eq."ncep") text30 = "1985-2014"

if (run.eq."ipsl") begin30 = 170 ! 2070
if (run.eq."ipsl") end30 = 199   ! 2099
if (run.eq."ipsl") text30 = "2070-2099"

do modelnum = 3,3

! Assign model details
! SEIB
    if (modelnum.eq.0) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"seib_cruncep_lai_month_1901_2014.nc4"
            outfile = trim(ncoutpath)//"seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "seib_cruncep_lai_month_1901_2014.nc4"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\seib_ipsl_lai_month_1901_2099.nc4"
            outfile = trim(ncoutpath)//"seib_ipsl_lai_annual_1901_2099_forest_mask_v4.nc"
            source = "seib_ipsl_lai_annual_1901_2099_forest_mask_v4.nc"
        endif

        veg_name = "vegtype"
        lai_name = "lai"
        latname = "latitude"
        lonname = "longitude"
        source = "seib_cruncep_lai_month_1901_2014.nc4"
        nclat = 360
        nclon = 720

! JULES
    else if (modelnum.eq.1) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"JULESC2_CRUNCEP_lai_Monthly_1901_2014.nc"
            outfile = trim(ncoutpath)//"JULESC2_CRUNCEP_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "JULESC2_CRUNCEP_lai_Monthly_1901_2014.nc"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\JULESC2_IPSL_lai_Monthly_1901_2005.nc4"
            infile2 = trim(ncpath)//"ipsl\JULESC2_IPSL_lai_Monthly_2005_2099.nc4"
            outfile = trim(ncoutpath)//"JULESC2_IPSL_lai_annual_1901_2099_forest_mask_v4.nc"
            source = "JULESC2_IPSL_lai_annual_1901_2099_forest_mask_v4.nc"
        end if

        veg_name = "pft"
        lai_name = "lai"
        latname = "lat"
        lonname = "lon"
        nclat = 112
        nclon = 192

! LPJmL
    else if (modelnum.eq.2) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"lpjml_cruncep_lai_annual_1901_2014.nc"
            outfile = trim(ncoutpath)//"lpjml_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "lpjml_cruncep_lai_annual_1901_2014.nc"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\lpjml_ipsl-cm5a_lr_lai_annual_1901_2099.nc"
            outfile = trim(ncoutpath)//"lpjml_ipsl-cm5a_lr_lai_annual_1901_2099_forest_mask_v4.nc"
            source = "lpjml_ipsl-cm5a_lr_lai_annual_1901_2099.nc"
        end if

        veg_name = "npft"
        lai_name = "lai"
        latname = "latitude"
        lonname = "longitude"
        nclat = 279
        nclon = 720

! LPJ-GUESS
    else if (modelnum.eq.3) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"lpj-guess_cruncep_lai_monthly_1901_2014.nc4"
            outfile = trim(ncoutpath)//"lpj-guess_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "lpj-guess_cruncep_lai_monthly_1901_2014.nc4"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\lpj-guess_ipsl_lai_monthly_1901_2014.nc4"
            outfile = trim(ncoutpath)//"lpj-guess_ipsl_lai_annual_1901_2100_forest_mask_v4.nc"
            source = "lpj-guess_ipsl_lai_monthly_1901_2014.nc4"
        end if

        veg_name = "vegtype"
        lai_name = "lai"
        latname = "latitude"
        lonname = "longitude"
        nclat = 360
        nclon = 720

! LPJ-wsl
    else if (modelnum.eq.4) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"LPJ_pft_mlai_BONE.nc"
            fpcfile = trim(ncpath)//"LPJ_landCoverFrac.nc"
            outfile = trim(ncoutpath)//"lpj-wsl_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "LPJ_pft_mlai_[pft name].nc, LPJ_landCoverFrac.nc"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\LPJ_pft_mlai_BONE.nc"
            fpcfile = trim(ncpath)//"ipsl\LPJ_landCoverFrac.nc"
            outfile = trim(ncoutpath)//"lpj-wsl_ipsl_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "LPJ_pft_mlai_[pft name].nc, LPJ_landCoverFrac.nc"
        end if

        lai_name = "lai"
        latname = "lat"
        lonname = "lon"
        nclat = 360
        nclon = 720

! CABLE-POP
    else if (modelnum.eq.5) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"CABLE-POP_cruncep_lai_month_1901_2015.nc4"
            outfile = trim(ncoutpath)//"CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc"
            source = "CABLE-POP_cruncep_lai_month_1901_2015.nc"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl/CABLE-POP_ipsl-cm5a-lr_lai_month_1901_2099.nc4"
            outfile = trim(ncoutpath)//"CABLE-POP_ipsl-cm5a-lr_lai_annual_1901_2099_forest_mask_v4.nc"
            source = "CABLE-POP_ipsl-cm5a-lr_lai_month_1901_2099.nc"
        end if

        lai_name = "lai"
        latname = "y"
        lonname = "x"
        nclat = 360
        nclon = 720

! ORCHIDEE
    else if (modelnum.eq.6) then
        if (run.eq."ncep") then
            infile = trim(ncpath)//"ORCHIDEEr3085_cruncep_lai_13pft_year_1901_2014.nc"
            fpcfile = trim(ncpath)//"ORCHIDEEr3085_cruncep_veget_max_13pft_year_1901_2014.nc"
            outfile = trim(ncoutpath)//"ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"
            source = "ORCHIDEEr3085_cruncep_lai_13pft_year_1901_2014.nc, ORCHIDEEr3085_cruncep_veget_max_13pft_year_1901_2014.nc"
        else if (run.eq."ipsl") then
            infile = trim(ncpath)//"ipsl\ORCHIDEEr3085_ipslcm5alr_lai_13pft_year_1901_2099.nc"
            fpcfile = trim(ncpath)//"ipsl\ORCHIDEEr3085_ipslcm5alr_veget_max_13pft_year_1901_2099.nc"
            outfile = trim(ncoutpath)//"ORCHIDEE_ipslcm5alr_lai_annual_1901_2099_forest_mask_v4.nc"
            source = "ORCHIDEEr3085_ipslcm5alr_lai_13pft_year_1901_2099.nc, ORCHIDEEr3085_cruncep_veget_max_13pft_year_1901_2099.nc"
        end if

        lai_name = "lai"
        latname = "lat"
        lonname = "lon"
        nclat = 360
        nclon = 720
    end if

    allocate (clat(nclat),clon(nclon))
    if (modelnum.eq.1.or.modelnum.eq.2.or.modelnum.eq.5.or.modelnum.eq.6) allocate (clat_float(nclat),clon_float(nclon))

    write (*,*) "Reading input data..."

    ! Get latitudes and longitudes
    call netcdf_open(infile,0,ncid)
    if (modelnum.eq.0.or.modelnum.eq.3.or.modelnum.eq.4) then
        call netcdf_get_1d_double(ncid,nclon,trim(lonname),clon)
        call netcdf_get_1d_double(ncid,nclat,trim(latname),clat)
    else if (modelnum.eq.1.or.modelnum.eq.2.or.modelnum.eq.5.or.modelnum.eq.6) then
        call netcdf_get_1d_real(ncid,nclon,trim(lonname),clon_float)
        call netcdf_get_1d_real(ncid,nclat,trim(latname),clat_float)
        clat = dble(clat_float)
        clon = dble(clon_float)
        deallocate (clat_float,clon_float)
    end if

    ! Get time variable
    if (modelnum.lt.6) then
        call netcdf_get_dimlen(ncid,"time",nmonths)
        allocate (time_months(nmonths))
        if (modelnum.eq.1) then 
            allocate (time_months_float(nmonths))
            call netcdf_get_1d_real(ncid,nmonths,"time",time_months_float)
            time_months = dble(time_months_float)
            deallocate (time_months_float)
        else if (modelnum.eq.2.or.modelnum.eq.4) then
            allocate (time_months_int(nmonths))
            call netcdf_get_1d_int(ncid,nmonths,"time",time_months_int)
            time_months = dble(time_months_int)
            deallocate (time_months_int)
        else if (modelnum.eq.0.or.modelnum.eq.3.or.modelnum.eq.5) then
            call netcdf_get_1d_double(ncid,nmonths,"time",time_months)
        end if
    else if (modelnum.eq.6) then
        call netcdf_get_dimlen(ncid,"time_counter",ntime) ! Annual data
        allocate (time(ntime))
        call netcdf_get_1d_double(ncid,ntime,"time_counter",time)
    endif

    if (modelnum.lt.4) call netcdf_get_dimlen(ncid,veg_name,nveg)
    if (modelnum.eq.4) nveg = 9   ! LPJ-wsl
    if (modelnum.eq.5) nveg = 10  ! CABLE-POP
    if (modelnum.eq.6) nveg = 13  ! ORCHIDEE

    if (modelnum.eq.0) then ! SEIB
        allocate (lai_month_model(nclon,nclat,nmonths,nveg))
        allocate (lai_month(nclon,nclat,nveg,nmonths))
        call netcdf_get_4d_double(ncid,nclon,nclat,nmonths,nveg,trim(lai_name),lai_month_model)
        call netcdf_close(infile,ncid)

        ! Restructure array
        do x = 1,nclon
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
            ! Read data from BONE file, which is already open
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

        call netcdf_open(fpcfile,0,ncid) ! veget_max = Maximum vegetation fraction (LAI -> infinity)
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
        ntime = nmonths ! Annual data
        allocate (time(ntime))

    ! LPJ-GUESS
    else if (modelnum.eq.3) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! LPJ-wsl
    else if (modelnum.eq.4) then
        ntime = nmonths / 12
        allocate (time(ntime))

        allocate (fpc(nclon,nclat,ntime,nveg))

        call netcdf_open(fpcfile,0,ncid)
            call netcdf_get_4d_double(ncid,nclon,nclat,ntime,nveg,"landCoverFrac",fpc)
        call netcdf_close(fpcfile,ncid)

    ! CABLE-POP
    else if (modelnum.eq.5) then
        ntime = nmonths / 12
        allocate (time(ntime))

    ! ORCHIDEE
    else if (modelnum.eq.6) then
        ! ORCHIDEE data are annual, not monthly, so time array read directly from file (see above code).
        ! This else-if statement is intentionally blank.
    end if

    do t = 1,ntime
        time(t) = (dble(t) * 1.0d0) + 1900.0d0
    end do

    if (modelnum.ne.2.and.modelnum.ne.6) then
        allocate (lai_annual(nclon,nclat,nveg,ntime))
        lai_annual = 0.0d0
    end if

    ! Determine annual maximum monthly LAI for each vegetation type.
    if (modelnum.ne.2.and.modelnum.ne.6) then
        do y = 1,nclat
            if (mod(y,100).eq.0) write (*,*) "row = ",y
            do x = 1,nclon
                year = 0
                do t = 1,nmonths,12
                    year = year + 1
                    do i = 1,nveg
                        fill_check = 0
                        lai_max = fillvalue
                        do m = 0,11
                            if (lai_month(x,y,i,t+m).ge.0.00000000d0) fill_check = 1
                            if (lai_month(x,y,i,t+m).gt.lai_max) lai_max = lai_month(x,y,i,t+m)
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
                        ! LPJmL fill value = -99999.
                        if (lai_annual(x,y,i,t).lt.-99990.0d0) lai_annual(x,y,i,t) = fillvalue_dble
                    end do
                end do
            end do
        end do
    else if (modelnum.eq.6) then ! ORCHIDEE
        ! For ORCHIDEE no additional modification needed.
    endif

    allocate (mask(nclon,nclat))
    mask = 0

    allocate (total_tree_lai(nclon,nclat,ntime),tree_max_lai(nclon,nclat,ntime),tree_max_lai_pft(nclon,nclat,ntime))
    total_tree_lai = 0.0d0

    write (*,*) "Calculating total tree lai..."

    do y = 1,nclat
        do x = 1,nclon
            do t = 1,ntime
                tree_max_lai(x,y,t) = -9999.0d0
                tree_max_lai_pft(x,y,t) = -9

                if (modelnum.eq.0) then ! SEIB
                    do z = 1,12
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.1) then ! JULES
                    do z = 1,5
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.2) then ! LPJmL
                    do z = 1,7
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.3) then ! LPJ-GUESS
                    do z = 1,9
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.4) then ! LPJ-wsl
                ! LAI multiplied by FPC (B. Poulter, personal communication)
                    do z = 1,7
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + (lai_annual(x,y,z,t) * fpc(x,y,t,z))
                            if ((lai_annual(x,y,z,t)*fpc(x,y,t,z)).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t) * fpc(x,y,t,z)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.5) then ! CABLE-POP
                    do z = 1,4
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs

                else if (modelnum.eq.6) then ! ORCHIDEE
                    do z = 2,9
                        if (lai_annual(x,y,z,t).ne.fillvalue_dble) then
                            total_tree_lai(x,y,t) = total_tree_lai(x,y,t) + lai_annual(x,y,z,t)
                            if (lai_annual(x,y,z,t).gt.tree_max_lai(x,y,t)) then
                                tree_max_lai(x,y,t) = lai_annual(x,y,z,t)
                                tree_max_lai_pft(x,y,t) = z
                            endif
                        endif
                    end do ! Tree PFTs
                end if

                do i = 1,nveg
                    if (lai_annual(x,y,i,t).ne.fillvalue_dble) then
                        mask(x,y) = 1
                    endif
                end do
            end do
        end do
    end do
    if (modelnum.eq.4) deallocate (fpc)

    deallocate (lai_annual)
    allocate (forest(nclon,nclat,ntime))
    forest = fillvalue_int

    write (*,*) "Assigning forest..."
    do y = 1,nclat
        write (*,*) y,clat(y)
        do x = 1,nclon
            do t = 1,ntime

                if (mask(x,y).eq.1) then

                    if (modelnum.ne.-9) then
                        if (total_tree_lai(x,y,t).gt.2.50d0) forest(x,y,t) = 1
                        if (total_tree_lai(x,y,t).le.2.50d0) forest(x,y,t) = 2
                    endif

                    ! Boreal forest
                    if (modelnum.eq.0.and.total_tree_lai(x,y,t).gt.0.50d0) then !Boreal forest
                        if ((tree_max_lai_pft(x,y,t).gt.9.and.tree_max_lai_pft(x,y,t).lt.13)) forest(x,y,t) = 1

                    else if (modelnum.eq.1.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if (tree_max_lai_pft(x,y,t).gt.2.and.tree_max_lai_pft(x,y,t).lt.6) forest(x,y,t) = 1

                    else if (modelnum.eq.2.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if ((tree_max_lai_pft(x,y,t).gt.5.and.tree_max_lai_pft(x,y,t).lt.8)) forest(x,y,t) = 1

                    else if (modelnum.eq.3.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if (tree_max_lai_pft(x,y,t).lt.4) forest(x,y,t) = 1

                    else if (modelnum.eq.4.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if (tree_max_lai_pft(x,y,t).eq.6.or.tree_max_lai_pft(x,y,t).eq.7) forest(x,y,t) = 1

                    else if (modelnum.eq.5.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if (tree_max_lai_pft(x,y,t).lt.5.and.tree_max_lai_pft(x,y,t).ne.2) forest(x,y,t) = 1

                    else if (modelnum.eq.6.and.total_tree_lai(x,y,t).gt.0.50d0) then
                        if (tree_max_lai_pft(x,y,t).gt.6.and.tree_max_lai_pft(x,y,t).lt.10) forest(x,y,t) = 1
                    endif
                else
                    total_tree_lai(x,y,t) = fillvalue_dble
                end if ! Mask
            end do ! Time loop
        end do ! Longitude loop
    end do ! Latitude loop

    allocate (forest30_10(nclon,nclat))
    forest30_10 = fillvalue_int

    do x = 1,nclon
        do y = 1,nclat
            forest_count = 0
            grass30 = 0

            do t = begin30,end30
                if (forest(x,y,t).eq.1) forest_count = forest_count + 1
                if (forest(x,y,t).eq.2) grass30 = 2
            end do

            if (forest_count.ge.10) forest30_10(x,y) = 1
            if (forest_count.lt.10.and.grass30.eq.2) forest30_10(x,y) = 2
            if (forest_count.lt.10.and.grass30.eq.0.and.mask(x,y).eq.1) forest30_10(x,y) = 3

        end do
    end do

    ! Create netCDF output file
    call netcdf_create(outfile,nf90_netcdf4,ncid)
    call netcdf_defdim(ncid,"lon",nclon,londimid)
    call netcdf_defdim(ncid,"lat",nclat,latdimid)
    call netcdf_defdim(ncid,"time",ntime,timedimid)

    call netcdf_defcoord(ncid,"lon",nf90_double,"longitude","degrees_east","X",londimid,lonid)
    call netcdf_put_att_text(ncid,lonid,"standard_name","longitude")
    call netcdf_put_att_text(ncid,lonid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,lonid,"_CoordinateAxisType","Lon")

    call netcdf_defcoord(ncid,"lat",nf90_double,"latitude","degrees_north","Y",latdimid,latid)
    call netcdf_put_att_text(ncid,latid,"standard_name","latitude")
    call netcdf_put_att_text(ncid,latid,"coodinate_defines","gridcell_center")
    call netcdf_put_att_text(ncid,latid,"_CoordinateAxisType","Lat")

    call netcdf_deftime(ncid,"time",nf90_double,"time","time","year","T","noleap",timedimid,timeid)

    ! Define netCDF output file variables
    if (modelnum.gt.-1) then
        call netcdf_defvar_3d(ncid,"forest",nf90_int,"forest","1",londimid,latdimid,timedimid,forest_id)
        call netcdf_put_att_int(ncid,forest_id,"_FillValue",fillvalue_int)
        call netcdf_put_att_text(ncid,forest_id,"comment","1=forest, 2=non-forest, 3=unassigned")

        call netcdf_defvar_2d(ncid,"forest_30yr_any_10_years",nf90_int,"forest any 10 years, "//trim(text30),"1",londimid,latdimid,forest30_10_id)
        call netcdf_put_att_int(ncid,forest30_10_id,"_FillValue",fillvalue_int)
        call netcdf_put_att_text(ncid,forest30_10_id,"comment","1=forest, 2=non-forest, 3=unassigned")
        call netcdf_put_att_text(ncid,forest30_10_id,"comment2","1=forest simulated for at least 10 of the 30 years")
    endif

    ! Add global attributes to netCDF output file
    if (modelnum.ne.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","30-minute")
    if (modelnum.eq.1) call netcdf_put_att_text(ncid,nf90_global,"grid_spacing","30-minute")

    call netcdf_put_att_text(ncid,nf90_global,"projection","geographic")
    call current_time(time_stamp)
    history = adjustl(trim(time_stamp))//" -- File created by S. Shafer"
    call netcdf_put_att_text(ncid,nf90_global,"history",history)
    call netcdf_put_att_text(ncid,nf90_global,"source_file",source)
    call netcdf_put_att_text(ncid,nf90_global,"references","LAI data: Pugh et al. (in prep.)")
    call netcdf_put_att_text(ncid,nf90_global,"contact","Sarah L. Shafer (sshafer@usgs.gov)")
    call netcdf_put_att_text(ncid,nf90_global,"institution","U.S. Geological Survey")

    if (modelnum.ne.-9.and.modelnum.ne.1.and.modelnum.ne.5) comment = "Forest is assigned if 1) annual total tree LAI > 2.5 or 2) &
    annual total tree LAI > 0.5 and the PFT with the maximum LAI is a boreal tree PFT."

    if (modelnum.eq.1.or.modelnum.eq.5) comment = "Forest is assigned if 1) annual total tree LAI > 2.5 or 2) annual total tree &
    LAI > 0.5 and the PFT with the maximum LAI is a needleleaf evergreen, needleleaf deciduous, or broadleaf deciduous tree PFT &
    (i.e., a PFT that could be boreal)."

    if (modelnum.eq.-9) comment = "Forest is assigned if annual total tree LAI is greater than the LAI &
    threshold specified in the variable comment field."
    call netcdf_put_att_text(ncid,nf90_global,"comment",comment)

    if (modelnum.eq.1) call netcdf_put_att_text(ncid, nf90_global,"comment2","JULES LAI data originally simulated on a 1.25 degrees latitude x &
    1.875 degrees longitude grid were regridded to a 0.5-degree grid before assigning forest.")

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
    call netcdf_put_1d_double(ncid,timeid,ntime,time)

    if (modelnum.ne.-9) then
        call netcdf_put_3d_int(ncid,forest_id,nclon,nclat,ntime,forest)
        call netcdf_put_2d_int(ncid,forest30_10_id,nclon,nclat,forest30_10)
    endif

    ! Close the netCDF output file
    call netcdf_close(outfile,ncid)

    ! Deallocate arrays
    deallocate (clat,clon,time)
    if (modelnum.lt.6) deallocate (time_months)
    deallocate (mask)
    deallocate (forest)
    deallocate (total_tree_lai)
    deallocate (tree_max_lai,tree_max_lai_pft)
    deallocate (forest30_10)
 
end do ! Model loop
end
