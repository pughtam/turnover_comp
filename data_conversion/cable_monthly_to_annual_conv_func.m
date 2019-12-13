%Script to read in CABLE-POP monthly by pft data and convert it to annual for
%manageability
%
%T. Pugh
%25.04.16

indir='/home/pughtam/data/turnover/CABLE/CRUNCEP/';
outdir='/home/pughtam/data/turnover/CABLE/CRUNCEP/';

file_in='CABLE-POP_cruncep_gpp_month_1901_2014.nc';
file_out='CABLE-POP_cruncep_gpp_annual_1901_2014.nc4';

varname='gpp';

ifpft=false; %Does the file have a PFT dimension?
npft=11; %Number of PFTs (if ifpft==true)

%---
time_in=ncread([indir,file_in],'time');
lon=ncread([indir,file_in],'x');
lat=ncread([indir,file_in],'y');
nlat=length(lat);
nlon=length(lon);
nmon=length(time_in);
nyear=nmon/12;

if ifpft
    data_out=NaN(nlon,nlat,nyear,npft);
    for yy=1:nyear
        mon_s=(yy*12)-11;
        data_in=ncread([indir,file_in],varname,[1 1 1 mon_s],[Inf Inf Inf 12]);
        data_out(:,:,yy,:)=squeeze(nansum(data_in,4)/12);
    end
else
    data_out=NaN(nlon,nlat,nyear);
    for yy=1:nyear
        mon_s=(yy*12)-11;
        data_in=ncread([indir,file_in],varname,[1 1 mon_s],[Inf Inf 12]);
        data_out(:,:,yy,:)=nansum(data_in,3)/12;
    end
end

%Write netcdf
ncid = netcdf.create([outdir,'/',file_out], 'NETCDF4');

dimid_lon=netcdf.defDim(ncid,'Longitude',nlon);
dimid_lat=netcdf.defDim(ncid,'Latitude',nlat);
dimid_time=netcdf.defDim(ncid,'Time',nyear);
varid_lon=netcdf.defVar(ncid,'Longitude','double',dimid_lon);
varid_lat=netcdf.defVar(ncid,'Latitude','double',dimid_lat);
varid_time=netcdf.defVar(ncid,'Time','double',dimid_time);
netcdf.putVar(ncid,varid_lon,lon)
netcdf.putVar(ncid,varid_lat,lat)

netcdf.putVar(ncid,varid_time,1:nyear)
netcdf.putAtt(ncid,varid_time,'Units','Calendar year')

if ifpft
    dimid_pft=netcdf.defDim(ncid,'vegtype',npft);
    varid_pft=netcdf.defVar(ncid,'vegtype','double',dimid_pft);
    netcdf.putVar(ncid,varid_pft,1:npft)
    varid=netcdf.defVar(ncid,varname,'double',[dimid_lon dimid_lat dimid_time dimid_pft]);
else
    varid=netcdf.defVar(ncid,varname,'double',[dimid_lon dimid_lat dimid_time]);
end
netcdf.putAtt(ncid,varid,'Note',['Created from ',file_in,' using cable_monthly_to_annual_conv_func.m ',date])
netcdf.defVarDeflate(ncid,varid,true,true,9)
netcdf.close(ncid)

ncwrite([outdir,file_out],varname,data_out)


