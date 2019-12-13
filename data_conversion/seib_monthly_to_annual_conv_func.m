%Script to read in SEIB-DGVM monthly by pft data and convert it to annual for
%manageability
%
%T. Pugh
%25.04.16

indir='/home/pughtam/data/turnover/SEIB/IPSL/netcdf_v6/';
outdir='/home/pughtam/data/turnover/SEIB/IPSL/netcdf_v6/';

file_in='seib_ipsl_ffire_month_1901_2099.nc4';
file_out='seib_ipsl_ffire_annual_1901_2099.nc4';

varname='ffire';

ifpft=false; %Does the file have a PFT dimension?
npft=14; %Number of PFTs (if ifpft==true)

%---
time_in=ncread([indir,file_in],'time');
lon=ncread([indir,file_in],'longitude');
lat=ncread([indir,file_in],'latitude');
nlat=length(lat);
nlon=length(lon);
nmon=length(time_in);
nyear=nmon/12;

if ifpft
    data_out=NaN(nlon,nlat,nyear,npft);
    for yy=1:nyear
        mon_s=(yy*12)-11;
        data_in=ncread([indir,file_in],varname,[1 1 mon_s 1],[Inf Inf 12 Inf]);
        data_out(:,:,yy,:)=nansum(data_in,3)/12;
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
netcdf.putAtt(ncid,varid,'Note',['Created from ',file_in,' using seib_monthly_to_annual_conv_func.m ',date])
netcdf.defVarDeflate(ncid,varid,true,true,9)
netcdf.close(ncid)

ncwrite([outdir,file_out],varname,data_out)


