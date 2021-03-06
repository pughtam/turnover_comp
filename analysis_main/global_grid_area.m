function area=global_grid_area()
%Script to calculate area of every gridcell on a 0.5 x 0.5 degree global grid in
%m2
%
%T. Pugh
%14.10.13

minlat=-90.0;
maxlat=89.5;
minlon=-180.0;
maxlon=179.5;
grid=0.5;
offset=0.25;
surf_earth = 510108933.5E6; %in m2

lon_map=minlon:grid:maxlon;
lat_map=minlat:grid:maxlat;

nx=length(lon_map); % have to know this (0.5-degree-grid)
ny=length(lat_map);

area=NaN(ny,nx);
for x=1:nx
    for y=1:ny
        lon=lon_map(x)+offset; lat=lat_map(y)+offset; % needed if lat/lon gives SW corner, not centre
        x1 = round((lon+180)/360*nx)+1;
        x2 = x1;
        % compute area per output grid cell taking into account multiple cells
        xsin = sin(180/ny/2*pi/180);
        area(y,x)=surf_earth/nx*cos(lat*pi/180)*xsin*(mod(x2-x1+nx,nx)+1);
        clear lat lon x1 x2 xsin y
    end
end

