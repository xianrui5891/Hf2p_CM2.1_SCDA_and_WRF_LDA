function ens_v(domains)
% MATLAB Script to read, interpolate, and write NetCDF data

% Define file paths
%  if nargin==0 || isempty(domains)
%         error('用法: ens_u(''d01'') 或 ens_u({''d01'',''d02'',''d03''})');
% end
% if ischar(domains); domains = {domains}; end
%     ok = ismember(domains, {'d01','d02','d03'});
% if ~all(ok)
%         error('域名必须是 d01/d02/d03，其它值: %s', strjoin(domains(~ok), ', '));
% end

% inputFile  = '/public2/wrf_teach/ens_read_write/wrfinput_d01';
inputFile = sprintf('/public2/wrf_teach/ens_read_write/wrfinput_%s', domains);
outputFile = '../wrfinput_d01';

% Read the output NetCDF file
lat_write = ncread(outputFile, 'XLAT_V');
lon_write = ncread(outputFile, 'XLONG_V');
znu_write = ncread(outputFile, 'ZNU');
t_write   = ncread(outputFile, 'V');
dims_write = size(t_write);

% Read the input NetCDF file
lat_read  = ncread(inputFile, 'XLAT_V');
lon_read  = ncread(inputFile, 'XLONG_V');
znu_read  = ncread(inputFile, 'ZNU');
t_read    = ncread(inputFile, 'V');
% time_read = ncread(inputFile, 'XTIME');
% Get dimensions
dims_read  = size(t_read);
% ntime_read = size(time_read, 1);

% Prepare output arrays
lat_4D = zeros( 1,dims_read(1), dims_read(2), dims_read(3));
lon_4D = zeros( 1,dims_read(1), dims_read(2), dims_read(3));
znu_4D = zeros( 1,dims_read(1), dims_read(2), dims_read(3));
lat_4D_write = zeros( dims_write(1), dims_write(2), dims_write(3));
lon_4D_write = zeros( dims_write(1), dims_write(2), dims_write(3));
znu_4D_write = zeros( dims_write(1), dims_write(2), dims_write(3));

% Fill in lat, lon, and znu arrays
for i = 1:dims_read(1)
    for j = 1:dims_read(2)
        for k = 1:dims_read(3)
            lat_4D(1, i, j, k) = lat_read(i, j);
            lon_4D(1, i, j, k) = lon_read(i, j);
            znu_4D(1, i, j, k) = znu_read(k);
        end
    end
end
for i = 1:dims_write(1)
    for j = 1:dims_write(2)
        for k = 1:dims_write(3)
            lat_4D_write( i, j, k) = lat_write(i, j);
            lon_4D_write( i, j, k) = lon_write(i, j);
            znu_4D_write( i, j, k) = znu_write(k,1);
        end
    end
end

% Reshape arrays to 1D
lat_1D = reshape(lat_4D, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
lon_1D = reshape(lon_4D, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
znu_1D = reshape(znu_4D, [dims_read(1) * dims_read(2) * dims_read(3), 1]);

% Initialize output temperature arrays
% t_out = zeros(1, dims_read(1), dims_read(2), dims_read(3));
lon_1D = double(lon_1D);
lat_1D = double(lat_1D);
znu_1D = double(znu_1D);
lon_4D_write = double(lon_4D_write);
lat_4D_write = double(lat_4D_write);
znu_4D_write = double(znu_4D_write);

% Loop over the temperature variables and interpolate

%======================= 5DAY T ==========================
t_out   = zeros( dims_write(1), dims_write(2), dims_write(3) , 20);
t_write = zeros( dims_write(1), dims_write(2), dims_write(3) , 1);
parpool(20);
parfor i = 1:20
    i
    T = ncread(inputFile, sprintf('ENS5DAYV%d', i));
    T_1D = reshape(T, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
    T_1D = double(T_1D);

    % Perform interpolation (example: griddata can be used)
    t_out( :, :, :,i) = griddata(lon_1D(:), lat_1D(:), znu_1D(:), T_1D(:), lon_4D_write, lat_4D_write, znu_4D_write ,  'nearest' );
end
delete(gcp);

% Write to output NetCDF file
ntime=size(1)
for i=1:20
    % nccreate(outputFile, sprintf('ENS5DAYV%d', i), 'Dimensions', {'west_east', dims_write(1) , 'south_north_stag', dims_write(2), 'bottom_top', dims_write(3),'Time', ntime(1) });
    t_write(:,:,:,1) = t_out( :, :, :,i);
    ncwrite(outputFile, sprintf('ENS5DAYV%d', i), t_write ) ;
    i
end

%======================= DAY T ==========================
t_out   = zeros( dims_write(1), dims_write(2), dims_write(3) , 20);
t_write = zeros( dims_write(1), dims_write(2), dims_write(3) , 1);
parpool(20);
parfor i = 1:20
    i
    T = ncread(inputFile, sprintf('ENS1DAYV%d', i));
    T_1D = reshape(T, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
    T_1D = double(T_1D);

    % Perform interpolation (example: griddata can be used)
    t_out( :, :, :,i) = griddata(lon_1D(:), lat_1D(:), znu_1D(:), T_1D(:), lon_4D_write, lat_4D_write, znu_4D_write ,  'nearest' );
end
delete(gcp);

% Write to output NetCDF file
ntime=size(1)
for i=1:20
    % nccreate(outputFile, sprintf('ENS1DAYV%d', i), 'Dimensions', {'west_east', dims_write(1) , 'south_north_stag', dims_write(2), 'bottom_top', dims_write(3),'Time', ntime(1) });
    t_write(:,:,:,1) = t_out( :, :, :,i);
    ncwrite(outputFile, sprintf('ENS1DAYV%d', i), t_write ) ;
    i
end

%======================= 6h T ==========================
t_out   = zeros( dims_write(1), dims_write(2), dims_write(3) , 20);
t_write = zeros( dims_write(1), dims_write(2), dims_write(3) , 1);
parpool(20);
parfor i = 1:20
    i
    T = ncread(inputFile, sprintf('ENS6HV%d', i));
    T_1D = reshape(T, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
    T_1D = double(T_1D);

    % Perform interpolation (example: griddata can be used)
    t_out( :, :, :,i) = griddata(lon_1D(:), lat_1D(:), znu_1D(:), T_1D(:), lon_4D_write, lat_4D_write, znu_4D_write ,  'nearest' );
end
delete(gcp);

% Write to output NetCDF file
ntime=size(1)
for i=1:20
    % nccreate(outputFile, sprintf('ENS6HV%d', i), 'Dimensions', {'west_east', dims_write(1) , 'south_north_stag', dims_write(2), 'bottom_top', dims_write(3),'Time', ntime(1) });
    t_write(:,:,:,1) = t_out( :, :, :,i);
    ncwrite(outputFile, sprintf('ENS6HV%d', i), t_write ) ;
    i
end


%======================= 3h T ==========================
t_out   = zeros( dims_write(1), dims_write(2), dims_write(3) , 20);
t_write = zeros( dims_write(1), dims_write(2), dims_write(3) , 1);
parpool(20);
parfor i = 1:20
    i
    T = ncread(inputFile, sprintf('ENS3HV%d', i));
    T_1D = reshape(T, [dims_read(1) * dims_read(2) * dims_read(3), 1]);
    T_1D = double(T_1D);

    % Perform interpolation (example: griddata can be used)
    t_out( :, :, :,i) = griddata(lon_1D(:), lat_1D(:), znu_1D(:), T_1D(:), lon_4D_write, lat_4D_write, znu_4D_write ,  'nearest' );
end
delete(gcp);

% Write to output NetCDF file
ntime=size(1)
for i=1:20
    % nccreate(outputFile, sprintf('ENS3HV%d', i), 'Dimensions', {'west_east', dims_write(1) , 'south_north_stag', dims_write(2), 'bottom_top', dims_write(3),'Time', ntime(1) });
    t_write(:,:,:,1) = t_out( :, :, :,i);
    ncwrite(outputFile, sprintf('ENS3HV%d', i), t_write ) ;
    i
end

disp('Data processing complete.');
end 

