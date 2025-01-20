function nested_bulk_any(parent_grd,child_grd,parent_blk,child_blk)

extrapmask=1;
%
% Title
%
title=['bulk file for the grid :',child_grd,...
' using parent bulk file: ',parent_blk];
disp(' ')
disp(title)

parent_grid_nc = netcdf(parent_grd);
child_grid_nc  = netcdf(child_grd);
parent_bulk_nc = netcdf(parent_blk);

Lp=length(parent_grid_nc('xi_rho'));
Mp=length(parent_grid_nc('eta_rho'));
Lpc=length(child_grid_nc('xi_rho'));
Mpc=length(child_grid_nc('eta_rho'));
if extrapmask==1
  mask=parent_grid_nc{'mask_rho'}(:);
else
  mask=[];
end
bulkt = parent_bulk_nc{'bulk_time'}(:);
bulkc = parent_bulk_nc{'bulk_time'}.cycle_length(:);

%
% Create the bulk file
%
disp(' ')
child_bulk_nc = netcdf(child_blk,'clobber');
redef(child_bulk_nc);
%
%  Create dimensions
%
child_bulk_nc('xi_rho') = Lpc;
child_bulk_nc('eta_rho') = Mpc;
child_bulk_nc('xi_psi') = Lpc-1;
child_bulk_nc('eta_psi') = Mpc-1;
child_bulk_nc('xi_u') = Lpc-1;
child_bulk_nc('eta_u') = Mpc;
child_bulk_nc('xi_v') = Lpc;
child_bulk_nc('eta_v') = Mpc-1;
child_bulk_nc('bulk_time') = length(bulkt);
%
%  Create variables and attributes
%
child_bulk_nc{'bulk_time'}              = ncdouble('bulk_time');
child_bulk_nc{'bulk_time'}.long_name    = ncchar('bulk formulation execution time');
child_bulk_nc{'bulk_time'}.long_name    = 'bulk formulation execution time';
child_bulk_nc{'bulk_time'}.units        = ncchar('days');
child_bulk_nc{'bulk_time'}.units        = 'days';
child_bulk_nc{'bulk_time'}.cycle_length = bulkc;

child_bulk_nc{'tair'}             = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'tair'}.long_name   = ncchar('surface air temperature');
child_bulk_nc{'tair'}.long_name   = 'surface air temperature';
child_bulk_nc{'tair'}.units       = ncchar('Celsius');
child_bulk_nc{'tair'}.units       = 'Celsius';

child_bulk_nc{'rhum'}             = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'rhum'}.long_name   = ncchar('relative humidity');
child_bulk_nc{'rhum'}.long_name   = 'relative humidity';
child_bulk_nc{'rhum'}.units       = ncchar('fraction');
child_bulk_nc{'rhum'}.units       = 'fraction';

child_bulk_nc{'prate'}            = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'prate'}.long_name  = ncchar('precipitation rate');
child_bulk_nc{'prate'}.long_name  = 'precipitation rate';
child_bulk_nc{'prate'}.units      = ncchar('cm day-1');
child_bulk_nc{'prate'}.units      = 'cm day-1';

child_bulk_nc{'wspd'}             = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'wspd'}.long_name   = ncchar('wind speed 10m');
child_bulk_nc{'wspd'}.long_name   = 'wind speed 10m';
child_bulk_nc{'wspd'}.units       = ncchar('m s-1');
child_bulk_nc{'wspd'}.units       = 'm s-1';

child_bulk_nc{'radlw'}            = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'radlw'}.long_name  = ncchar('net outgoing longwave radiation');
child_bulk_nc{'radlw'}.long_name  = 'net outgoing longwave radiation';
child_bulk_nc{'radlw'}.units      = ncchar('Watts meter-2');
child_bulk_nc{'radlw'}.units      = 'Watts meter-2';
child_bulk_nc{'radlw'}.positive   = ncchar('upward flux, cooling water');
child_bulk_nc{'radlw'}.positive   = 'upward flux, cooling water';

child_bulk_nc{'radlw_in'}            = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'radlw_in'}.long_name  = ncchar('downward longwave radiation');
child_bulk_nc{'radlw_in'}.long_name  = 'downward longwave radiation';
child_bulk_nc{'radlw_in'}.units      = ncchar('Watts meter-2');
child_bulk_nc{'radlw_in'}.units      = 'Watts meter-2';
child_bulk_nc{'radlw_in'}.positive   = ncchar('downward flux, cooling water');
child_bulk_nc{'radlw_in'}.positive   = 'downward flux, cooling water';

child_bulk_nc{'radsw'}            = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
child_bulk_nc{'radsw'}.long_name  = ncchar('solar shortwave radiation');
child_bulk_nc{'radsw'}.long_name  = 'shortwave radiation';
child_bulk_nc{'radsw'}.units      = ncchar('Watts meter-2');
child_bulk_nc{'radsw'}.units      = 'Watts meter-2';
child_bulk_nc{'radsw'}.positive   = ncchar('downward flux, heating water');
child_bulk_nc{'radsw'}.positive   = 'downward flux, heating water';

child_bulk_nc{'sustr'} = ncdouble('bulk_time', 'eta_u', 'xi_u');
child_bulk_nc{'sustr'}.long_name = ncchar('surface u-momentum stress');
child_bulk_nc{'sustr'}.long_name = 'surface u-momentum stress';
child_bulk_nc{'sustr'}.units = ncchar('Newton meter-2');
child_bulk_nc{'sustr'}.units = 'Newton meter-2';

child_bulk_nc{'svstr'} = ncdouble('bulk_time', 'eta_v', 'xi_v');
child_bulk_nc{'svstr'}.long_name = ncchar('surface v-momentum stress');
child_bulk_nc{'svstr'}.long_name = 'surface v-momentum stress';
child_bulk_nc{'svstr'}.units = ncchar('Newton meter-2');
child_bulk_nc{'svstr'}.units = 'Newton meter-2';

child_bulk_nc{'uwnd'} = ncdouble('bulk_time', 'eta_u', 'xi_u');
child_bulk_nc{'uwnd'}.long_name = ncchar('10m u-wind component');
child_bulk_nc{'uwnd'}.long_name = 'u-wind';
child_bulk_nc{'uwnd'}.units = ncchar('meter second-1');
child_bulk_nc{'uwnd'}.units = 'm/s';

child_bulk_nc{'vwnd'} = ncdouble('bulk_time', 'eta_v', 'xi_v');
child_bulk_nc{'vwnd'}.long_name = ncchar('10m v-wind component');
child_bulk_nc{'vwnd'}.long_name = 'v-wind';
child_bulk_nc{'vwnd'}.units = ncchar('meter second-1');
child_bulk_nc{'vwnd'}.units = 'm/s';

result = endef(child_bulk_nc);

%
% Create global attributes
%

child_bulk_nc.title = ncchar(title);
child_bulk_nc.title = title;
child_bulk_nc.date = ncchar(date);
child_bulk_nc.date = date;
child_bulk_nc.type = ncchar('CROCO bulk file');
child_bulk_nc.type = 'CROCO bulk file';
child_bulk_nc.grd_file = ncchar(child_grd);
child_bulk_nc.grd_file = child_grd;
child_bulk_nc.parent_grid = ncchar(parent_grd);
child_bulk_nc.parent_grid = parent_grd;
child_bulk_nc.parent_file = ncchar(parent_blk);
child_bulk_nc.parent_file = parent_blk;


%
% Write time variables
%

child_bulk_nc{'bulk_time'}(:) = bulkt;



%
% parent indices
%
[igrd_r,jgrd_r]=meshgrid((1:1:Lp),(1:1:Mp));
[igrd_p,jgrd_p]=meshgrid((1:1:Lp-1),(1:1:Mp-1));
[igrd_u,jgrd_u]=meshgrid((1:1:Lp-1),(1:1:Mp));
[igrd_v,jgrd_v]=meshgrid((1:1:Lp),(1:1:Mp-1));
%
% the children indices
%
[ichildgrd_r,jchildgrd_r]=meshgrid((1:1:Lpc),(1:1:Mpc));
[ichildgrd_p,jchildgrd_p]=meshgrid((1:1:Lpc-1),(1:1:Mpc-1));
[ichildgrd_u,jchildgrd_u]=meshgrid((1:1:Lpc-1),(1:1:Mpc));
[ichildgrd_v,jchildgrd_v]=meshgrid((1:1:Lpc),(1:1:Mpc-1));
%
% interpolations
%
disp(' ')
disp(' Do the interpolations...')

disp('tair...')
curvar = parent_bulk_nc{'tair'}(:);
varchild = zeros(size(child_bulk_nc{'tair'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'tair',mask,tindex);
end
child_bulk_nc{'tair'}(:)=varchild;

disp('rhum...')
curvar = parent_bulk_nc{'rhum'}(:);
varchild = zeros(size(child_bulk_nc{'rhum'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'rhum',mask,tindex);
end
child_bulk_nc{'rhum'}(:)=varchild;

disp('prate...')
curvar = parent_bulk_nc{'prate'}(:);
varchild = zeros(size(child_bulk_nc{'prate'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'prate',mask,tindex);
end
child_bulk_nc{'prate'}(:)=varchild;

disp('wspd...')
curvar = parent_bulk_nc{'wspd'}(:);
varchild = zeros(size(child_bulk_nc{'wspd'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'wspd',mask,tindex);
end
child_bulk_nc{'wspd'}(:)=varchild;

disp('radlw...')
curvar = parent_bulk_nc{'radlw'}(:);
varchild = zeros(size(child_bulk_nc{'radlw'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'radlw',mask,tindex);
end
child_bulk_nc{'radlw'}(:)=varchild;

disp('radlw_in...')
curvar = parent_bulk_nc{'radlw_in'}(:);
varchild = zeros(size(child_bulk_nc{'radlw_in'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'radlw_in',mask,tindex);
end
child_bulk_nc{'radlw_in'}(:)=varchild;

disp('radsw...')
curvar = parent_bulk_nc{'radsw'}(:);
varchild = zeros(size(child_bulk_nc{'radsw'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,...
      'radsw',mask,tindex);
end
child_bulk_nc{'radsw'}(:)=varchild;

disp('uwnd...')
curvar = parent_bulk_nc{'uwnd'}(:);
varchild = zeros(size(child_bulk_nc{'uwnd'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_u,jgrd_u,ichildgrd_u,jchildgrd_u,...
      'uwnd',mask,tindex);
end
child_bulk_nc{'uwnd'}(:)=varchild;

disp('vwnd...')
curvar = parent_bulk_nc{'vwnd'}(:);
varchild = zeros(size(child_bulk_nc{'vwnd'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_v,jgrd_v,ichildgrd_v,jchildgrd_v,...
      'vwnd',mask,tindex);
end
child_bulk_nc{'vwnd'}(:)=varchild;

disp('sustr...')
curvar = parent_bulk_nc{'sustr'}(:);
varchild = zeros(size(child_bulk_nc{'sustr'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_u,jgrd_u,ichildgrd_u,jchildgrd_u,...
      'sustr',mask,tindex);
end
child_bulk_nc{'sustr'}(:)=varchild;

disp('svstr...')
curvar = parent_bulk_nc{'svstr'}(:);
varchild = zeros(size(child_bulk_nc{'svstr'}));
parfor tindex=1:length(bulkt)
  varchild(tindex,:,:) = interpvar3d_par(curvar,...
      igrd_v,jgrd_v,ichildgrd_v,jchildgrd_v,...
      'svstr',mask,tindex);
end
child_bulk_nc{'svstr'}(:)=varchild;

close(parent_bulk_nc);
close(child_bulk_nc);
close(parent_grid_nc);
close(child_grid_nc);

disp(' ')
disp(' Done ')
%
% Make a plot
%
disp(' ')
disp(' Make a plot...')
figure(1)
plot_nestbulk(child_blk,'tair',[1 6],1)
figure(2)
plot_nestbulk(child_blk,'wspd',[1 6],1)
