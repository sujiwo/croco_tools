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
[imin,imax,jmin,jmax,refinecoeff] = detect_grid(parent_grid_nc, child_grid_nc);
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
create_nestedbulk (child_blk, parent_blk, child_grd, title, bulkt, bulkc);
child_bulk_nc = netcdf(child_blk, 'write');

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
[ichildgrd_p,jchildgrd_p,...
ichildgrd_r,jchildgrd_r,...
ichildgrd_u,jchildgrd_u,...
ichildgrd_v,jchildgrd_v] = create_query_grid(imin,imax,jmin,jmax,refinecoeff,Lpc,Mpc);
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
