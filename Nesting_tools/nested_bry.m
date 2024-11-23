function nested_bry( parent_grd, child_grd, parent_bry, child_bry,...
                     vertical_correc, extrapmask )

% Identify title
[pb1,pb2,pb3] = fileparts(parent_bry);
[cg1,cg2,cg3] = fileparts(child_grd);
title=['Boundary file for embedded grid: ', cg2,'.',cg3,', parent: ',pb2,'.',pb3];

parentgridnc = netcdf(parent_grd);
childgridnc = netcdf(child_grd);
imin=childgridnc{'grd_pos'}(1);
imax=childgridnc{'grd_pos'}(2);
jmin=childgridnc{'grd_pos'}(3);
jmax=childgridnc{'grd_pos'}(4);
refinecoeff=childgridnc{'refine_coef'}(:);
Lp = length(parentgridnc('xi_rho'));
Mp = length(parentgridnc('eta_rho'));
if extrapmask==1
  mask=parentgridnc{'mask_rho'}(:);
else
  mask=[];
end
L=Lp-1;
M=Mp-1;
N=length(parentgridnc{'s_rho'});
Np=N+1;

% Read parent bry file
parent_bry_nc = netcdf(parent_bry);
time = parent_bry_nc{'bry_time'}(:);
disp(' ')
disp(' Read in the parent boundary file...')
nc = netcdf(parent_clim);
theta_s = nc{'theta_s'}(:);
theta_b = nc{'theta_b'}(:);
Tcline = nc{'Tcline'}(:);
N = length(nc('s_rho'));
vtransform=nc{'Vtransform'}(:);
hc = nc{'hc'}(:);
disp([' Use parent VTRANSFORM = ',num2str(vtransform)])
if ~exist('vtransform') | isempty(vtransform)
    disp([' No VTRANSFORM parameter found'])
    disp([' Use the default one VTRANSFORM = 1'])
    vtransform=1;
end
ttime = nc{'tclm_time'}(:);
tcycle = nc{'tclm_time'}.cycle_length(:);
stime = nc{'sclm_time'}(:);
scycle = nc{'sclm_time'}.cycle_length(:);
utime = nc{'ssh_time'}(:);
ucycle = nc{'ssh_time'}.cycle_length(:);
vtime = nc{'ssh_time'}(:);
vcycle = nc{'ssh_time'}.cycle_length(:);
sshtime = nc{'ssh_time'}(:);
sshcycle = nc{'ssh_time'}.cycle_length(:);


% create bry file
%% obc, hc & vtransform are taken from crocotools_param
if  ~exist('vtransform')
  vtransform=1; %Old Vtransform
  disp([' NO VTRANSFORM parameter found'])
  disp([' USE TRANSFORM default value vtransform = 1'])
end
create_bryfile(child_bry, child_grd, title, obc, ...
  theta_s, theta_b, hc, N, ...
  time, parent_bry_nc{'bry_time'}.cycle_length, ...
  'clobber', vt)
child_bry_nc = netcdf(child_bry);

disp(' ')
disp(' Vertical interpolations')
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
ipchild=(imin:1/refinecoeff:imax);
jpchild=(jmin:1/refinecoeff:jmax);
irchild=(imin+0.5-0.5/refinecoeff:1/refinecoeff:imax+0.5+0.5/refinecoeff);
jrchild=(jmin+0.5-0.5/refinecoeff:1/refinecoeff:jmax+0.5+0.5/refinecoeff);
[ichildgrd_p,jchildgrd_p]=meshgrid(ipchild,jpchild);
[ichildgrd_r,jchildgrd_r]=meshgrid(irchild,jrchild);
[ichildgrd_u,jchildgrd_u]=meshgrid(ipchild,jrchild);
[ichildgrd_v,jchildgrd_v]=meshgrid(irchild,jpchild);

% Interpolations
disp("south")
disp("temp")
for tindex=1:length(bry_time)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(bry_time))])

end

return
