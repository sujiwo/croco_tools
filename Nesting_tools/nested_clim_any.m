function nested_clim_any(parent_grd,child_grd,parent_clim,child_clim,...
                     vertical_correc,extrapmask)

%
% Title
%
title=['Climatology file for the new grid :',child_clim,...
       ' using parent forcing file: ',parent_clim];
disp(' ')
disp(title)
if extrapmask==1
  disp('Extrapolation under mask is on')
  disp('====================')
end
%
if vertical_correc==1
  disp('Vertical correction is on')
  disp('===============')
end

parent_grid_nc = netcdf(parent_grd);
child_grid_nc  = netcdf(child_grd);
parent_clim_nc = netcdf(parent_clim);

Lp=length(parent_grid_nc('xi_rho'));
Mp=length(parent_grid_nc('eta_rho'));
Np=length(parent_grid_nc('s_rho'));
Lpc=length(child_grid_nc('xi_rho'));
Mpc=length(child_grid_nc('eta_rho'));
Npc=length(child_grid_nc('s_rho'));
[imin,imax,jmin,jmax,refinecoeff] = detect_grid(parent_grid_nc, child_grid_nc);
if extrapmask==1
  mask=parent_grid_nc{'mask_rho'}(:);
else
  mask=[];
end
%
% Read in the parent climatology file
%
disp(' ')
disp(' Read in the parent climatology file...')
%nc = netcdf(parent_clim);
theta_s = parent_clim_nc{'theta_s'}(:);
theta_b = parent_clim_nc{'theta_b'}(:);
Tcline = parent_clim_nc{'Tcline'}(:);
N = length(parent_clim_nc('s_rho'));
vtransform=parent_clim_nc{'Vtransform'}(:);
hc = parent_clim_nc{'hc'}(:);
disp([' Use parent VTRANSFORM = ',num2str(vtransform)])
if ~exist('vtransform') | isempty(vtransform)
    disp([' No VTRANSFORM parameter found'])
    disp([' Use the default one VTRANSFORM = 1'])
    vtransform=1;
end
ttime = parent_clim_nc{'tclm_time'}(:);
tcycle = parent_clim_nc{'tclm_time'}.cycle_length(:);
stime = parent_clim_nc{'sclm_time'}(:);
scycle = parent_clim_nc{'sclm_time'}.cycle_length(:);
utime = parent_clim_nc{'ssh_time'}(:);
ucycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
vtime = parent_clim_nc{'ssh_time'}(:);
vcycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
sshtime = parent_clim_nc{'ssh_time'}(:);
sshcycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
climtime=ttime;
if stime~=climtime | utime~=climtime | vtime~=climtime | sshtime~=climtime
  error('Nested_clim: different times... ')
end

%
% Create the climatology file
%
disp(' ')
disp(' Create the climatology file...')
child_clim_nc=create_nestedclim(child_clim,child_grd,parent_clim,title,...
			 theta_s,theta_b,Tcline,N,...
			 ttime,stime,utime,vtime,sshtime,...
			 tcycle,scycle,ucycle,vcycle,sshcycle,...
             0,0,0, 0, ...
             'clobber',...
			 0,0,...
             0, 0,0,0,...
             0,0,0,0,hc,vtransform);



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

disp('u...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar4d(parent_clim_nc,child_clim_nc,igrd_u,jgrd_u,ichildgrd_u,jchildgrd_u,'u',mask,tindex,N)
end
disp('v...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar4d(parent_clim_nc,child_clim_nc,igrd_v,jgrd_v,ichildgrd_v,jchildgrd_v,'v',mask,tindex,N)
end
disp('zeta...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar3d(parent_clim_nc,child_clim_nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'SSH',mask,tindex)
end
disp('ubar...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar3d(parent_clim_nc,child_clim_nc,igrd_u,jgrd_u,ichildgrd_u,jchildgrd_u,'ubar',mask,tindex)
end
disp('vbar...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar3d(parent_clim_nc,child_clim_nc,igrd_v,jgrd_v,ichildgrd_v,jchildgrd_v,'vbar',mask,tindex)
end
disp('temp...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar4d(parent_clim_nc,child_clim_nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'temp',mask,tindex,N)
end
disp('salt...')
for tindex=1:length(climtime)
  disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
  interpvar4d(parent_clim_nc,child_clim_nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'salt',mask,tindex,N)
end
close(parent_clim_nc);
close(child_clim_nc);

%
%
%  Vertical corrections
%
if (vertical_correc==1)
    disp('Process variable physical variables')
    for tindex=1:length(climtime)
        disp([' Time index : ',num2str(tindex),' of ',num2str(length(climtime))])
        vert_correc(child_clim,tindex,0,0,namebiol,namepisces)
    end
end

%
% Make a plot
%
disp(' ')
disp(' Make a plot...')
figure(1)
plot_nestclim(child_clim,child_grd,'temp',4)
return


close(child_clim_nc);
