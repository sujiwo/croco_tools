function nested_forcing_any(parent_grd,child_grd,parent_frc,child_frc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  compute the forcing file of the embedded grid
%
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2004-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extrapmask=1;
%
% Title
%
title=['Forcing file for the embedded grid :',child_frc,...
' using parent forcing file: ',parent_frc];
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
parent_frc_nc = netcdf(parent_frc);

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
% Read in the parent forcing file
%
disp(' ')
disp(' Read in the parent forcing file...')

smst = parent_frc_nc{'sms_time'}(:);
smsc = parent_frc_nc{'sms_time'}.cycle_length(:);
shft = parent_frc_nc{'shf_time'}(:);
shfc = parent_frc_nc{'shf_time'}.cycle_length(:);
swft = parent_frc_nc{'swf_time'}(:);
swfc = parent_frc_nc{'swf_time'}.cycle_length(:);
sstt = parent_frc_nc{'sst_time'}(:);
sstc = parent_frc_nc{'sst_time'}.cycle_length(:);
ssst = parent_frc_nc{'sss_time'}(:);
sssc = parent_frc_nc{'sss_time'}.cycle_length(:);
srft = parent_frc_nc{'srf_time'}(:);
srfc = parent_frc_nc{'srf_time'}.cycle_length(:);
tide_period = parent_frc_nc{'tide_period'}(:)

%
% Create the forcing file
%
disp(' ')
disp(' Create the forcing file...')
create_nestedforcing(parent_frc_nc, ...
                     child_frc,parent_frc,child_grd,title,smst,...
                     shft,swft,srft,sstt,ssst,smsc,...
                     shfc,swfc,srfc,sstc,sssc)

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
%
% interpolations
% 
disp(' ')
disp(' Do the interpolations...')                 
np=netcdf(parent_frc);
nc=netcdf(child_frc,'write');
disp('the wind...')
for tindex=1:length(smst)
  interpvar3d(np,nc,igrd_u,jgrd_u,ichildgrd_u,jchildgrd_u,'sustr',mask,tindex)
  interpvar3d(np,nc,igrd_v,jgrd_v,ichildgrd_v,jchildgrd_v,'svstr',mask,tindex)
end
disp('heat flux...')
for tindex=1:length(shft)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'shflux',mask,tindex)
end
disp('freshwater flux...')
for tindex=1:length(swft)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'swflux',mask,tindex)
end
disp('SST and sst correction...')
for tindex=1:length(sstt)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'SST',mask,tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'dQdSST',mask,tindex)
end
disp('SSS...')
for tindex=1:length(ssst)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'SSS',mask,tindex)
end
disp('shortwave radiation...')
for tindex=1:length(srft)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'swrad',mask,tindex)
end
disp('tides...')
for tindex=1:length(tide_period)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Ephase',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Eamp',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Cmin',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Cmax',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Cangle',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Cphase',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Pamp',mask, tindex)
  interpvar3d(np,nc,igrd_r,jgrd_r,ichildgrd_r,jchildgrd_r,'tide_Pphase',mask, tindex)
end

close(parent_frc_nc);
result=close(nc);
disp(' ')
disp(' Done ')
%
% Make a plot
%
disp(' ')
disp(' Make a plot...')
figure(1)
plot_nestforcing(child_frc,'SSS',[1 6],1)
