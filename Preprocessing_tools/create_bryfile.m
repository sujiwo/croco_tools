function create_bryfile(bryname,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,cycle,clobber,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer 
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
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
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')
if nargin < 12
    disp(' NO VTRANSFORM parameter found')
    disp(' USE TRANSFORM default value vtransform = 1')
    vtransform = 1; 
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file and check the topography
%
nc = nc4.netcdf(grdname, 'nowrite');
h=nc.var('h').get();
maskr=nc.var('mask_rho').get();
Lp=nc.dim('xi_rho').length();
Mp=nc.dim('eta_rho').length();
nc.close();
hmin=min(min(h(maskr==1)));
if vtransform ==1
  if hc > hmin
    error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
  end
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'CROCO' ;
nc = nc4.netcdf.create(bryname);
%%result = redef(nc);
%
%  Create dimensions
%
nc.createDimension('xi_u', L);
nc.createDimension('xi_v', Lp);
nc.createDimension('xi_rho', Lp);
nc.createDimension('eta_u', Mp);
nc.createDimension('eta_v', M);
nc.createDimension('eta_rho', Mp);
nc.createDimension('s_rho', N);
nc.createDimension('s_w', Np);
nc.createDimension('tracer', 2);
nc.createDimension('bry_time', length(time));
nc.createDimension('tclm_time', length(time));
nc.createDimension('temp_time', length(time));
nc.createDimension('sclm_time', length(time));
nc.createDimension('salt_time', length(time));
nc.createDimension('uclm_time', length(time));
nc.createDimension('vclm_time', length(time));
nc.createDimension('v2d_time', length(time));
nc.createDimension('v3d_time', length(time));
nc.createDimension('ssh_time', length(time));
nc.createDimension('zeta_time', length(time));
nc.createDimension('one', 1);
%
%  Create variables and attributes
%
nc.createVariable('spherical', nc4.ncchar('one') );
nc.var('spherical').createAttribute('long_name', 'grid type logical switch');
nc.var('spherical').createAttribute('flag_values', 'T, F');
nc.var('spherical').createAttribute('flag_meanings', 'spherical Cartesian');
%
nc.createVariable('Vtransform', nc4.ncint('one') );
nc.var('Vtransform').createAttribute('long_name', 'vertical terrain-following transformation equation');
%
nc.createVariable('Vstretching', nc4.ncint('one') );
nc.var('Vstretching').createAttribute('long_name', 'vertical terrain-following stretching function');
%
nc.createVariable('tstart', nc4.ncdouble('one') );
nc.var('tstart').createAttribute('long_name', 'start processing day');
nc.var('tstart').createAttribute('units', 'day');
%
nc.createVariable('tend', nc4.ncdouble('one') );
nc.var('tend').createAttribute('long_name', 'end processing day');
nc.var('tend').createAttribute('units', 'day');
%
nc.createVariable('theta_s', nc4.ncdouble('one') );
nc.var('theta_s').createAttribute('long_name', 'S-coordinate surface control parameter');
nc.var('theta_s').createAttribute('units', 'nondimensional');
%
nc.createVariable('theta_b', nc4.ncdouble('one') );
nc.var('theta_b').createAttribute('long_name', 'S-coordinate bottom control parameter');
nc.var('theta_b').createAttribute('units', 'nondimensional');
%
nc.createVariable('Tcline', nc4.ncdouble('one') );
nc.var('Tcline').createAttribute('long_name', 'S-coordinate surface/bottom layer width');
nc.var('Tcline').createAttribute('units', 'meter');
%
nc.createVariable('hc', nc4.ncdouble('one') );
nc.var('hc').createAttribute('long_name', 'S-coordinate parameter, critical depth');
nc.var('hc').createAttribute('units', 'meter');
%
nc.createVariable('s_rho', nc4.ncdouble('s_rho') );
nc.var('s_rho').createAttribute('long_name', 'S-coordinate at RHO-points');
nc.var('s_rho').createAttribute('valid_min', -1.);
nc.var('s_rho').createAttribute('valid_max', 0.);
nc.var('s_rho').createAttribute('positive', 'up');
if (vtransform == 1)
nc.var('s_rho').createAttribute('standard_name', 'ocean_s_coordinate_g1');
elseif (vtransform == 2)
nc.var('s_rho').createAttribute('standard_name', 'ocean_s_coordinate_g2');
end
nc.var('s_rho').createAttribute('formula_terms', 's: s_rho C: Cs_rho eta: zeta depth: h depth_c: hc');
%
nc.createVariable('s_w', nc4.ncdouble('s_w') );
nc.var('s_w').createAttribute('long_name', 'S-coordinate at W-points');
nc.var('s_w').createAttribute('valid_min', -1. );
nc.var('s_w').createAttribute('valid_max', 0. );
nc.var('s_w').createAttribute('positive', 'up');
if (vtransform == 1)
nc.var('s_w').createAttribute('standard_name', 'ocean_s_coordinate_g1');
elseif (vtransform == 2)
nc.var('s_w').createAttribute('standard_name', 'ocean_s_coordinate_g2');
end
nc.var('s_w').createAttribute('formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc');
%
nc.createVariable('Cs_rho', nc4.ncdouble('s_rho') );
nc.var('Cs_rho').createAttribute('long_name', 'S-coordinate stretching curves at RHO-points');
nc.var('Cs_rho').createAttribute('units', 'nondimensional');
nc.var('Cs_rho').createAttribute('valid_min', -1);
nc.var('Cs_rho').createAttribute('valid_max', 0);
%
nc.createVariable('Cs_w', nc4.ncdouble('s_w') );
nc.var('Cs_w').createAttribute('long_name', 'S-coordinate stretching curves at W-points');
nc.var('Cs_w').createAttribute('units', 'nondimensional');
nc.var('Cs_w').createAttribute('valid_min', -1);
nc.var('Cs_w').createAttribute('valid_max', 0);
%
nc.createVariable('bry_time', nc4.ncdouble('bry_time') );
nc.var('bry_time').createAttribute('long_name', 'time for boundary climatology');
nc.var('bry_time').createAttribute('units', 'day');
nc.var('bry_time').createAttribute('calendar', '360.0 days in every year');
nc.var('bry_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('tclm_time', nc4.ncdouble('tclm_time') );
nc.var('tclm_time').createAttribute('long_name', 'time for temperature climatology');
nc.var('tclm_time').createAttribute('units', 'day');
nc.var('tclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('tclm_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('temp_time', nc4.ncdouble('temp_time') );
nc.var('temp_time').createAttribute('long_name', 'time for temperature climatology');
nc.var('temp_time').createAttribute('units', 'day');
nc.var('temp_time').createAttribute('calendar', '360.0 days in every year');
nc.var('temp_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('sclm_time', nc4.ncdouble('sclm_time') );
nc.var('sclm_time').createAttribute('long_name', 'time for salinity climatology');
nc.var('sclm_time').createAttribute('units', 'day');
nc.var('sclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('sclm_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('salt_time', nc4.ncdouble('salt_time') );
nc.var('salt_time').createAttribute('long_name', 'time for salinity climatology');
nc.var('salt_time').createAttribute('units', 'day');
nc.var('salt_time').createAttribute('calendar', '360.0 days in every year');
nc.var('salt_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('uclm_time', nc4.ncdouble('uclm_time') );
nc.var('uclm_time').createAttribute('long_name', 'time climatological u');
nc.var('uclm_time').createAttribute('units', 'day');
nc.var('uclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('uclm_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('vclm_time', nc4.ncdouble('vclm_time') );
nc.var('vclm_time').createAttribute('long_name', 'time climatological v');
nc.var('vclm_time').createAttribute('units', 'day');
nc.var('vclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('vclm_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('v2d_time', nc4.ncdouble('v2d_time') );
nc.var('v2d_time').createAttribute('long_name', 'time for 2D velocity climatology');
nc.var('v2d_time').createAttribute('units', 'day');
nc.var('v2d_time').createAttribute('calendar', '360.0 days in every year');
nc.var('v2d_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('v3d_time', nc4.ncdouble('v3d_time') );
nc.var('v3d_time').createAttribute('long_name', 'time for 3D velocity climatology');
nc.var('v3d_time').createAttribute('units', 'day');
nc.var('v3d_time').createAttribute('calendar', '360.0 days in every year');
nc.var('v3d_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('ssh_time', nc4.ncdouble('ssh_time') );
nc.var('ssh_time').createAttribute('long_name', 'time for sea surface height');
nc.var('ssh_time').createAttribute('units', 'day');
nc.var('ssh_time').createAttribute('calendar', '360.0 days in every year');
nc.var('ssh_time').createAttribute('cycle_length', cycle);
%
nc.createVariable('zeta_time', nc4.ncdouble('zeta_time') );
nc.var('zeta_time').createAttribute('long_name', 'time for sea surface height');
nc.var('zeta_time').createAttribute('units', 'day');
nc.var('zeta_time').createAttribute('calendar', '360.0 days in every year');
nc.var('zeta_time').createAttribute('cycle_length', cycle);
%
if obc(1)==1
%
%   Southern boundary
%
  nc.createVariable('temp_south', nc4.ncdouble('temp_time','s_rho','xi_rho') );
nc.var('temp_south').createAttribute('long_name', 'southern boundary potential temperature');
nc.var('temp_south').createAttribute('units', 'Celsius');
nc.var('temp_south').createAttribute('coordinates', 'lon_rho s_rho temp_time');
%
  nc.createVariable('salt_south', nc4.ncdouble('salt_time','s_rho','xi_rho') );
nc.var('salt_south').createAttribute('long_name', 'southern boundary salinity');
nc.var('salt_south').createAttribute('units', 'PSU');
nc.var('salt_south').createAttribute('coordinates', 'lon_rho s_rho salt_time');
%
  nc.createVariable('u_south', nc4.ncdouble('v3d_time','s_rho','xi_u') );
nc.var('u_south').createAttribute('long_name', 'southern boundary u-momentum component');
nc.var('u_south').createAttribute('units', 'meter second-1');
nc.var('u_south').createAttribute('coordinates', 'lon_u s_rho u_time');
%
  nc.createVariable('v_south', nc4.ncdouble('v3d_time','s_rho','xi_rho') );
nc.var('v_south').createAttribute('long_name', 'southern boundary v-momentum component');
nc.var('v_south').createAttribute('units', 'meter second-1');
nc.var('v_south').createAttribute('coordinates', 'lon_v s_rho vclm_time');
%
  nc.createVariable('ubar_south', nc4.ncdouble('v2d_time','xi_u') );
nc.var('ubar_south').createAttribute('long_name', 'southern boundary vertically integrated u-momentum component');
nc.var('ubar_south').createAttribute('units', 'meter second-1');
nc.var('ubar_south').createAttribute('coordinates', 'lon_u uclm_time');
%
  nc.createVariable('vbar_south', nc4.ncdouble('v2d_time','xi_rho') );
nc.var('vbar_south').createAttribute('long_name', 'southern boundary vertically integrated v-momentum component');
nc.var('vbar_south').createAttribute('units', 'meter second-1');
nc.var('vbar_south').createAttribute('coordinates', 'lon_v vclm_time');
%
  nc.createVariable('zeta_south', nc4.ncdouble('zeta_time','xi_rho') );
nc.var('zeta_south').createAttribute('long_name', 'southern boundary sea surface height');
nc.var('zeta_south').createAttribute('units', 'meter');
nc.var('zeta_south').createAttribute('coordinates', 'lon_rho zeta_time');
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  nc.createVariable('temp_east', nc4.ncdouble('temp_time','s_rho','eta_rho') );
nc.var('temp_east').createAttribute('long_name', 'eastern boundary potential temperature');
nc.var('temp_east').createAttribute('units', 'Celsius');
nc.var('temp_east').createAttribute('coordinates', 'lat_rho s_rho temp_time');
%
  nc.createVariable('salt_east', nc4.ncdouble('salt_time','s_rho','eta_rho') );
nc.var('salt_east').createAttribute('long_name', 'eastern boundary salinity');
nc.var('salt_east').createAttribute('units', 'PSU');
nc.var('salt_east').createAttribute('coordinates', 'lat_rho s_rho salt_time');
%
  nc.createVariable('u_east', nc4.ncdouble('v3d_time','s_rho','eta_rho') );
nc.var('u_east').createAttribute('long_name', 'eastern boundary u-momentum component');
nc.var('u_east').createAttribute('units', 'meter second-1');
nc.var('u_east').createAttribute('coordinates', 'lat_u s_rho u_time');
%
  nc.createVariable('v_east', nc4.ncdouble('v3d_time','s_rho','eta_v') );
nc.var('v_east').createAttribute('long_name', 'eastern boundary v-momentum component');
nc.var('v_east').createAttribute('units', 'meter second-1');
nc.var('v_east').createAttribute('coordinates', 'lat_v s_rho vclm_time');
%
  nc.createVariable('ubar_east', nc4.ncdouble('v2d_time','eta_rho') );
nc.var('ubar_east').createAttribute('long_name', 'eastern boundary vertically integrated u-momentum component');
nc.var('ubar_east').createAttribute('units', 'meter second-1');
nc.var('ubar_east').createAttribute('coordinates', 'lat_u uclm_time');
%
  nc.createVariable('vbar_east', nc4.ncdouble('v2d_time','eta_v') );
nc.var('vbar_east').createAttribute('long_name', 'eastern boundary vertically integrated v-momentum component');
nc.var('vbar_east').createAttribute('units', 'meter second-1');
nc.var('vbar_east').createAttribute('coordinates', 'lat_v vclm_time');
%
  nc.createVariable('zeta_east', nc4.ncdouble('zeta_time','eta_rho') );
nc.var('zeta_east').createAttribute('long_name', 'eastern boundary sea surface height');
nc.var('zeta_east').createAttribute('units', 'meter');
nc.var('zeta_east').createAttribute('coordinates', 'lat_rho zeta_time');
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  nc.createVariable('temp_north', nc4.ncdouble('temp_time','s_rho','xi_rho') );
nc.var('temp_north').createAttribute('long_name', 'northern boundary potential temperature');
nc.var('temp_north').createAttribute('units', 'Celsius');
nc.var('temp_north').createAttribute('coordinates', 'lon_rho s_rho temp_time');
%
  nc.createVariable('salt_north', nc4.ncdouble('salt_time','s_rho','xi_rho') );
nc.var('salt_north').createAttribute('long_name', 'northern boundary salinity');
nc.var('salt_north').createAttribute('units', 'PSU');
nc.var('salt_north').createAttribute('coordinates', 'lon_rho s_rho salt_time');
%
  nc.createVariable('u_north', nc4.ncdouble('v3d_time','s_rho','xi_u') );
nc.var('u_north').createAttribute('long_name', 'northern boundary u-momentum component');
nc.var('u_north').createAttribute('units', 'meter second-1');
nc.var('u_north').createAttribute('coordinates', 'lon_u s_rho u_time');
%
  nc.createVariable('v_north', nc4.ncdouble('v3d_time','s_rho','xi_rho') );
nc.var('v_north').createAttribute('long_name', 'northern boundary v-momentum component');
nc.var('v_north').createAttribute('units', 'meter second-1');
nc.var('v_north').createAttribute('coordinates', 'lon_v s_rho vclm_time');
%
  nc.createVariable('ubar_north', nc4.ncdouble('v2d_time','xi_u') );
nc.var('ubar_north').createAttribute('long_name', 'northern boundary vertically integrated u-momentum component');
nc.var('ubar_north').createAttribute('units', 'meter second-1');
nc.var('ubar_north').createAttribute('coordinates', 'lon_u uclm_time');
%
  nc.createVariable('vbar_north', nc4.ncdouble('v2d_time','xi_rho') );
nc.var('vbar_north').createAttribute('long_name', 'northern boundary vertically integrated v-momentum component');
nc.var('vbar_north').createAttribute('units', 'meter second-1');
nc.var('vbar_north').createAttribute('coordinates', 'lon_v vclm_time');

  nc.createVariable('zeta_north', nc4.ncdouble('zeta_time','xi_rho') );
nc.var('zeta_north').createAttribute('long_name', 'northern boundary sea surface height');
nc.var('zeta_north').createAttribute('units', 'meter');
nc.var('zeta_north').createAttribute('coordinates', 'lon_rho zeta_time');
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc.createVariable('temp_west', nc4.ncdouble('temp_time','s_rho','eta_rho') );
nc.var('temp_west').createAttribute('long_name', 'western boundary potential temperature');
nc.var('temp_west').createAttribute('units', 'Celsius');
nc.var('temp_west').createAttribute('coordinates', 'lat_rho s_rho temp_time');
%
  nc.createVariable('salt_west', nc4.ncdouble('salt_time','s_rho','eta_rho') );
nc.var('salt_west').createAttribute('long_name', 'western boundary salinity');
nc.var('salt_west').createAttribute('units', 'PSU');
nc.var('salt_west').createAttribute('coordinates', 'lat_rho s_rho salt_time');
%
  nc.createVariable('u_west', nc4.ncdouble('v3d_time','s_rho','eta_rho') );
nc.var('u_west').createAttribute('long_name', 'western boundary u-momentum component');
nc.var('u_west').createAttribute('units', 'meter second-1');
nc.var('u_west').createAttribute('coordinates', 'lat_u s_rho u_time');
%
  nc.createVariable('v_west', nc4.ncdouble('v3d_time','s_rho','eta_v') );
nc.var('v_west').createAttribute('long_name', 'western boundary v-momentum component');
nc.var('v_west').createAttribute('units', 'meter second-1');
nc.var('v_west').createAttribute('coordinates', 'lat_v s_rho vclm_time');
%
  nc.createVariable('ubar_west', nc4.ncdouble('v2d_time','eta_rho') );
nc.var('ubar_west').createAttribute('long_name', 'western boundary vertically integrated u-momentum component');
nc.var('ubar_west').createAttribute('units', 'meter second-1');
nc.var('ubar_west').createAttribute('coordinates', 'lat_u uclm_time');
%
  nc.createVariable('vbar_west', nc4.ncdouble('v2d_time','eta_v') );
nc.var('vbar_west').createAttribute('long_name', 'western boundary vertically integrated v-momentum component');
nc.var('vbar_west').createAttribute('units', 'meter second-1');
nc.var('vbar_west').createAttribute('coordinates', 'lat_v vclm_time');
%
  nc.createVariable('zeta_west', nc4.ncdouble('zeta_time','eta_rho') );
nc.var('zeta_west').createAttribute('long_name', 'western boundary sea surface height');
nc.var('zeta_west').createAttribute('units', 'meter');
nc.var('zeta_west').createAttribute('coordinates', 'lat_rho zeta_time');
%
end
%
%
% Create global attributes
%
nc.createAttribute('title', title);
nc.createAttribute('date', date);
nc.createAttribute('clim_file', bryname);
nc.createAttribute('grd_file', grdname);
nc.createAttribute('type', type);
nc.createAttribute('history', history);
%
% Leave define mode
%
%%result = endef(nc);
%
% Compute S coordinates
%
[s_rho,Cs_rho,s_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%disp(['vtransform=',num2str(vtransform)])
%
% Write variables
%
nc.var('spherical').set('T');
nc.var('Vtransform').set(vtransform);
nc.var('Vstretching').set(1);
nc.var('tstart').set(min([min(time) min(time) min(time)]));
nc.var('tend').set(max([max(time) max(time) max(time)]));
nc.var('theta_s').set(theta_s);
nc.var('theta_b').set(theta_b);
nc.var('Tcline').set(hc);
nc.var('hc').set(hc);
nc.var('s_rho').set(s_rho);
nc.var('s_w').set(s_w);
nc.var('Cs_rho').set(Cs_rho );
nc.var('Cs_w').set(Cs_w);
nc.var('tclm_time').set(time);
nc.var('temp_time').set(time);
nc.var('sclm_time').set(time);
nc.var('salt_time').set(time);
nc.var('uclm_time').set(time);
nc.var('vclm_time').set(time);
nc.var('v2d_time').set(time);
nc.var('v3d_time').set(time);
nc.var('ssh_time').set(time);
nc.var('zeta_time').set(time);
nc.var('bry_time').set(time);
if obc(1)==1
  nc.var('u_south').set(0);
  nc.var('v_south').set(0);
  nc.var('ubar_south').set(0);
  nc.var('vbar_south').set(0);
  nc.var('zeta_south').set(0);
  nc.var('temp_south').set(0);
  nc.var('salt_south').set(0);
end 
if obc(2)==1
  nc.var('u_east').set(0);
  nc.var('v_east').set(0);
  nc.var('ubar_east').set(0);
  nc.var('vbar_east').set(0);
  nc.var('zeta_east').set(0);
  nc.var('temp_east').set(0);
  nc.var('salt_east').set(0);
end 
if obc(3)==1
  nc.var('u_north').set(0);
  nc.var('v_north').set(0);
  nc.var('ubar_north').set(0);
  nc.var('vbar_north').set(0);
  nc.var('zeta_north').set(0);
  nc.var('temp_north').set(0);
  nc.var('salt_north').set(0);
end 
if obc(4)==1
  nc.var('u_west').set(0);
  nc.var('v_west').set(0);
  nc.var('ubar_west').set(0);
  nc.var('vbar_west').set(0);
  nc.var('zeta_west').set(0);
  nc.var('temp_west').set(0);
  nc.var('salt_west').set(0);
end 
nc.close();
return


