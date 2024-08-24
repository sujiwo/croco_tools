function create_climfile(clmname,grdname,title,...
                         theta_s,theta_b,hc,N,...
                         time,cycle,clobber,vtransform);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_climfile(clmname,grdname,title,...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology
%   file.
%
%   Input:
%
%   clmname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time        Temperature climatology time.(vector)
%   time        Salinity climatology time.(vector)
%   time        Velocity climatology time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
%
%   Output
%
%   nc       Output netcdf object.
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
disp([' Creating the file : ',clmname])
disp(' ')
if nargin < 11
   disp([' NO VTRANSFORM parameter found'])
   disp([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1;
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file
%
nc = nc4.netcdf(grdname, 'nowrite');
h=nc.var('h').get;
maskr=nc.var('mask_rho').get();
Lp=length(nc.dim('xi_rho'));
Mp=length(nc.dim('eta_rho'));
nc.close();
hmin=min(min(h(maskr==1)));
if vtransform ==1;
  if hc > hmin
    error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
  end
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the climatology file
%
type = 'CLIMATOLOGY file' ;
history = 'CROCO' ;
nc = nc4.netcdf.create(clmname);
nc.redef();
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
%  Create variables
%
nc.createVariable('spherical', nc4.ncchar('one'));
nc.createVariable('Vtransform', nc4.ncint('one'));
nc.createVariable('Vstretching', nc4.ncint('one'));
nc.createVariable('tstart', nc4.ncdouble('one'));
nc.createVariable('tend', nc4.ncdouble('one'));
nc.createVariable('theta_s', nc4.ncdouble('one'));
nc.createVariable('theta_b', nc4.ncdouble('one'));
nc.createVariable('Tcline', nc4.ncdouble('one'));
nc.createVariable('hc', nc4.ncdouble('one'));
nc.createVariable('s_rho', nc4.ncdouble('s_rho'));
nc.createVariable('s_w', nc4.ncdouble('s_w'));
nc.createVariable('Cs_rho', nc4.ncdouble('s_rho'));
nc.createVariable('Cs_w', nc4.ncdouble('s_w'));
nc.createVariable('tclm_time', nc4.ncdouble('tclm_time'));
nc.createVariable('temp_time', nc4.ncdouble('temp_time'));
nc.createVariable('sclm_time', nc4.ncdouble('sclm_time'));
nc.createVariable('salt_time', nc4.ncdouble('salt_time'));
nc.createVariable('uclm_time', nc4.ncdouble('uclm_time'));
nc.createVariable('vclm_time', nc4.ncdouble('vclm_time'));
nc.createVariable('v2d_time', nc4.ncdouble('v2d_time'));
nc.createVariable('v3d_time', nc4.ncdouble('v3d_time'));
nc.createVariable('ssh_time', nc4.ncdouble('ssh_time'));
nc.createVariable('zeta_time', nc4.ncdouble('zeta_time'));
nc.createVariable('temp', nc4.ncdouble('tclm_time','s_rho','eta_rho','xi_rho'));
nc.createVariable('salt', nc4.ncdouble('sclm_time','s_rho','eta_rho','xi_rho'));
nc.createVariable('u', nc4.ncdouble('uclm_time','s_rho','eta_u','xi_u'));
nc.createVariable('v', nc4.ncdouble('vclm_time','s_rho','eta_v','xi_v'));
nc.createVariable('ubar', nc4.ncdouble('uclm_time','eta_u','xi_u'));
nc.createVariable('vbar', nc4.ncdouble('vclm_time','eta_v','xi_v'));
nc.createVariable('SSH', nc4.ncdouble('ssh_time','eta_rho','xi_rho'));
nc.createVariable('zeta', nc4.ncdouble('zeta_time','eta_rho','xi_rho'));
%
%  Create attributes
%
nc.var('Vtransform').createAttribute('long_name', 'vertical terrain-following transformation equation');
%

nc.var('Vstretching').createAttribute('long_name', 'vertical terrain-following stretching function');
%
nc.var('spherical').createAttribute('long_name', 'grid type logical switch');
nc.var('spherical').createAttribute('flag_values', 'T, F');
nc.var('spherical').createAttribute('flag_meanings', 'spherical Cartesian');
%
nc.var('tstart').createAttribute('long_name', 'start processing day');
nc.var('tstart').createAttribute('units', 'day');
%
nc.var('tend').createAttribute('long_name', 'end processing day');
nc.var('tend').createAttribute('units', 'day');
%
nc.var('theta_s').createAttribute('long_name', 'S-coordinate surface control parameter');
nc.var('theta_s').createAttribute('units', 'nondimensional');
%
nc.var('theta_b').createAttribute('long_name', 'S-coordinate bottom control parameter');
nc.var('theta_b').createAttribute('units', 'nondimensional');
%
nc.var('Tcline').createAttribute('long_name', 'S-coordinate surface/bottom layer width');
nc.var('Tcline').createAttribute('units', 'meter');
%
nc.var('hc').createAttribute('long_name', 'S-coordinate parameter, critical depth');
nc.var('hc').createAttribute('units', 'meter');
%
nc.var('s_rho').createAttribute('long_name', 'S-coordinate at RHO-points');
nc.var('s_rho').createAttribute('valid_min', -1.);
nc.var('s_rho').createAttribute('valid_max', 0.);
nc.var('s_rho').createAttribute('positive', 'up');
if (vtransform ==1)
    nc.var('s_rho').createAttribute('standard_name', 'ocean_s_coordinate_g1');
elseif (vtransform ==2)
    nc.var('s_rho').createAttribute('standard_name', 'ocean_s_coordinate_g2');
end
nc.var('s_rho').createAttribute('formula_terms', 's: s_rho C: Cs_rho eta: zeta depth: h depth_c: hc');
%
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
nc.var('Cs_rho').createAttribute('long_name', 'S-coordinate stretching curves at RHO-points');
nc.var('Cs_rho').createAttribute('units', 'nondimensional');
nc.var('Cs_rho').createAttribute('valid_min', -1);
nc.var('Cs_rho').createAttribute('valid_max', 0);
%
nc.var('Cs_w').createAttribute('long_name', 'S-coordinate stretching curves at W-points');
nc.var('Cs_w').createAttribute('units', 'nondimensional');
nc.var('Cs_w').createAttribute('valid_min', -1);
nc.var('Cs_w').createAttribute('valid_max', 0);
%
nc.var('tclm_time').createAttribute('long_name', 'time for temperature climatology');
nc.var('tclm_time').createAttribute('units', 'day');
nc.var('tclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('tclm_time').createAttribute('cycle_length', cycle);
%
nc.var('temp_time').createAttribute('long_name', 'time for temperature climatology');
nc.var('temp_time').createAttribute('units', 'day');
nc.var('temp_time').createAttribute('calendar', '360.0 days in every year');
nc.var('temp_time').createAttribute('cycle_length', cycle);
%
nc.var('sclm_time').createAttribute('long_name', 'time for salinity climatology');
nc.var('sclm_time').createAttribute('units', 'day');
nc.var('sclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('sclm_time').createAttribute('cycle_length', cycle);
%
nc.var('salt_time').createAttribute('long_name', 'time for salinity climatology');
nc.var('salt_time').createAttribute('units', 'day');
nc.var('salt_time').createAttribute('calendar', '360.0 days in every year');
nc.var('salt_time').createAttribute('cycle_length', cycle);
%
nc.var('uclm_time').createAttribute('long_name', 'time climatological u');
nc.var('uclm_time').createAttribute('units', 'day');
nc.var('uclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('uclm_time').createAttribute('cycle_length', cycle);
%
nc.var('vclm_time').createAttribute('long_name', 'time climatological v');
nc.var('vclm_time').createAttribute('units', 'day');
nc.var('vclm_time').createAttribute('calendar', '360.0 days in every year');
nc.var('vclm_time').createAttribute('cycle_length', cycle);
%
nc.var('v2d_time').createAttribute('long_name', 'time for 2D velocity climatology');
nc.var('v2d_time').createAttribute('units', 'day');
nc.var('v2d_time').createAttribute('calendar', '360.0 days in every year');
nc.var('v2d_time').createAttribute('cycle_length', cycle);
%
nc.var('v3d_time').createAttribute('long_name', 'time for 3D velocity climatology');
nc.var('v3d_time').createAttribute('units', 'day');
nc.var('v3d_time').createAttribute('calendar', '360.0 days in every year');
nc.var('v3d_time').createAttribute('cycle_length', cycle);
%
nc.var('ssh_time').createAttribute('long_name', 'time for sea surface height');
nc.var('ssh_time').createAttribute('units', 'day');
nc.var('ssh_time').createAttribute('calendar', '360.0 days in every year');
nc.var('ssh_time').createAttribute('cycle_length', cycle);
%
nc.var('zeta_time').createAttribute('long_name', 'time for sea surface height');
nc.var('zeta_time').createAttribute('units', 'day');
nc.var('zeta_time').createAttribute('calendar', '360.0 days in every year');
nc.var('zeta_time').createAttribute('cycle_length', cycle);
%
nc.var('temp').createAttribute('long_name', 'potential temperature');
nc.var('temp').createAttribute('units', 'Celsius');
nc.var('temp').createAttribute('time', 'temp_time');
nc.var('temp').createAttribute('coordinates', 'lon_rho lat_rho s_rho temp_time');
%
nc.var('salt').createAttribute('long_name', 'salinity');
nc.var('salt').createAttribute('units', 'PSU');
nc.var('salt').createAttribute('time', 'salt_time');
nc.var('salt').createAttribute('coordinates', 'lon_rho lat_rho s_rho salt_time');
%
nc.var('u').createAttribute('long_name', 'u-momentum component');
nc.var('u').createAttribute('units', 'meter second-1');
nc.var('u').createAttribute('time', 'uclm_time');
nc.var('u').createAttribute('coordinates', 'lon_u lat_u s_rho u_time');
%
nc.var('v').createAttribute('long_name', 'v-momentum component');
nc.var('v').createAttribute('units', 'meter second-1');
nc.var('v').createAttribute('time', 'vclm_time');
nc.var('v').createAttribute('coordinates', 'lon_v lat_v s_rho vclm_time');
%
nc.var('ubar').createAttribute('long_name', 'vertically integrated u-momentum component');
nc.var('ubar').createAttribute('units', 'meter second-1');
nc.var('ubar').createAttribute('time', 'uclm_time');
nc.var('ubar').createAttribute('coordinates', 'lon_u lat_u uclm_time');
%
nc.var('vbar').createAttribute('long_name', 'vertically integrated v-momentum component');
nc.var('vbar').createAttribute('units', 'meter second-1');
nc.var('vbar').createAttribute('time', 'vclm_time');
nc.var('vbar').createAttribute('coordinates', 'lon_v lat_v vclm_time');
%
nc.var('SSH').createAttribute('long_name', 'sea surface height');
nc.var('SSH').createAttribute('units', 'meter');
nc.var('SSH').createAttribute('time', 'zeta_time');
nc.var('SSH').createAttribute('coordinates', 'lon_rho lat_rho zeta_time');
%
nc.var('zeta').createAttribute('long_name', 'sea surface height');
nc.var('zeta').createAttribute('units', 'meter');
nc.var('zeta').createAttribute('time', 'zeta_time');
nc.var('zeta').createAttribute('coordinates', 'lon_rho lat_rho zeta_time');
%
% Create global attributes
%
nc.createAttribute('title', title);
nc.createAttribute('date', date);
nc.createAttribute('clim_file', clmname);
nc.createAttribute('grd_file', grdname);
nc.createAttribute('type', type);
nc.createAttribute('history', history);
%
% Leave define mode
%
%result = endef(nc);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
[s_rho,Cs_rho,s_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%disp(['vtransform=',num2str(vtransform)])

% cff1=1./sinh(theta_s);
% cff2=0.5/tanh(0.5*theta_s);
% s_rho=((1:N)-N-0.5)/N;
% Cs_rho=(1.-theta_b)*cff1*sinh(theta_s*s_rho)...
%     +theta_b*(cff2*tanh(theta_s*(s_rho+0.5))-0.5);
% s_w=((0:N)-N)/N;
% Cs_w=(1.-theta_b)*cff1*sinh(theta_s*s_w)...
%     +theta_b*(cff2*tanh(theta_s*(s_w+0.5))-0.5);


%
% Write variables
%
nc.var('spherical').set('T');
nc.var('Vtransform').set(vtransform);
nc.var('Vstretching').set(1);
nc.var('tstart').set(min([min(time) min(time) min(time)]) );
nc.var('tend').set(max([max(time) max(time) max(time)]));
nc.var('theta_s').set(theta_s);
nc.var('theta_b').set(theta_b);
nc.var('Tcline').set(hc);
nc.var('hc').set(hc);
nc.var('s_rho').set(s_rho);
nc.var('s_w').set(s_w);
nc.var('Cs_rho').set(Cs_rho);
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
nc.var('u').set(0);
nc.var('v').set(0);
nc.var('ubar').set(0);
nc.var('vbar').set(0);
nc.var('SSH').set(0);
nc.var('zeta').set(0);
nc.var('temp').set(0);
nc.var('salt').set(0);
nc.close();
return




