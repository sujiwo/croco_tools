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
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
nc{'spherical'}.long_name = ncchar('grid type logical switch');
nc{'spherical'}.long_name = 'grid type logical switch';
nc{'spherical'}.flag_values = ncchar('T, F');
nc{'spherical'}.flag_values = 'T, F';
nc{'spherical'}.flag_meanings = ncchar('spherical Cartesian');
nc{'spherical'}.flag_meanings = 'spherical Cartesian';
%
nc{'tstart'}.long_name = ncchar('start processing day');
nc{'tstart'}.long_name = 'start processing day';
nc{'tstart'}.units = ncchar('day');
nc{'tstart'}.units = 'day';
%                   
nc{'tend'}.long_name = ncchar('end processing day');
nc{'tend'}.long_name = 'end processing day';
nc{'tend'}.units = ncchar('day');
nc{'tend'}.units = 'day';
%
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'s_rho'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'s_rho'}.long_name = 'S-coordinate at RHO-points';
nc{'s_rho'}.valid_min = -1.;
nc{'s_rho'}.valid_max = 0.;
nc{'s_rho'}.positive = ncchar('up');
nc{'s_rho'}.positive = 'up';
if (vtransform ==1)
    nc{'s_rho'}.standard_name = ncchar('ocean_s_coordinate_g1');
    nc{'s_rho'}.standard_name = 'ocean_s_coordinate_g1';
elseif (vtransform ==2)
    nc{'s_rho'}.standard_name = ncchar('ocean_s_coordinate_g2');
    nc{'s_rho'}.standard_name = 'ocean_s_coordinate_g2';     
end
nc{'s_rho'}.formula_terms = ncchar('s: s_rho C: Cs_rho eta: zeta depth: h depth_c: hc');
nc{'s_rho'}.formula_terms = 's: s_rho C: Cs_rho eta: zeta depth: h depth_c: hc';
%
nc{'s_w'}.long_name = ncchar('S-coordinate at W-points');
nc{'s_w'}.long_name = 'S-coordinate at W-points';
nc{'s_w'}.valid_min = -1. ;
nc{'s_w'}.valid_max = 0. ;
nc{'s_w'}.positive = ncchar('up');
nc{'s_w'}.positive = 'up';
if (vtransform == 1)
    nc{'s_w'}.standard_name = ncchar('ocean_s_coordinate_g1');
    nc{'s_w'}.standard_name = 'ocean_s_coordinate_g1';
elseif (vtransform == 2)
    nc{'s_w'}.standard_name = ncchar('ocean_s_coordinate_g2');
    nc{'s_w'}.standard_name = 'ocean_s_coordinate_g2';
end
nc{'s_w'}.formula_terms = ncchar('s: s_w C: Cs_w eta: zeta depth: h depth_c: hc');
nc{'s_w'}.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
%
nc{'Cs_rho'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_rho'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_rho'}.units = ncchar('nondimensional');
nc{'Cs_rho'}.units = 'nondimensional';
nc{'Cs_rho'}.valid_min = -1;
nc{'Cs_rho'}.valid_max = 0;
%
nc{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc{'Cs_w'}.long_name = 'S-coordinate stretching curves at W-points';
nc{'Cs_w'}.units = ncchar('nondimensional');
nc{'Cs_w'}.units = 'nondimensional';
nc{'Cs_w'}.valid_min = -1;
nc{'Cs_w'}.valid_max = 0;
%
nc{'tclm_time'}.long_name = ncchar('time for temperature climatology');
nc{'tclm_time'}.long_name = 'time for temperature climatology';
nc{'tclm_time'}.units = ncchar('day');
nc{'tclm_time'}.units = 'day';
nc{'tclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'tclm_time'}.calendar = '360.0 days in every year';
nc{'tclm_time'}.cycle_length = cycle;
%
nc{'temp_time'}.long_name = ncchar('time for temperature climatology');
nc{'temp_time'}.long_name = 'time for temperature climatology';
nc{'temp_time'}.units = ncchar('day');
nc{'temp_time'}.units = 'day';
nc{'temp_time'}.calendar = ncchar('360.0 days in every year');
nc{'temp_time'}.calendar = '360.0 days in every year';
nc{'temp_time'}.cycle_length = cycle;
%
nc{'sclm_time'}.long_name = ncchar('time for salinity climatology');
nc{'sclm_time'}.long_name = 'time for salinity climatology';
nc{'sclm_time'}.units = ncchar('day');
nc{'sclm_time'}.units = 'day';
nc{'sclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'sclm_time'}.calendar = '360.0 days in every year';
nc{'sclm_time'}.cycle_length = cycle;
%
nc{'salt_time'}.long_name = ncchar('time for salinity climatology');
nc{'salt_time'}.long_name = 'time for salinity climatology';
nc{'salt_time'}.units = ncchar('day');
nc{'salt_time'}.units = 'day';
nc{'salt_time'}.calendar = ncchar('360.0 days in every year');
nc{'salt_time'}.calendar = '360.0 days in every year';
nc{'salt_time'}.cycle_length = cycle;
%
nc{'uclm_time'}.long_name = ncchar('time climatological u');
nc{'uclm_time'}.long_name = 'time climatological u';
nc{'uclm_time'}.units = ncchar('day');
nc{'uclm_time'}.units = 'day';
nc{'uclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'uclm_time'}.calendar = '360.0 days in every year';
nc{'uclm_time'}.cycle_length = cycle;
%
nc{'vclm_time'}.long_name = ncchar('time climatological v');
nc{'vclm_time'}.long_name = 'time climatological v';
nc{'vclm_time'}.units = ncchar('day');
nc{'vclm_time'}.units = 'day';
nc{'vclm_time'}.calendar = ncchar('360.0 days in every year');
nc{'vclm_time'}.calendar = '360.0 days in every year';
nc{'vclm_time'}.cycle_length = cycle;
%
nc{'v2d_time'}.long_name = ncchar('time for 2D velocity climatology');
nc{'v2d_time'}.long_name = 'time for 2D velocity climatology';
nc{'v2d_time'}.units = ncchar('day');
nc{'v2d_time'}.units = 'day';
nc{'v2d_time'}.calendar = ncchar('360.0 days in every year');
nc{'v2d_time'}.calendar = '360.0 days in every year';
nc{'v2d_time'}.cycle_length = cycle;
%
nc{'v3d_time'}.long_name = ncchar('time for 3D velocity climatology');
nc{'v3d_time'}.long_name = 'time for 3D velocity climatology';
nc{'v3d_time'}.units = ncchar('day');
nc{'v3d_time'}.units = 'day';
nc{'v3d_time'}.calendar = ncchar('360.0 days in every year');
nc{'v3d_time'}.calendar = '360.0 days in every year';
nc{'v3d_time'}.cycle_length = cycle;
%
nc{'ssh_time'}.long_name = ncchar('time for sea surface height');
nc{'ssh_time'}.long_name = 'time for sea surface height';
nc{'ssh_time'}.units = ncchar('day');
nc{'ssh_time'}.units = 'day';
nc{'ssh_time'}.calendar = ncchar('360.0 days in every year');
nc{'ssh_time'}.calendar = '360.0 days in every year';
nc{'ssh_time'}.cycle_length = cycle;
%
nc{'zeta_time'}.long_name = ncchar('time for sea surface height');
nc{'zeta_time'}.long_name = 'time for sea surface height';
nc{'zeta_time'}.units = ncchar('day');
nc{'zeta_time'}.units = 'day';
nc{'zeta_time'}.calendar = ncchar('360.0 days in every year');
nc{'zeta_time'}.calendar = '360.0 days in every year';
nc{'zeta_time'}.cycle_length = cycle;
%
nc{'temp'}.long_name = ncchar('potential temperature');
nc{'temp'}.long_name = 'potential temperature';
nc{'temp'}.units = ncchar('Celsius');
nc{'temp'}.units = 'Celsius';
nc{'temp'}.time = ncchar('temp_time');
nc{'temp'}.time = 'temp_time';
nc{'temp'}.coordinates = ncchar('lon_rho lat_rho s_rho temp_time');
nc{'temp'}.coordinates = 'lon_rho lat_rho s_rho temp_time';
%
nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';
nc{'salt'}.time = ncchar('salt_time');
nc{'salt'}.time = 'salt_time';
nc{'salt'}.coordinates = ncchar('lon_rho lat_rho s_rho salt_time');
nc{'salt'}.coordinates = 'lon_rho lat_rho s_rho salt_time';
%
nc{'u'}.long_name = ncchar('u-momentum component');
nc{'u'}.long_name = 'u-momentum component';
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.units = 'meter second-1';
nc{'u'}.time = ncchar('uclm_time');
nc{'u'}.time = 'uclm_time';
nc{'u'}.coordinates = ncchar('lon_u lat_u s_rho uclm_time');
nc{'u'}.coordinates = 'lon_u lat_u s_rho u_time';
%
nc{'v'}.long_name = ncchar('v-momentum component');
nc{'v'}.long_name = 'v-momentum component';
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.units = 'meter second-1';
nc{'v'}.time = ncchar('vclm_time');
nc{'v'}.time = 'vclm_time';
nc{'v'}.coordinates = ncchar('lon_v lat_v s_rho vclm_time');
nc{'v'}.coordinates = 'lon_v lat_v s_rho vclm_time';
%
nc{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc{'ubar'}.long_name = 'vertically integrated u-momentum component';
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.units = 'meter second-1';
nc{'ubar'}.time = ncchar('uclm_time');
nc{'ubar'}.time = 'uclm_time';
nc{'ubar'}.coordinates = ncchar('lon_u lat_u uclm_time');
nc{'ubar'}.coordinates = 'lon_u lat_u uclm_time';
%
nc{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc{'vbar'}.long_name = 'vertically integrated v-momentum component';
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.units = 'meter second-1';
nc{'vbar'}.time = ncchar('vclm_time');
nc{'vbar'}.time = 'vclm_time';
nc{'vbar'}.coordinates = ncchar('lon_v lat_v vclm_time');
nc{'vbar'}.coordinates = 'lon_v lat_v vclm_time';
%
nc{'SSH'}.long_name = ncchar('sea surface height');
nc{'SSH'}.long_name = 'sea surface height';
nc{'SSH'}.units = ncchar('meter');
nc{'SSH'}.units = 'meter';
nc{'SSH'}.time = ncchar('zeta_time');
nc{'SSH'}.time = 'zeta_time';
nc{'SSH'}.coordinates = ncchar('lon_rho lat_rho zeta_time');
nc{'SSH'}.coordinates = 'lon_rho lat_rho zeta_time';
%
nc{'zeta'}.long_name = ncchar('sea surface height');
nc{'zeta'}.long_name = 'sea surface height';
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.units = 'meter';
nc{'zeta'}.time = ncchar('zeta_time');
nc{'zeta'}.time = 'zeta_time';
nc{'zeta'}.coordinates = ncchar('lon_rho lat_rho zeta_time');
nc{'zeta'}.coordinates = 'lon_rho lat_rho zeta_time';
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(clmname);
nc.clim_file = clmname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
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
nc{'spherical'}(:)='T';
nc{'Vtransform'}(:)=vtransform;
nc{'Vstretching'}(:)=1;
nc{'tstart'}(:) =  min([min(time) min(time) min(time)]); 
nc{'tend'}(:) =  max([max(time) max(time) max(time)]); 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'s_rho'}(:) = s_rho;
nc{'s_w'}(:) = s_w;
nc{'Cs_rho'}(:) =  Cs_rho; 
nc{'Cs_w'}(:) = Cs_w;
nc{'tclm_time'}(:) =  time; 
nc{'temp_time'}(:) =  time; 
nc{'sclm_time'}(:) =  time; 
nc{'salt_time'}(:) =  time; 
nc{'uclm_time'}(:) =  time; 
nc{'vclm_time'}(:) =  time; 
nc{'v2d_time'}(:) =   time; 
nc{'v3d_time'}(:) =   time; 
nc{'ssh_time'}(:) =   time;
nc{'zeta_time'}(:) =  time;
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'SSH'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
close(nc)
return


