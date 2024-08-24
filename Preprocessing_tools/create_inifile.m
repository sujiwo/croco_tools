function create_inifile(inifile,gridfile,title,...
                         theta_s,theta_b,hc,N,time,clobber,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function nc=create_inifile(inifile,gridfile,theta_s,...
%                  theta_b,hc,N,ttime,stime,utime,... 
%                  cycle,clobber)
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input: 
% 
%   inifile      Netcdf initial file name (character string).
%   gridfile     Netcdf grid file name (character string).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)  
%   time         Initial time.(Real) 
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
disp([' Creating the file : ',inifile])
if nargin < 10
   disp([' NO VTRANSFORM parameter found'])
   disp([' USE TRANSFORM default value vtransform = 1'])
   vtransform = 1; 
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file
%
nc=nc4.netcdf(gridfile);
h=nc.var('h').get();
mask=nc.var('mask_rho').get();
nc.close();
hmin=min(min(h(mask==1)));
if vtransform ==1;
    if hc > hmin
        error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
    end
end
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'CROCO' ;
nc = netcdf(inifile,clobber);
%result = redef(nc);
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
nc.createDimension('time', 0);
nc.createDimension('one', 1);
%
%  Create variables
%
nc.createVariable('spherical', nc4.ncchar('one') );
nc.createVariable('Vtransform', nc4.ncint('one') );
nc.createVariable('Vstretching', nc4.ncint('one') );
nc.createVariable('tstart', nc4.ncdouble('one') );
nc.createVariable('tend', nc4.ncdouble('one') );
nc.createVariable('theta_s', nc4.ncdouble('one') );
nc.createVariable('theta_b', nc4.ncdouble('one') );
nc.createVariable('Tcline', nc4.ncdouble('one') );
nc.createVariable('hc', nc4.ncdouble('one') );
nc.createVariable('s_rho', nc4.ncdouble('s_rho') );
nc.createVariable('Cs_rho', nc4.ncdouble('s_rho') );
nc.createVariable('ocean_time', nc4.ncdouble('time') );
nc.createVariable('scrum_time', nc4.ncdouble('time') );
nc.createVariable('u', nc4.ncdouble('time','s_rho','eta_u','xi_u') );
nc.createVariable('v', nc4.ncdouble('time','s_rho','eta_v','xi_v') );
nc.createVariable('ubar', nc4.ncdouble('time','eta_u','xi_u') );
nc.createVariable('vbar', nc4.ncdouble('time','eta_v','xi_v') );
nc.createVariable('zeta', nc4.ncdouble('time','eta_rho','xi_rho') );
nc.createVariable('temp', nc4.ncdouble('time','s_rho','eta_rho','xi_rho') );
nc.createVariable('salt', nc4.ncdouble('time','s_rho','eta_rho','xi_rho') );
%
%  Create attributes
%
nc.var('Vtransform').createAttribute('long_name', 'vertical terrain-following transformation equation');
%
nc.var('Vstretching').createAttribute('long_name', 'vertical terrain-following stretching function');
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
nc.var('s_rho').createAttribute('units', 'nondimensional');
nc.var('s_rho').createAttribute('valid_min', -1);
nc.var('s_rho').createAttribute('valid_max', 0);
%
nc.var('Cs_rho').createAttribute('long_name', 'S-coordinate stretching curves at RHO-points');
nc.var('Cs_rho').createAttribute('units', 'nondimensional');
nc.var('Cs_rho').createAttribute('valid_min', -1);
nc.var('Cs_rho').createAttribute('valid_max', 0);
%
nc.var('ocean_time').createAttribute('long_name', 'time since initialization');
nc.var('ocean_time').createAttribute('units', 'second');
%
nc.var('scrum_time').createAttribute('long_name', 'time since initialization');
nc.var('scrum_time').createAttribute('units', 'second');
%
nc.var('u').createAttribute('long_name', 'u-momentum component');
nc.var('u').createAttribute('units', 'meter second-1');
%
nc.var('v').createAttribute('long_name', 'v-momentum component');
nc.var('v').createAttribute('units', 'meter second-1');
%
nc.var('ubar').createAttribute('long_name', 'vertically integrated u-momentum component');
nc.var('ubar').createAttribute('units', 'meter second-1');
%
nc.var('vbar').createAttribute('long_name', 'vertically integrated v-momentum component');
nc.var('vbar').createAttribute('units', 'meter second-1');
%
nc.var('zeta').createAttribute('long_name', 'free-surface');
nc.var('zeta').createAttribute('units', 'meter');
%
nc.var('temp').createAttribute('long_name', 'potential temperature');
nc.var('temp').createAttribute('units', 'Celsius');
%
nc.var('salt').createAttribute('long_name', 'salinity');
nc.var('salt').createAttribute('units', 'PSU');
%
% Create global attributes
%
nc.createAttribute('title', title);
nc.createAttribute('date', date);
nc.createAttribute('clim_file', inifile);
nc.createAttribute('grd_file', gridfile);
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
nc.var('tstart').set(time);
nc.var('tend').set(time);
nc.var('theta_s').set(theta_s);
nc.var('theta_b').set(theta_b);
nc.var('Tcline').set(hc);
nc.var('hc').set(hc);
nc.var('s_rho').set(s_rho);
nc.var('Cs_rho').set(Cs_rho);
nc.var('scrum_time').set(time*24*3600);
nc.var('ocean_time').set(time*24*3600);
nc.var('u').set(0);
nc.var('v').set(0);
nc.var('zeta').set(0);
nc.var('ubar').set(0);
nc.var('vbar').set(0);
nc.var('temp').set(0);
nc.var('salt').set(0);
%
% Synchronize on disk
%
nc.close();
return


