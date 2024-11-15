function ncini=create_nestedinitial(inifile,gridfile,parentfile,title,...
				    theta_s,theta_b,Tcline,N,time,clobber,...
				    biol,pisces,namebiol,namepisces,...
				    unitbiol,unitpisces,hc,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function ncini=create_inifile(inifile,gridfile,theta_s,... 
%                  theta_b,Tcline,N,time,clobber) 
%
%   This function create the header of a Netcdf initial 
%   file.
%
%   Input: 
%
%   inifile     Netcdf initial file name (character string)
%   gridfile     Netcdf grid file name (character string).
%   parentfile  Netcdf parent initial file name (character string).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   Tcline       Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time        
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
%
%   Output
%
%   ncini       Output netcdf object.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp(['Creating the file : ',inifile])
disp(' ')
if nargin < 17
    disp([' NO VTRANSFORM parameter found'])
    disp([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1; 
    hc=[];
elseif nargin <18
    disp([' NO VTRANSFORM parameter found'])
    disp([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1; 
end
%
%  Read the grid file
%
nc=netcdf(gridfile, 'nowrite');
h=nc{'h'}(:); 
close(nc);
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'CROCO' ;
ncini = netcdf(inifile,clobber);
%%result = redef(ncini);
%
%  Create dimensions
%
ncini('xi_u') = L;
ncini('xi_v') = Lp;
ncini('xi_rho') = Lp;
ncini('eta_u') = Mp;
ncini('eta_v') = M;
ncini('eta_rho') = Mp;
ncini('s_rho') = N;
ncini('s_w') = Np;
ncini('tracer') = 2;
ncini('time') = length(time);
ncini('one') = 1;
%
%  Create variables
%
ncini{'spherical'} = ncchar('one') ;
ncini{'Vtransform'} = ncint('one') ;
ncini{'Vstretching'} = ncint('one') ;
ncini{'tstart'} = ncdouble('one') ;
ncini{'tend'} = ncdouble('one') ;
ncini{'theta_s'} = ncdouble('one') ;
ncini{'theta_b'} = ncdouble('one') ;
ncini{'Tcline'} = ncdouble('one') ;
ncini{'hc'} = ncdouble('one') ;
ncini{'s_rho'} = ncdouble('s_rho') ;
ncini{'Cs_rho'} = ncdouble('s_rho') ;
ncini{'scrum_time'} = ncdouble('time') ;
ncini{'u'} = ncdouble('time','s_rho','eta_u','xi_u') ;
ncini{'v'} = ncdouble('time','s_rho','eta_v','xi_v') ;
ncini{'ubar'} = ncdouble('time','eta_u','xi_u') ;
ncini{'vbar'} = ncdouble('time','eta_v','xi_v') ;
ncini{'zeta'} = ncdouble('time','eta_rho','xi_rho') ;
ncini{'temp'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
ncini{'salt'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
ncini{'hbl'} = ncdouble('time', 'eta_rho', 'xi_rho');
ncini{'hbbl'} = ncdouble('time', 'eta_rho', 'xi_rho');
ncini{'rufrc'} = ncdouble('time', 'eta_rho', 'xi_u');
ncini{'rvfrc'} = ncdouble('time', 'eta_v', 'xi_rho');
ncini{'dRdx'} = ncdouble('time', 's_rho', 'eta_rho', 'xi_u');
ncini{'dRde'} = ncdouble('time', 's_rho', 'eta_v', 'xi_rho');

if biol==1
    for k=1:length(namebiol)
        ncini{char(namebiol(k))}= ncdouble('time','s_rho','eta_rho','xi_rho');
    end
end
if pisces==1
    for k=1:length(namepisces)
        ncini{char(namepisces(k))}= ncdouble('time','s_rho','eta_rho','xi_rho');
    end
end
%
%  Create attributes
%
ncini{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
ncini{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
ncini{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
ncini{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
ncini{'tstart'}.long_name = ncchar('start processing day');
ncini{'tstart'}.long_name = 'start processing day';
ncini{'tstart'}.units = ncchar('day');
ncini{'tstart'}.units = 'day';
%
ncini{'tend'}.long_name = ncchar('end processing day');
ncini{'tend'}.long_name = 'end processing day';
ncini{'tend'}.units = ncchar('day');
ncini{'tend'}.units = 'day';
%
ncini{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
ncini{'theta_s'}.long_name = 'S-coordinate surface control parameter';
ncini{'theta_s'}.units = ncchar('nondimensional');
ncini{'theta_s'}.units = 'nondimensional';
%
ncini{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
ncini{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
ncini{'theta_b'}.units = ncchar('nondimensional');
ncini{'theta_b'}.units = 'nondimensional';
%
ncini{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
ncini{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
ncini{'Tcline'}.units = ncchar('meter');
ncini{'Tcline'}.units = 'meter';
%
ncini{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
ncini{'hc'}.long_name = 'S-coordinate parameter, critical depth';
ncini{'hc'}.units = ncchar('meter');
ncini{'hc'}.units = 'meter';
%
ncini{'s_rho'}.long_name = ncchar('S-coordinate at RHO-points');
ncini{'s_rho'}.long_name = 'S-coordinate at RHO-points';
ncini{'s_rho'}.units = ncchar('nondimensional');
ncini{'s_rho'}.units = 'nondimensional';
ncini{'s_rho'}.valid_min = -1;
ncini{'s_rho'}.valid_max = 0;
ncini{'s_rho'}.field = ncchar('s_rho, scalar');
ncini{'s_rho'}.field = 's_rho, scalar';
%
ncini{'Cs_rho'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
ncini{'Cs_rho'}.long_name = 'S-coordinate stretching curves at RHO-points';
ncini{'Cs_rho'}.units = ncchar('nondimensional');
ncini{'Cs_rho'}.units = 'nondimensional';
ncini{'Cs_rho'}.valid_min = -1;
ncini{'Cs_rho'}.valid_max = 0;
ncini{'Cs_rho'}.field = ncchar('Cs_rho, scalar');
ncini{'Cs_rho'}.field = 'Cs_rho, scalar';
%
ncini{'scrum_time'}.long_name = ncchar('time since initialization');
ncini{'scrum_time'}.long_name = 'time since initialization';
ncini{'scrum_time'}.units = ncchar('second');
ncini{'scrum_time'}.units = 'second';
ncini{'scrum_time'}.field = ncchar('time, scalar, series');
ncini{'scrum_time'}.field = 'time, scalar, series';
%
ncini{'u'}.long_name = ncchar('u-momentum component');
ncini{'u'}.long_name = 'u-momentum component';
ncini{'u'}.units = ncchar('meter second-1');
ncini{'u'}.units = 'meter second-1';
ncini{'u'}.field = ncchar('u-velocity, scalar, series');
ncini{'u'}.field = 'u-velocity, scalar, series';
%
ncini{'v'}.long_name = ncchar('v-momentum component');
ncini{'v'}.long_name = 'v-momentum component';
ncini{'v'}.units = ncchar('meter second-1');
ncini{'v'}.units = 'meter second-1';
ncini{'v'}.field = ncchar('v-velocity, scalar, series');
ncini{'v'}.field = 'v-velocity, scalar, series';
%
ncini{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
ncini{'ubar'}.long_name = 'vertically integrated u-momentum component';
ncini{'ubar'}.units = ncchar('meter second-1');
ncini{'ubar'}.units = 'meter second-1';
ncini{'ubar'}.field = ncchar('ubar-velocity, scalar, series');
ncini{'ubar'}.field = 'ubar-velocity, scalar, series';
%
ncini{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
ncini{'vbar'}.long_name = 'vertically integrated v-momentum component';
ncini{'vbar'}.units = ncchar('meter second-1');
ncini{'vbar'}.units = 'meter second-1';
ncini{'vbar'}.field = ncchar('vbar-velocity, scalar, series');
ncini{'vbar'}.field = 'vbar-velocity, scalar, series';
%
ncini{'zeta'}.long_name = ncchar('free-surface');
ncini{'zeta'}.long_name = 'free-surface';
ncini{'zeta'}.units = ncchar('meter');
ncini{'zeta'}.units = 'meter';
ncini{'zeta'}.field = ncchar('free-surface, scalar, series');
ncini{'zeta'}.field = 'free-surface, scalar, series';
%
ncini{'temp'}.long_name = ncchar('potential temperature');
ncini{'temp'}.long_name = 'potential temperature';
ncini{'temp'}.units = ncchar('Celsius');
ncini{'temp'}.units = 'Celsius';
ncini{'temp'}.field = ncchar('temperature, scalar, series');
ncini{'temp'}.field = 'temperature, scalar, series';
%
ncini{'salt'}.long_name = ncchar('salinity');
ncini{'salt'}.long_name = 'salinity';
ncini{'salt'}.units = ncchar('PSU');
ncini{'salt'}.units = 'PSU';
ncini{'salt'}.field = ncchar('salinity, scalar, series');
ncini{'salt'}.field = 'salinity, scalar, series';
%
ncini{'hbl'}.long_name = 'depth of planetary boundary layer';
ncini{'hbl'}.units = 'meter';
ncini{'hbl'}.field = 'hbl, scalar, series';
ncini{'hbl'}.standard_name = 'ocean_mixed_layer_thickness_defined_by_mixing_scheme';
ncini{'hbl'}.coordinates = 'lat_rho lon_rho';
%
ncini{'hbbl'}.long_name = 'depth of bottom boundary layer';
ncini{'hbbl'}.units = 'meter';
ncini{'hbbl'}.field = 'hbbl, scalar, series';
ncini{'hbbl'}.coordinates = 'lat_rho lon_rho';
%
ncini{'rufrc'}.long_name = 'barotropic forcing terms';
ncini{'rufrc'}.units = 'meter second-2';
ncini{'rufrc'}.field = 'rufrc, scalar, series';
ncini{'rufrc'}.standard_name = 'barotropic_forcing_terms';
ncini{'rufrc'}.coordinates = 'lat_u lon_u';
%
ncini{'rvfrc'}.long_name = 'barotropic forcing terms';
ncini{'rvfrc'}.units = 'meter second-2';
ncini{'rvfrc'}.field = 'rvfrc, scalar, series';
ncini{'rvfrc'}.standard_name = 'barotropic_forcing_terms';
ncini{'rvfrc'}.coordinates = 'lat_v lon_v';
%
ncini{'dRdx'}.long_name = 'Density gradient';
ncini{'dRdx'}.units = 'kg m-4';
ncini{'dRdx'}.field = 'dRdx, scalar, series';
ncini{'dRdx'}.standard_name = 'Density_gradient';
ncini{'dRdx'}.coordinates = 'lat_u lon_u';
%
ncini{'dRde'}.long_name = 'Density gradient';
ncini{'dRde'}.units = 'kg m-4';
ncini{'dRde'}.field = 'dRde, scalar, series';
ncini{'dRde'}.standard_name = 'Density_gradient';
ncini{'dRde'}.coordinates = 'lat_v lon_v';
%
if biol
  for k=1:length(namebiol)
%    disp(['K=',num2str(k)])
    ncini{char(namebiol(k))}.long_name = ncchar(char(namebiol(k)));
    ncini{char(namebiol(k))}.long_name = char(namebiol(k));
    ncini{char(namebiol(k))}.units = ncchar(char(unitbiol(k)));
    ncini{char(namebiol(k))}.units = char(unitbiol(k));
    ncini{char(namebiol(k))}.field = ncchar([char(namebiol(k)),', scalar, series']);
    ncini{char(namebiol(k))}.field = [char(namebiol(k)),', scalar, series'];
  end
end;
%
if pisces
  for k=1:length(namepisces)
 %   disp(['K=',num2str(k)])
    ncini{char(namepisces(k))}.long_name = ncchar(char(namepisces(k)));
    ncini{char(namepisces(k))}.long_name = char(namepisces(k));
    ncini{char(namepisces(k))}.units = ncchar(char(unitpisces(k)));
    ncini{char(namepisces(k))}.units = char(unitpisces(k));
    ncini{char(namepisces(k))}.field = ncchar([char(namepisces(k)),', scalar, series']);
    ncini{char(namepisces(k))}.field = [char(namepisces(k)),', scalar, series'];
  end
end;
%
% Create global attributes
%
disp('Create global attribute')
ncini.title = ncchar(title);
ncini.title = title;
ncini.date = ncchar(date);
ncini.date = date;
ncini.clim_file = ncchar(inifile);
ncini.clim_file = inifile;
ncini.grd_file = ncchar(gridfile);
ncini.grd_file = gridfile;
ncini.parent_file = ncchar(parentfile);
ncini.parent_file = parentfile;
ncini.type = ncchar(type);
ncini.type = type;
ncini.history = ncchar(history);
ncini.history = history;
ncini.creation_method = ncchar('Nestgui @ CROCOTOOLS');
ncini.creation_method = 'Nestgui @ CROCOTOOLS';
%
% Leave define mode
%
%%result = endef(ncini);
[s_rho,Cs_rho,s_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%
% Write variables
%
ncini{'spherical'}(:)='T';
ncini{'Vtransform'}(:)=vtransform;
ncini{'tstart'}(:)=time/(24*3600); 
ncini{'tend'}(:)=time/(24*3600); 
ncini{'theta_s'}(:)=theta_s; 
ncini{'theta_b'}(:)=theta_b; 
ncini{'Tcline'}(:)=Tcline; 
ncini{'hc'}(:)=hc; 
ncini{'s_rho'}(:)=s_rho; 
ncini{'Cs_rho'}(:)=Cs_rho; 
ncini{'scrum_time'}(1)=time; 
ncini{'u'}(:)=0; 
ncini{'v'}(:)=0; 
ncini{'zeta'}(:)=0; 
ncini{'ubar'}(:)=0; 
ncini{'vbar'}(:)=0; 
ncini{'temp'}(:)=0; 
ncini{'salt'}(:)=0; 
%
if biol==1
   % disp('Write variable biology')
    for k=1:length(namebiol)
        ncini{char(namebiol(k))}(:) = 0;
    end
end
%
if pisces==1
   % disp('Write variable pisces')
    for k=1:length(namepisces)
        ncini{char(namepisces(k))}(:) = 0;
    end
end
%
%
% Synchronize on disk
%
sync(ncini);
return


