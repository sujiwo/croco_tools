function  create_nestedforcing(parentnc, ...
                         frcname,parentname,grdname,title,smst,...
                         shft,swft,srft,sstt,ssst,smsc,...
                         shfc,swfc,srfc,sstc,sssc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf forcing file
%       frcname: name of the forcing file
%       grdname: name of the grid file
%       title: title in the netcdf file  
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
nc=netcdf(grdname);
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
close(nc);
Lp=L+1;
Mp=M+1;
nw = netcdf(frcname, 'clobber');
redef(nw);

tide_period = parentnc{'tide_period'}(:);
tide_Ephase = parentnc{'tide_Ephase'}(:);
tide_Eamp = parentnc{'tide_Eamp'}(:);
tide_Cmin = parentnc{'tide_Cmin'}(:);
tide_Cmax = parentnc{'tide_Cmax'}(:);
tide_Cangle = parentnc{'tide_Cangle'}(:);
tide_Cphase = parentnc{'tide_Cphase'}(:);
tide_Pamp = parentnc{'tide_Pamp'}(:);
tide_Pphase = parentnc{'tide_Pphase'}(:);

%
%  Create dimensions
%
nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('sms_time') = length(smst);
nw('shf_time') = length(shft);
nw('swf_time') = length(swft);
nw('sst_time') = length(sstt);
nw('srf_time') = length(srft);
nw('sss_time') = length(ssst);
nw('tide_period') = length(tide_period);

%
%  Create variables and attributes
%
nw{'sms_time'} = ncdouble('sms_time');
nw{'sms_time'}.long_name = ncchar('surface momentum stress time');
nw{'sms_time'}.long_name = 'surface momentum stress time';
nw{'sms_time'}.units = ncchar('days');
nw{'sms_time'}.units = 'days';
nw{'sms_time'}.cycle_length = smsc;
nw{'sms_time'}.field = ncchar('time, scalar, series');
nw{'sms_time'}.field = 'time, scalar, series';

nw{'shf_time'} = ncdouble('shf_time');
nw{'shf_time'}.long_name = ncchar('surface heat flux time');
nw{'shf_time'}.long_name = 'surface heat flux time';
nw{'shf_time'}.units = ncchar('days');
nw{'shf_time'}.units = 'days';
nw{'shf_time'}.cycle_length =shfc ;
nw{'shf_time'}.field = ncchar('time, scalar, series');
nw{'shf_time'}.field = 'time, scalar, series';

nw{'swf_time'} = ncdouble('swf_time');
nw{'swf_time'}.long_name = ncchar('surface freshwater flux time');
nw{'swf_time'}.long_name = 'surface freshwater flux time';
nw{'swf_time'}.units = ncchar('days');
nw{'swf_time'}.units = 'days';
nw{'swf_time'}.cycle_length = swfc;
nw{'swf_time'}.field = ncchar('time, scalar, series');
nw{'swf_time'}.field = 'time, scalar, series';


nw{'sst_time'} = ncdouble('sst_time');
nw{'sst_time'}.long_name = ncchar('sea surface temperature time');
nw{'sst_time'}.long_name = 'sea surface temperature time';
nw{'sst_time'}.units = ncchar('days');
nw{'sst_time'}.units = 'days';
nw{'sst_time'}.cycle_length = sstc;
nw{'sst_time'}.field = ncchar('time, scalar, series');
nw{'sst_time'}.field = 'time, scalar, series';

nw{'sss_time'} = ncdouble('sss_time');
nw{'sss_time'}.long_name = ncchar('sea surface salinity time');
nw{'sss_time'}.long_name = 'sea surface salinity time';
nw{'sss_time'}.units = ncchar('days');
nw{'sss_time'}.units = 'days';
nw{'sss_time'}.cycle_length = sssc;
nw{'sss_time'}.field = ncchar('time, scalar, series');
nw{'sss_time'}.field = 'time, scalar, series';

nw{'srf_time'} = ncdouble('srf_time');
nw{'srf_time'}.long_name = ncchar('solar shortwave radiation time');
nw{'srf_time'}.long_name = 'solar shortwave radiation time';
nw{'srf_time'}.units = ncchar('days');
nw{'srf_time'}.units = 'days';
nw{'srf_time'}.cycle_length = srfc;
nw{'srf_time'}.field = ncchar('time, scalar, series');
nw{'srf_time'}.field = 'time, scalar, series';

nw{'sustr'} = ncdouble('sms_time', 'eta_u', 'xi_u');
nw{'sustr'}.long_name = ncchar('surface u-momentum stress');
nw{'sustr'}.long_name = 'surface u-momentum stress';
nw{'sustr'}.units = ncchar('Newton meter-2');
nw{'sustr'}.units = 'Newton meter-2';
nw{'sustr'}.field = ncchar('surface u-mometum stress, scalar, series');
nw{'sustr'}.field = 'surface u-mometum stress, scalar, series';

nw{'svstr'} = ncdouble('sms_time', 'eta_v', 'xi_v');
nw{'svstr'}.long_name = ncchar('surface v-momentum stress');
nw{'svstr'}.long_name = 'surface v-momentum stress';
nw{'svstr'}.units = ncchar('Newton meter-2');
nw{'svstr'}.units = 'Newton meter-2';
nw{'svstr'}.field = ncchar('surface v-momentum stress, scalar, series');
nw{'svstr'}.field = 'surface v-momentum stress, scalar, series';

nw{'shflux'} = ncdouble('shf_time', 'eta_rho', 'xi_rho');
nw{'shflux'}.long_name = ncchar('surface net heat flux');
nw{'shflux'}.long_name = 'surface net heat flux';
nw{'shflux'}.units = ncchar('Watts meter-2');
nw{'shflux'}.units = 'Watts meter-2';
nw{'shflux'}.field = ncchar('surface heat flux, scalar, series');
nw{'shflux'}.field = 'surface heat flux, scalar, series';

nw{'swflux'} = ncdouble('swf_time', 'eta_rho', 'xi_rho');
nw{'swflux'}.long_name = ncchar('surface freshwater flux (E-P)');
nw{'swflux'}.long_name = 'surface freshwater flux (E-P)';
nw{'swflux'}.units = ncchar('centimeter day-1');
nw{'swflux'}.units = 'centimeter day-1';
nw{'swflux'}.field = ncchar('surface freshwater flux, scalar, series');
nw{'swflux'}.field = 'surface freshwater flux, scalar, series';
nw{'swflux'}.positive = ncchar('net evaporation');
nw{'swflux'}.positive = 'net evaporation';
nw{'swflux'}.negative = ncchar('net precipitation');
nw{'swflux'}.negative = 'net precipitation';

nw{'SST'} = ncdouble('sst_time', 'eta_rho', 'xi_rho');
nw{'SST'}.long_name = ncchar('sea surface temperature');
nw{'SST'}.long_name = 'sea surface temperature';
nw{'SST'}.units = ncchar('Celsius');
nw{'SST'}.units = 'Celsius';
nw{'SST'}.field = ncchar('sea surface temperature, scalar, series');
nw{'SST'}.field = 'sea surface temperature, scalar, series';

nw{'SSS'} = ncdouble('sss_time', 'eta_rho', 'xi_rho');
nw{'SSS'}.long_name = ncchar('sea surface salinity');
nw{'SSS'}.long_name = 'sea surface salinity';
nw{'SSS'}.units = ncchar('PSU');
nw{'SSS'}.units = 'PSU';
nw{'SSS'}.field = ncchar('sea surface salinity, scalar, series');
nw{'SSS'}.field = 'sea surface salinity, scalar, series';

nw{'dQdSST'} = ncdouble('sst_time', 'eta_rho', 'xi_rho');
nw{'dQdSST'}.long_name = ncchar('surface net heat flux sensitivity to SST');
nw{'dQdSST'}.long_name = 'surface net heat flux sensitivity to SST';
nw{'dQdSST'}.units = ncchar('Watts meter-2 Celsius-1');
nw{'dQdSST'}.units = 'Watts meter-2 Celsius-1';
nw{'dQdSST'}.field = ncchar('dQdSST, scalar, series');
nw{'dQdSST'}.field = 'dQdSST, scalar, series';

nw{'swrad'} = ncdouble('srf_time', 'eta_rho', 'xi_rho');
nw{'swrad'}.long_name = ncchar('solar shortwave radiation');
nw{'swrad'}.long_name = 'solar shortwave radiation';
nw{'swrad'}.units = ncchar('Watts meter-2');
nw{'swrad'}.units = 'Watts meter-2';
nw{'swrad'}.field = ncchar('shortwave radiation, scalar, series');
nw{'swrad'}.field = 'shortwave radiation, scalar, series';
nw{'swrad'}.positive = ncchar('downward flux, heating');
nw{'swrad'}.positive = 'downward flux, heating';
nw{'swrad'}.negative = ncchar('upward flux, cooling');
nw{'swrad'}.negative = 'upward flux, cooling';

nw{'tide_period'} = ncdouble('tide_period');
nw{'tide_period'}.long_name = ncchar('Tide angular period');
nw{'tide_period'}.long_name = 'Tide angular period';
nw{'tide_period'}.units = ncchar('Hours');
nw{'tide_period'}.units = 'Hours';

nw{'tide_Ephase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Ephase'}.long_name = ncchar('Tidal elevation phase angle');
nw{'tide_Ephase'}.long_name = 'Tidal elevation phase angle';
nw{'tide_Ephase'}.units = ncchar('Degrees');
nw{'tide_Ephase'}.units = 'Degrees';

nw{'tide_Eamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Eamp'}.long_name = ncchar('Tidal elevation amplitude');
nw{'tide_Eamp'}.long_name = 'Tidal elevation amplitude';
nw{'tide_Eamp'}.units = ncchar('Meter');
nw{'tide_Eamp'}.units = 'Meter';

nw{'tide_Cmin'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cmin'}.long_name = ncchar('Tidal current ellipse semi-minor axis');
nw{'tide_Cmin'}.long_name = 'Tidal current ellipse semi-minor axis';
nw{'tide_Cmin'}.units = ncchar('Meter second-1');
nw{'tide_Cmin'}.units = 'Meter second-1';

nw{'tide_Cmax'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cmax'}.long_name = ncchar('Tidal current, ellipse semi-major axis');
nw{'tide_Cmax'}.long_name = 'Tidal current, ellipse semi-major axis';
nw{'tide_Cmax'}.units = ncchar('Meter second-1');
nw{'tide_Cmax'}.units = 'Meter second-1';

nw{'tide_Cangle'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cangle'}.long_name = ncchar('Tidal current inclination angle');
nw{'tide_Cangle'}.long_name = 'Tidal current inclination angle';
nw{'tide_Cangle'}.units = ncchar('Degrees between semi-major axis and East');
nw{'tide_Cangle'}.units = 'Degrees between semi-major axis and East';

nw{'tide_Cphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cphase'}.long_name = ncchar('Tidal current phase angle');
nw{'tide_Cphase'}.long_name = 'Tidal current phase angle';
nw{'tide_Cphase'}.units = ncchar('Degrees');
nw{'tide_Cphase'}.units = 'Degrees';

nw{'tide_Pamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Pamp'}.long_name = ncchar('Tidal potential amplitude');
nw{'tide_Pamp'}.long_name = 'Tidal potential amplitude';
nw{'tide_Pamp'}.units = ncchar('Meter');
nw{'tide_Pamp'}.units = 'Meter';

nw{'tide_Pphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Pphase'}.long_name = ncchar('Tidal potential phase angle');
nw{'tide_Pphase'}.long_name = 'Tidal potential phase angle';
nw{'tide_Pphase'}.units = ncchar('Degrees');
nw{'tide_Pphase'}.units = 'Degrees';

endef(nw);

%
% Create global attributes
%

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.parent_file = ncchar(parentname);
nw.parent_file = parentname;

%
% Write time variables
%

nw{'sms_time'}(:) = smst;
nw{'shf_time'}(:) = shft;
nw{'swf_time'}(:) = swft;
nw{'sst_time'}(:) = sstt;
nw{'srf_time'}(:) = srft;
nw{'sss_time'}(:) = ssst;

nw{'tide_period'}(:) = tide_period;

close(nw);
return
