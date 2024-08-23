function download_mercator_frcst_python(pathCMC,user,password,mercator_type, ...
                                        raw_mercator_name,product_id, ...
                                        lh,lf, ...
                                        lonmin,lonmax,latmin,latmax,zmax, ...
                                        FRCST_dir,FRCST_prefix,date_frcst,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract a subgrid from mercator to get a CROCO forcing
% Store that into monthly files.
% Take care of the Greenwitch Meridian.
% 
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
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    06-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Update     23-Oct-2020 by Gildas Cambon
%  Update     06-May-2023 by Efrain Rodriguez-Rubio & P. Marchesiello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathMotu is a deprecated parameter ! 
download_raw_data=1;
convert_raw2crocotools=1; % convert -> crocotools format data
%
disp(['Making output data directory ',FRCST_dir]) % create directory
eval(['!mkdir ',FRCST_dir, ' 2> /dev/null'])
%
% Set variable names according to mercator type data
%
if mercator_type==5
  vars = {'zos_detided', ...
          'uo', ...
          'vo', ...
          'thetao', ...
          'so'};
else
  vars = {'zos' ...
          'uo' ...
          'vo' ...
          'thetao' ...
          'so'};
end
%
% Get dates
%
rundate_str=date;
date_run=datenum(rundate_str);
rundate=date_run-datenum(Yorig,1,1);
for i=1:lh+1
    time1(i)=datenum(rundate_str)-(lh+2-i);
end
time2=datenum(rundate_str);
for j=1:lf+2
    time3(j)=datenum(rundate_str)+j-1;
end
time=cat(2,time1,time2,time3);
tiempo_inicial = time1(1);
tiempo_final = time3(end);
tiempo_p0 = time1(lh+1);
tiempo_p1 = date_run;
%
if (lonmin > 180)
    lonmin = lonmin - 360;
end
if (lonmax > 180)
    lonmax = lonmax - 360;
end
%
disp([' '])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude:  ',num2str(latmin)])
disp(['Maximum Latitude:  ',num2str(latmax)])
disp([' '])
%
if download_raw_data
    %
    disp(['--> Extraction Forecasts (assuming local time is after time dissemination : 12hUTC)  :  ', datestr(tiempo_p1,'yyyy-mm-dd')])
    if lh ~= 0
      for i=1:lh
          disp(['Get Hindcasts NRT ana.: ',datestr(time1(i),'yyyy-mm-dd')])
      end
      disp(['Get Hindcasts P0      : ',datestr(tiempo_p0,'yyyy-mm-dd')])
    end
    disp(['Get ',num2str(lf+2),' Forecasts       : from ',datestr(tiempo_p1,'yyyy-mm-dd'),' to ',datestr(tiempo_final,  'yyyy-mm-dd')])
    %
    %
    % Get data  (now splitted in 4 different files)
    get_file_python_mercator(pathCMC, ...                    % SSH
                            product_id{1}, ...
                            vars(1), ...
                            [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                            {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                            datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                            {user password}, ...
                            [raw_mercator_name(1:end-3),'_z.nc']); 
    %
    get_file_python_mercator(pathCMC, ...                    % U/V 
                            product_id{2}, ...
                            vars(2:3), ...
                            [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                            {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                            datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                            {user password}, ...
                            [raw_mercator_name(1:end-3),'_u.nc']); 
    %
    get_file_python_mercator(pathCMC, ...                    % T 
                            product_id{3}, ...
                            vars(4), ...
                            [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                            {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                            datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                            {user password}, ...
                            [raw_mercator_name(1:end-3),'_t.nc']);                         
    %
    get_file_python_mercator(pathCMC, ...                    % S
                            product_id{4}, ...
                            vars(5), ...
                            [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                            {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                            datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                            {user password}, ...
                            [raw_mercator_name(1:end-3),'_s.nc']); 
%
end  % download_raw_data
%
if convert_raw2crocotools 
    %
    % Convert data format and write in a more CROCOTOOLS 
    % compatible input file 
    %
    write_mercator_multi(FRCST_dir,FRCST_prefix,raw_mercator_name, ...
                         mercator_type,vars,time,num2str(rundate),Yorig); % write data
    eval(['! rm -f ',[raw_mercator_name(1:end-3),'*nc '] ]);
end
%
end
