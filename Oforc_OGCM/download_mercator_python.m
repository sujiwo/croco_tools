function download_mercator_python(pathCMC,user,password,mercator_type, ...
                                  product_id, ...
                                  Y,M, ...
                                  lonmin,lonmax,latmin,latmax,zmax, ...
                                  OGCM_dir,OGCM_prefix,thedatemonth,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
%  Updated    6-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Update     23-Oct-2020 by Gildas Cambon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
download_raw_data=1;
convert_raw2crocotools=1; % convert -> crocotools format data
%
% Set variable names according to mercator type data
%
vars = {'zos' ...
        'uo' ...
        'vo' ...
        'thetao' ...
        'so'};
%
disp([' '])
disp(['Get data for Y',num2str(Y),'M',num2str(M)])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])
%
% Create the directory
%
disp(['Making output data directory ',OGCM_dir])
mkdir(OGCM_dir);

%
% Get files and dates
%
% for interanual monthly files
time1=datenum(Y,M,01);
%time2=datenum(Y,M+1,01) - 1;
time2=datenum(Y,M+1,01);
%time2=datenum(Y,M,02);  %for debug
time=cat(2,time1,time2);
tiempo_inicial = time(1);
tiempo_final = time(end);

if (lonmin > 180)
    lonmin = lonmin - 360;
end
if (lonmax > 180)
    lonmax = lonmax - 360;
end
raw_mercator_name=[OGCM_dir,'raw_',OGCM_prefix,thedatemonth,'.nc'];

if download_raw_data
    %
    % Get data
    if mercator_type==1
      get_file_python_mercator(pathCMC, ...
                             product_id, ...
                             vars, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             raw_mercator_name);

    elseif mercator_type==2
      %zos
      get_file_python_mercator(pathCMC,...
                             product_id{1}, ...
                             vars{1}, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             [raw_mercator_name(1:end-3),'_z.nc']);

      %uo/vo
      get_file_python_mercator(pathCMC, ...
                             product_id{2}, ...
                             {vars{2:3}}, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             [raw_mercator_name(1:end-3),'_u.nc']);

      %thetao
      get_file_python_mercator(pathCMC, ...
                             product_id{3}, ...
                             vars{4}, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             [raw_mercator_name(1:end-3),'_t.nc']);
      %so
      get_file_python_mercator(pathCMC, ...
                             product_id{4}, ...
                             vars{5}, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             [raw_mercator_name(1:end-3),'_s.nc']);
    end
%
end  %end download_raw_data
%
if convert_raw2crocotools
    %
    % Convert data format and write in a more CROCOTOOLS
    % compatible input file
    %
    mercator_name=[OGCM_dir,'raw_',OGCM_prefix,thedatemonth,'.nc'];
    if exist(mercator_name)
        disp('Mercator file already exist => overwrite it')
    end
    if mercator_type==1
      write_mercator(OGCM_dir,OGCM_prefix,raw_mercator_name, ...
                     mercator_type,vars,time,thedatemonth,Yorig); % write data
%      eval(['! rm -f ',raw_mercator_name]);
    elseif mercator_type==2
      write_mercator_multi(OGCM_dir,OGCM_prefix,raw_mercator_name, ...
                           mercator_type,vars,time,thedatemonth,Yorig); % write data
%      eval(['! rm -f ',raw_mercator_name(1:end-3),'*.nc']);
    end


end

end

