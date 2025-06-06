%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Create and fill CROCO clim and bry files with OGCM
%  MERCATOR data.
%
%  On crocotools_param.m, available datasets:
%    SODA     -> Description :  1/4 deg SODA v2.2.4 monthly means
%                                                           (187101 -> 201012)
%    mercator -> 2 options for mercator_type parameter:
%
%                1 -->  1/12 deg Mercator global reanalysis (1993   -> 202101)
%                                           (Monthly or Daily product)
%                2 -->  1/12 deg Mercator global analysis   (202011 -> )
%
%            Note for 2:
%                    (7days from now -> now-1day : Near Real Time analysis)
%                    (more than 15days from now  : Best analysis)
%
%  Online reference to MERCATOR is at http://marine.copernicus.eu
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
%  Contributions of P. Marchesiello (IRD), J. Lefevre (IRD),
%                   and F. Colberg (UCT)
%
%  Updated    6-Sep-2006 by Pierrick Penven
%  Update     13-Sep-2009 by Gildas Cambon (IRD)
%  Update     14-March-2011 by Gildas Cambon & Serena Illig (IRD)
%  Update     24-November-2020 by Gildas Cambon(IRD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start % to be used in batch mode %
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param
%
disp('=============================================')
disp('Take care of OGCM defined in crocotools_param')
disp('=============================================')
%
%  MERCATOR : see get_file_python_mercator.m for python-motu, url, and others options
%
itolap_tot=itolap_a + itolap_p;
disp(['Overlap before =',num2str(itolap_a)])
disp(['Overlap after =',num2str(itolap_p)])
disp(['Total overlap =',num2str(itolap_tot)])
disp(['...'])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(['==================='])
disp([' Compute time axis '])
disp(['==================='])
if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)
%
% Extract data over the internet
%
if Download_data==1
    %
    % Get the model limits
    %
    lonmin=min(min(lon));
    lonmax=max(max(lon));
    latmin=min(min(lat));
    latmax=max(max(lat));
    %
    % Download data with DODS (the download matlab routine depends on the OGCM)
    %
    disp('Download data...')
    % Loop on the years
    %
    for Y=Ymin:Ymax
        disp(['Processing year: ',num2str(Y)])
        if Y==Ymin
            mo_min=Mmin;
        else
            mo_min=1;
        end
        if Y==Ymax
            mo_max=Mmax;
        else
            mo_max=12;
        end
        for M=mo_min:mo_max
            thedatemonth=['Y',num2str(Y),'M',num2str(M)];
            download_mercator_python(pathCMC,user,password,mercator_type,...
                                     product_id, ...
                                     Y,M, ...
                                     lonmin,lonmax,latmin,latmax,hmax, ...
                                     OGCM_dir,OGCM_prefix,thedatemonth,Yorig)
        end
    end  %End loop over month and years for the downloading with python/motu client
end %End loop for Download_data
%
%------------------------------------------------------------------------------------
%
% Get the OGCM grid
%
nc=nc4.netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),'.cdf']);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=-nc{'depth'}(:);
NZ=length(Z);
NZ=NZ-rmdepth;
Z=Z(1:NZ);
close(nc)
%
% Initial file
% (the strategy is to start at the begining of a month)
% it is possible to do some temporal interpolation...
% but I am too lazy. lets start the first day of
% month Mmin of year Ymin... with the first data available.
%
if makeini==1
    if  ~exist('vtransform')
        vtransform=1; %Old Vtransform
        disp([' NO VTRANSFORM parameter found'])
        disp([' USE vtransform default value  Vtransfor = 1'])
    end
  ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,Mmin)),nc_suffix];
  %
  % Process the time in Yorig time (i.e days since Yorig-01-01)
  %
  tini=datenum(Ymin,Mmin,1)-datenum(Yorig,1,1);
  disp([' '])
  disp(['=> Create an initial file for ',datestr(tini+datenum(Yorig,1,1));])
  create_inifile(ininame,grdname,CROCO_title,...
		 theta_s,theta_b,hc,N,...
		 tini,'clobber', vtransform);
  nc_ini=netcdf(ininame,'write');
  interp_OGCM2(OGCM_dir,OGCM_prefix,Ymin,Mmin,Roa,interp_method,...
	            lonU,latU,lonV,latV,lonT,latT,Z,1,...
	            nc_ini,[],lon,lat,angle,h,1,obc,vtransform)
  nc_add_globatt(ininame,Yorig,Mmin,Dmin,Hmin,Min_min,Smin,OGCM);
  close(nc_ini)
end
%
% Clim and Bry files
%
if makeclim==1 | makebry==1
    if  ~exist('vtransform')
        vtransform=1; %Old Vtransform
        disp([' NO VTRANSFORM parameter found'])
        disp([' USE vtransform default value  Vtransfor = 1'])
    end
  %
  % Loop on the years and the months
  %
  for Y=Ymin:Ymax
    if Y==Ymin
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
      mo_max=12;
    end
    for M=mo_min:mo_max
    tic
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
	    ' - month ',num2str(M)])
      disp(' ')
      %
      Mm=M-1;Ym=Y;
      if Mm==0
	      Mm=12;
	      Ym=Y-1;
      end
      Mp=M+1;Yp=Y;
      if Mp==13
	      Mp=1;
	      Yp=Y+1;
      end
      %
      % Add 2 times step in the CROCO files: 1 at the beginning and 1 at the end
      %
      nc=netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(Y),'M',num2str(M),'.cdf']);
      OGCM_time=nc{'time'}(:);
      ntimes=length(OGCM_time);
      if ntimes==1
	dt=30; % monthly files
        itolap_a=1; itolap_p=1;
        itolap_tot=itolap_a + itolap_p;
        disp(['Reduced overlap for monthly files'])
        disp(['...'])
      else
	      dt=max(gradient(OGCM_time));
      end
      %
      %% Fill the time axis
      %
      croco_time=0*(1:ntimes+itolap_tot);
      %Current month
      croco_time(itolap_a+1:end-itolap_p)=OGCM_time;
      %
      %Previous  month
      %
      disp(['==================================='])
      for aa= 1:itolap_a
	      disp(['Compute beginning overlap, time index:',num2str(aa)])
	      disp(['Add ',num2str(-(itolap_a + 1 - aa)), ' timestep dt'])
	      disp(['--------'])
	      croco_time(aa) = croco_time(itolap_a+1) - ((itolap_a + 1 - aa).* dt);
      end
      %
      %Next month
      %
      disp(['==================================='])
      for aa= 1:itolap_p
	      disp(['Compute end overlap, time index:',num2str(ntimes+itolap_tot - itolap_p + aa)])
	      disp(['Add ',num2str(aa), ' timestep dt'])
        disp(['--------'])
	      croco_time(end - itolap_p +  aa   ) = croco_time(end - itolap_p) +  aa.* dt;
      end
      disp(['==================================='])
      close(nc)
      %-----------------------------------------------------
      %
      % Create and open the CROCO files
      %
      if makebry==1
	      bryname=[bry_prefix,'Y',num2str(Y),...
		            'M',num2str(sprintf(Mth_format,M)),nc_suffix];
	      disp([' '])
          % only create when file does not exist
          if ~exist(bryname, 'file')
              disp(['=> Create bry file : ',  bryname ])
	          create_bryfile(bryname,grdname,CROCO_title,[1 1 1 1],...
		           theta_s,theta_b,hc,N,...
		           croco_time,0,'clobber',vtransform);
          end
          nc_add_globatt(bryname,Yorig,Mmin,Dmin,Hmin,Min_min,Smin,OGCM);
	      nc_bry=netcdf(bryname,'write');
          % XXX: cara akses attribute
          % nc_bry.hasInterpPrev = 2;
          % disp(nc_bry.hasInterpPrev(1))
      else
	      nc_bry=[];
      end
      if makeclim==1
	      clmname=[clm_prefix,'Y',num2str(Y),...
	            	 'M',num2str(sprintf(Mth_format,M)),nc_suffix];
          % only create when file does not exist
          if ~exist(clmname, 'file')
	          create_climfile(clmname,grdname,CROCO_title,...
	                  		    theta_s,theta_b,hc,N,...
			  croco_time,0,'clobber',vtransform);
          end
          nc_add_globatt(clmname,Yorig,Mmin,Dmin,Hmin,Min_min,Smin,OGCM);
          nc_clm=netcdf(clmname,'write');
      else
	      nc_clm=[];
      end
      %
      % Check if there are OGCM files for the previous Month
      %
      fname=[OGCM_dir,OGCM_prefix,'Y',num2str(Ym),'M',num2str(Mm),'.cdf'];
      if exist(fname, 'file')==0
	      disp(['   No data for the previous month: using current month'])
        Mm=M;
	      Ym=Y;
	      tndx_OGCM=ones(itolap_a,1);
      else
	      nc=netcdf(fname);
	      tndx_OGCM=[(length(nc('time'))- (itolap_a -1) ):1: (length(nc('time')))];
	      close(nc)
      end
      %
      % Perform the interpolations for the previous month
      %
      disp(' Previous month :')
      disp('=================')

      if ( (~isempty(nc_clm) & nc_clm.hasInterpPrev(1)==0) | (~isempty(nc_bry) & nc_bry.hasInterpPrev(1)==0) )
          for aa=1:itolap_a
	          disp(['Beg overlap # ', num2str(aa),' -> tindex ',num2str(aa)])
	          disp(['It. of prev month used for it= ',num2str(tndx_OGCM(aa))])
	          interp_OGCM2(OGCM_dir,OGCM_prefix,Ym,Mm,Roa,interp_method,...
	            lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM(aa),...
		        nc_clm,nc_bry,lon,lat,angle,h,aa,obc,vtransform)
          end
          if ~isempty(nc_bry), nc_bry.hasInterpPrev=1; end
          if ~isempty(nc_clm), nc_clm.hasInterpPrev=1; end
      else
          disp(" Interpolation for previous month already done, skipped")
      end

      %
      % Perform the interpolations for the current month
      %
      disp(' Current month :')
      disp('================')
      if ( (~isempty(nc_clm) & nc_clm.hasInterpCur(1)==0) | (~isempty(nc_bry) & nc_bry.hasInterpCur(1)==0) )
          for tndx_OGCM=1:ntimes
	          disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(ntimes),' :'])
	          interp_OGCM2(OGCM_dir,OGCM_prefix,Y,M,Roa,interp_method,...
                lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM,...
		        nc_clm,nc_bry,lon,lat,angle,h,tndx_OGCM+itolap_a,obc,vtransform)
          end
          if ~isempty(nc_bry), nc_bry.hasInterpCur=1; end
          if ~isempty(nc_clm), nc_clm.hasInterpCur=1; end
      else
          disp(" Interpolation for current month already done, skipped")
      end

      %
      % Read the OGCM file for the next month
      %
      fname=[OGCM_dir,OGCM_prefix,'Y',num2str(Yp),'M',num2str(Mp),'.cdf'];
      if exist(fname,'file')==0
	      disp(['   No data for the next month: using current month'])
	      Mp=M;
	      Yp=Y;
	      for aa=1:itolap_p
	        tndx_OGCM(aa)=ntimes;
	      end
      else
	    for aa=1:itolap_p
	      tndx_OGCM(aa)=aa;
	    end;
      end
      %
      % Perform the interpolations for the next month
      %
      disp(' Next month :')
      disp('=============')
      if ( (~isempty(nc_clm) & nc_clm.hasInterpNext(1)==0) | (~isempty(nc_bry) & nc_bry.hasInterpNext(1)==0) )
          for aa=1:itolap_p
	          disp(['End Overlap #',num2str(aa),' -> tindex ',num2str(ntimes+itolap_a+aa)])
	          disp(['It. of next month used for it= ',num2str(tndx_OGCM(aa))])
	          interp_OGCM2(OGCM_dir,OGCM_prefix,Yp,Mp,Roa,interp_method,...
	            	        lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM(aa),...
		        nc_clm,nc_bry,lon,lat,angle,h,ntimes+itolap_a+aa,obc,vtransform)
          end
          if ~isempty(nc_bry), nc_bry.hasInterpNext=1; end
          if ~isempty(nc_clm), nc_clm.hasInterpNext=1; end
      else
          disp(" Interpolation for next month already done, skipped")
      end
      %
      % Close the CROCO files
      %
      if ~isempty(nc_clm)
	      close(nc_clm);
      end
      if ~isempty(nc_bry)
	      close(nc_bry);
      end
      %
    end
    tx=toc;
    disp(['Elapsed comptime for month: ', num2str(tx), ' seconds'])
  end
end
%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%
if SPIN_Long>0
  %
  % Initial file
  %
  if makeini==1
    %
    % Copy the file
    %
    ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,Mmin)),nc_suffix];
    ininame2=[ini_prefix,'Y',num2str(Ymin-SPIN_Long),'M',num2str(sprintf(Mth_format,Mmin)),nc_suffix];
    disp(['Create ',ininame2])
    eval(['!cp ',ininame,' ',ininame2])
    %
    % Change the time
    %
    nc=netcdf(ininame2,'write');
    time=nc{'scrum_time'}(:)-365.*SPIN_Long*(24*3600);
    nc{'scrum_time'}(:)=time;
    close(nc)
  end
  %
  M=Mmin-1;
  Y=Ymin-SPIN_Long;
  for month=1:12*SPIN_Long
    M=M+1;
    if M==13
      M=1;
      Y=Y+1;
    end
    %
    % Climatology files
    %
    if makeclim==1
      %
      % Copy the file
      %
      clmname=[clm_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
      clmname2=[clm_prefix,'Y',num2str(Y),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
      disp(['Create ',clmname2])
      eval(['!cp ',clmname,' ',clmname2])
      %
      % Change the time
      %
      nc=netcdf(clmname2,'write');
      time=nc{'tclm_time'}(:)-365.*(Ymin-Y);
      nc{'tclm_time'}(:)=time;
      nc{'temp_time'}(:)=time;
      nc{'sclm_time'}(:)=time;
      nc{'salt_time'}(:)=time;
      nc{'uclm_time'}(:)=time;
      nc{'vclm_time'}(:)=time;
      nc{'v2d_time'}(:)=time;
      nc{'v3d_time'}(:)=time;
      nc{'ssh_time'}(:)=time;
      nc{'zeta_time'}(:)=time;
      close(nc)
    end
    %
    % Boundary files
    %
    if makebry==1
      %
      % Copy the file
      %
      bryname=[bry_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
      bryname2=[bry_prefix,'Y',num2str(Y),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
      disp(['Create ',bryname2])
      eval(['!cp ',bryname,' ',bryname2])
      %
      % Change the time
      %
      nc=netcdf(bryname2,'write');
      time=nc{'bry_time'}(:)-365.*(Ymin-Y);
      nc{'bry_time'}(:)=time;
      close(nc)
    end
  end
end
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp(' Make a few plots...')
  if makeini==1
    ininame=[ini_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,Mmin)),nc_suffix];
    figure
    test_clim(ininame,grdname,'temp',1,coastfileplot)
    figure
    test_clim(ininame,grdname,'salt',1,coastfileplot)
  end
  if makeclim==1
    clmname=[clm_prefix,'Y',num2str(Y),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
    figure
    test_clim(clmname,grdname,'temp',1,coastfileplot)
    figure
    test_clim(clmname,grdname,'salt',1,coastfileplot)
  end
  if makebry==1
    bryname=[bry_prefix,'Y',num2str(Y),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
    figure
    test_bry(bryname,grdname,'temp',1,obc)
    figure
    test_bry(bryname,grdname,'salt',1,obc)
    figure
    test_bry(bryname,grdname,'u',1,obc)
    figure
    test_bry(bryname,grdname,'v',1,obc)
  end
end
