%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a PISCES forcing file
%
%  Extrapole and interpole surface data to get surface boundary
%  conditions for PISCES (forcing netcdf file)
%
%  Data input format (netcdf):
%     ndep(T, Y, X)
%     T : time [Months]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO Climate Data Library 
%                (Atlas of Surface Marine Data 1994)
%
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.DASILVA/
%
%  Pierrick Penven, IRD, 2005.                                    %
%  Olivier Aumont the master, IRD, 2006.                          %
%  Patricio Marchesiello, chief, IRD, 2007.                       %
%  Christophe Eugene Raoul Menkes, the slave, IRD, 2007.          %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('Creating biology forcing file')
%
%  Title - Grid file name - Forcing file name
%
crocotools_param
%
% bioname
%
%bioname='croco_frcbio.nc'
%
% Nitrogen deposition file 
%
ndep_file=[woapisces_dir,'Ndep_CMIP_NCAR-CCMI-2-0_gn_199001-201012-clim.nc'];
ndep_name='ndep';
%time=woa_time;
time=[0.5:1:11.5];
cycle=woa_cycle;
disp(' ')
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname,'r');
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
L=Lp-1;
M=Mp-1;
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
close(nc);
%
% create dust forcing file
%
disp(' Creating file')
%nc = netcdf(bioname, 'clobber');
nc = netcdf(bioname, 'write');
%%result = redef(nc);
%
nc('xi_u') = L;
nc('xi_v') = Lp;
nc('xi_rho') = Lp;
nc('eta_u') = Mp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
%
nc('ndep_time') = length(time);
nc{'ndep_time'} = ncdouble('ndep_time') ;
nc{'ndep'} = ncdouble('ndep_time','eta_rho','xi_rho') ;
nc{'noyndepo'} = ncdouble('ndep_time','eta_rho','xi_rho') ;
nc{'nhxndepo'} = ncdouble('ndep_time','eta_rho','xi_rho') ;
%
nc{'ndep_time'}.long_name = ncchar('time for nitrogen deposition');
nc{'ndep_time'}.long_name = 'time for nitrogen deposition';
nc{'ndep_time'}.units = ncchar('day');
nc{'ndep_time'}.units = 'day';
if cycle~=0
  nc{'ndep_time'}.cycle_length = cycle;
end
%
nc{'ndep'}.long_name = ncchar('Nitrogen Deposition');
nc{'ndep'}.long_name = 'Nitrogen Deposition';
nc{'ndep'}.units = ncchar('KgN m-2 s-1');
nc{'ndep'}.units = 'KgN m-2 s-1';
nc{'ndep'}.fields = ncchar('ndep, scalar, series');
nc{'ndep'}.fields = 'ndep, scalar, series';
%
nc{'noyndepo'}.long_name = ncchar('NOy Deposition');
nc{'noyndepo'}.long_name = 'NOy Deposition';
nc{'noyndepo'}.units = ncchar('KgN m-2 s-1');
nc{'noyndepo'}.units = 'KgN m-2 s-1';
nc{'noyndepo'}.fields = ncchar('noyndepo, scalar, series');
nc{'noyndepo'}.fields = 'noyndepo, scalar, series';
%
nc{'nhxndepo'}.long_name = ncchar('NHx Deposition');
nc{'nhxndepo'}.long_name = 'NHx Deposition';
nc{'nhxndepo'}.units = ncchar('KgN m-2 s-1');
nc{'nhxndepo'}.units = 'KgN m-2 s-1';
nc{'nhxndepo'}.fields = ncchar('nhxndepo, scalar, series');
nc{'nhxndepo'}.fields = 'nhxndepo, scalar, series';
%
%%endef(nc);

% Create global attributes
nc.title = ncchar(CROCO_title);
nc.title = CROCO_title;
nc.date = ncchar(date);
nc.date = date;
nc.grd_file = ncchar(bioname);
nc.grd_file = grdname;
nc.type = ncchar('CROCO biology forcing file');
nc.type = 'CROCO biology forcing file';

% Write time variable
nc{'ndep_time'}(:)=time.*30; % if time in month in the dataset !!!

close(nc)

nc=netcdf(bioname,'write');
%
% Loop on time
for tindex=1:length(time)
  time=nc{'ndep_time'}(tindex);
  nc{'ndep'}(tindex,:,:)=ext_data(ndep_file,'ndep',tindex,...
             lon,lat,time,Roa,1);
  nc{'noyndepo'}(tindex,:,:)=ext_data(ndep_file,'noyndepo',tindex,...
             lon,lat,time,Roa,1);
  nc{'nhxndepo'}(tindex,:,:)=ext_data(ndep_file,'nhxndepo',tindex,...
             lon,lat,time,Roa,1);
end
%
close(nc)
if (makeplot)
%
% Make a few plots
%
disp(' Make a few plots...')
test_bioforcing(bioname,grdname,'ndep',[1 4 7 10],3,coastfileplot)
test_bioforcing(bioname,grdname,'noyndepo',[1 4 7 10],3,coastfileplot)
test_bioforcing(bioname,grdname,'nhxndepo',[1 4 7 10],3,coastfileplot)
end % if makeplot
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
