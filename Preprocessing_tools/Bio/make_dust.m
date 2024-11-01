%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a PISCES forcing file
%
%  Extrapole and interpole surface data to get surface boundary
%  conditions for PISCES (forcing netcdf file)
%
%  Data input format (netcdf):
%     irondep(T, Y, X)
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
% Dust deposition file 
%
%dust_file=[woapisces_dir,'dust.iron.cdf'];
%dust_name='irondep';
%dust_file=[woapisces_dir,'dust_seas.cdf'];
dust_file=[woapisces_dir,'DUST_INCA_new_r360x180.nc'];
dust_name='dust';
%time=woa_time;
time=[0.5:1:11.5];
cycle=woa_cycle;
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
nc = netcdf(bioname, 'clobber');
%nc = netcdf(bioname, 'write');
%%result = redef(nc);
%
nc('xi_u') = L;
nc('xi_v') = Lp;
nc('xi_rho') = Lp;
nc('eta_u') = Mp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
%
nc('dust_time') = length(time);
nc{'dust_time'} = ncdouble('dust_time') ;
nc{'dust'}      = ncdouble('dust_time','eta_rho','xi_rho') ;
nc{'dustfer'}   = ncdouble('dust_time','eta_rho','xi_rho') ;
nc{'dustpo4'}   = ncdouble('dust_time','eta_rho','xi_rho') ;
nc{'dustsi'}    = ncdouble('dust_time','eta_rho','xi_rho') ;
nc{'solubility2'}    = ncdouble('dust_time','eta_rho','xi_rho') ;
%
nc{'dust_time'}.long_name = ncchar('time for dust');
nc{'dust_time'}.long_name = 'time for dust';
nc{'dust_time'}.units = ncchar('day');
nc{'dust_time'}.units = 'day';
if cycle~=0
  nc{'dust_time'}.cycle_length = cycle;
end
%
nc{'dust'}.long_name = ncchar('Dust Deposition');
nc{'dust'}.long_name = 'Dust Deposition';
nc{'dust'}.units = ncchar('Kg m-2 s-1');
nc{'dust'}.units = 'Kg m-2 s-1';
nc{'dust'}.fields = ncchar('dust, scalar, series');
nc{'dust'}.fields = 'dust, scalar, series';
%
nc{'dustfer'}.long_name = ncchar('Fe Dust Deposition');
nc{'dustfer'}.long_name = 'Fe Dust Deposition';
nc{'dustfer'}.units = ncchar('Kg m-2 s-1');
nc{'dustfer'}.units = 'Kg m-2 s-1';
nc{'dustfer'}.fields = ncchar('dustfer, scalar, series');
nc{'dustfer'}.fields = 'dustfer, scalar, series';
%
nc{'dustpo4'}.long_name = ncchar('PO4 Dust Deposition');
nc{'dustpo4'}.long_name = 'PO4 Dust Deposition';
nc{'dustpo4'}.units = ncchar('Kg m-2 s-1');
nc{'dustpo4'}.units = 'Kg m-2 s-1';
nc{'dustpo4'}.fields = ncchar('dustpo4, scalar, series');
nc{'dustpo4'}.fields = 'dustpo4, scalar, series';
%
nc{'dustsi'}.long_name = ncchar('Si Dust Deposition');
nc{'dustsi'}.long_name = 'Si Dust Deposition';
nc{'dustsi'}.units = ncchar('Kg m-2 s-1');
nc{'dustsi'}.units = 'Kg m-2 s-1';
nc{'dustsi'}.fields = ncchar('dustsi, scalar, series');
nc{'dustsi'}.fields = 'dustsi, scalar, series';
%
nc{'solubility2'}.long_name = ncchar('Fe solubility from Mahowald');
nc{'solubility2'}.long_name = 'Fe solubility from Mahowald';
nc{'solubility2'}.units = ncchar('%');
nc{'solubility2'}.units = '%';
nc{'solubility2'}.fields = ncchar('solubility2, scalar, series');
nc{'solubility2'}.fields = 'solubility2, scalar, series';
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
nc{'dust_time'}(:)=time.*30; % if time in month in the dataset !!!

close(nc)

nc=netcdf(bioname,'write');
%
% Loop on time
%
for tindex=1:length(time)
  time=nc{'dust_time'}(tindex);
  nc{'dust'}(tindex,:,:)=ext_data(dust_file,'dust',tindex,...
             lon,lat,time,Roa,1);
  nc{'dustfer'}(tindex,:,:)=ext_data(dust_file,'dustfer',tindex,...
             lon,lat,time,Roa,1);
  nc{'dustpo4'}(tindex,:,:)=ext_data(dust_file,'dustpo4',tindex,...
             lon,lat,time,Roa,1);
  nc{'dustsi'}(tindex,:,:)=ext_data(dust_file,'dustsi',tindex,...
             lon,lat,time,Roa,1);
  nc{'solubility2'}(tindex,:,:)=ext_data(dust_file,'solubility2',tindex,...
             lon,lat,time,Roa,1);
end
close(nc)
if (makeplot)
%
% Make a few plots
%
disp(' Make a few plots...')
test_bioforcing(bioname,grdname,'dust',[1 4 7 10],3,coastfileplot)
test_bioforcing(bioname,grdname,'dustfer',[1 4 7 10],3,coastfileplot)
test_bioforcing(bioname,grdname,'dustpo4',[1 4 7 10],3,coastfileplot)
test_bioforcing(bioname,grdname,'dustsi',[1 4 7 10],3,coastfileplot)
end % if makeplot
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
