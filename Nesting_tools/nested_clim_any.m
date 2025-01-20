function nested_clim(parent_grd,child_grd,parent_clim,child_clim,...
                     vertical_correc,extrapmask)

%
% Title
%
title=['Climatology file for the new grid :',child_clim,...
       ' using parent forcing file: ',parent_clim];
disp(' ')
disp(title)
if extrapmask==1
  disp('Extrapolation under mask is on')
  disp('====================')
end
%
if vertical_correc==1
  disp('Vertical correction is on')
  disp('===============')
end

parent_grid_nc = netcdf(parent_grd);
child_grid_nc  = netcdf(child_grd);
parent_clim_nc = netcdf(parent_blk);

Lp=length(parent_grid_nc('xi_rho'));
Mp=length(parent_grid_nc('eta_rho'));
Np=length(parent_grid_nc('s_rho'));
Lpc=length(child_grid_nc('xi_rho'));
Mpc=length(child_grid_nc('eta_rho'));
Npc=length(child_grid_nc('s_rho'));
if extrapmask==1
  mask=parent_grid_nc{'mask_rho'}(:);
else
  mask=[];
end
%
% Read in the parent climatology file
%
disp(' ')
disp(' Read in the parent climatology file...')
%nc = netcdf(parent_clim);
theta_s = parent_clim_nc{'theta_s'}(:);
theta_b = parent_clim_nc{'theta_b'}(:);
Tcline = parent_clim_nc{'Tcline'}(:);
N = length(parent_clim_nc('s_rho'));
vtransform=parent_clim_nc{'Vtransform'}(:);
hc = parent_clim_nc{'hc'}(:);
disp([' Use parent VTRANSFORM = ',num2str(vtransform)])
if ~exist('vtransform') | isempty(vtransform)
    disp([' No VTRANSFORM parameter found'])
    disp([' Use the default one VTRANSFORM = 1'])
    vtransform=1;
end
ttime = parent_clim_nc{'tclm_time'}(:);
tcycle = parent_clim_nc{'tclm_time'}.cycle_length(:);
stime = parent_clim_nc{'sclm_time'}(:);
scycle = parent_clim_nc{'sclm_time'}.cycle_length(:);
utime = parent_clim_nc{'ssh_time'}(:);
ucycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
vtime = parent_clim_nc{'ssh_time'}(:);
vcycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
sshtime = parent_clim_nc{'ssh_time'}(:);
sshcycle = parent_clim_nc{'ssh_time'}.cycle_length(:);
climtime=ttime;
if stime~=climtime | utime~=climtime | vtime~=climtime | sshtime~=climtime
  error('Nested_clim: different times... ')
end

%
% Create the climatology file
%
disp(' ')
disp(' Create the climatology file...')
type = 'CLIMATOLOGY file' ;
history = 'CROCO' ;
child_clim_nc = netcdf(climfile,clobber);
redef(child_clim_nc);
child_clim_nc('xi_u') = Lpc-1;
child_clim_nc('xi_v') = Lpc;
child_clim_nc('xi_rho') = Lpc;
child_clim_nc('eta_u') = Mpc;
child_clim_nc('eta_v') = Mpc-1;
child_clim_nc('eta_rho') = Mpc;
child_clim_nc('s_rho') = Npc;
child_clim_nc('s_w') = Npc+1;
child_clim_nc('tracer') = 2;
child_clim_nc('tclm_time') = length(ttime);
child_clim_nc('sclm_time') = length(stime);
child_clim_nc('uclm_time') = length(utime);
child_clim_nc('vclm_time') = length(vtime);
child_clim_nc('ssh_time') = length(sshtime);
child_clim_nc('one') = 1;
%
%  Create variables
%
child_clim_nc{'spherical'} = ncchar('one') ;
child_clim_nc{'Vtransform'} = ncint('one') ;
child_clim_nc{'tstart'} = ncdouble('one') ;
child_clim_nc{'tend'} = ncdouble('one') ;
child_clim_nc{'theta_s'} = ncdouble('one') ;
child_clim_nc{'theta_b'} = ncdouble('one') ;
child_clim_nc{'Tcline'} = ncdouble('one') ;
child_clim_nc{'hc'} = ncdouble('one') ;
child_clim_nc{'s_rho'} = ncdouble('s_rho') ;
child_clim_nc{'Cs_rho'} = ncdouble('s_rho') ;
child_clim_nc{'tclm_time'} = ncdouble('tclm_time') ;
child_clim_nc{'sclm_time'} = ncdouble('sclm_time') ;
child_clim_nc{'uclm_time'} = ncdouble('uclm_time') ;
child_clim_nc{'vclm_time'} = ncdouble('vclm_time') ;
child_clim_nc{'ssh_time'} = ncdouble('ssh_time') ;
child_clim_nc{'temp'} = ncdouble('tclm_time','s_rho','eta_rho','xi_rho') ;
child_clim_nc{'salt'} = ncdouble('sclm_time','s_rho','eta_rho','xi_rho') ;
child_clim_nc{'u'} = ncdouble('uclm_time','s_rho','eta_u','xi_u') ;
child_clim_nc{'v'} = ncdouble('vclm_time','s_rho','eta_v','xi_v') ;
child_clim_nc{'ubar'} = ncdouble('uclm_time','eta_u','xi_u') ;
child_clim_nc{'vbar'} = ncdouble('vclm_time','eta_v','xi_v') ;
child_clim_nc{'SSH'} = ncdouble('ssh_time','eta_rho','xi_rho') ;
%
%
%  Create attributes
%
child_clim_nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
child_clim_nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
child_clim_nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
child_clim_nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
child_clim_nc{'tstart'}.long_name = ncchar('start processing day');
child_clim_nc{'tstart'}.long_name = 'start processing day';
child_clim_nc{'tstart'}.units = ncchar('day');
child_clim_nc{'tstart'}.units = 'day';
%
child_clim_nc{'tend'}.long_name = ncchar('end processing day');
child_clim_nc{'tend'}.long_name = 'end processing day';
child_clim_nc{'tend'}.units = ncchar('day');
child_clim_nc{'tend'}.units = 'day';
%
child_clim_nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
child_clim_nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
child_clim_nc{'theta_s'}.units = ncchar('nondimensional');
child_clim_nc{'theta_s'}.units = 'nondimensional';
%
child_clim_nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
child_clim_nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
child_clim_nc{'theta_b'}.units = ncchar('nondimensional');
child_clim_nc{'theta_b'}.units = 'nondimensional';
%
child_clim_nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
child_clim_nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
child_clim_nc{'Tcline'}.units = ncchar('meter');
child_clim_nc{'Tcline'}.units = 'meter';
%
child_clim_nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
child_clim_nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
child_clim_nc{'hc'}.units = ncchar('meter');
child_clim_nc{'hc'}.units = 'meter';
%
child_clim_nc{'s_rho'}.long_name = ncchar('S-coordinate at RHO-points');
child_clim_nc{'s_rho'}.long_name = 'S-coordinate at RHO-points';
child_clim_nc{'s_rho'}.units = ncchar('nondimensional');
child_clim_nc{'s_rho'}.units = 'nondimensional';
child_clim_nc{'s_rho'}.valid_min = -1;
child_clim_nc{'s_rho'}.valid_max = 0;
child_clim_nc{'s_rho'}.field = ncchar('s_rho, scalar');
child_clim_nc{'s_rho'}.field = 's_rho, scalar';
%
child_clim_nc{'Cs_rho'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
child_clim_nc{'Cs_rho'}.long_name = 'S-coordinate stretching curves at RHO-points';
child_clim_nc{'Cs_rho'}.units = ncchar('nondimensional');
child_clim_nc{'Cs_rho'}.units = 'nondimensional';
child_clim_nc{'Cs_rho'}.valid_min = -1;
child_clim_nc{'Cs_rho'}.valid_max = 0;
child_clim_nc{'Cs_rho'}.field = ncchar('Cs_rho, scalar');
child_clim_nc{'Cs_rho'}.field = 'Cs_rho, scalar';
%
child_clim_nc{'tclm_time'}.long_name = ncchar('time for temperature climatology');
child_clim_nc{'tclm_time'}.long_name = 'time for temperature climatology';
child_clim_nc{'tclm_time'}.units = ncchar('day');
child_clim_nc{'tclm_time'}.units = 'day';
child_clim_nc{'tclm_time'}.cycle_length = tcycle;
child_clim_nc{'tclm_time'}.field = ncchar('tclm_time, scalar, series');
child_clim_nc{'tclm_time'}.field = 'tclm_time, scalar, series'  ;
%
child_clim_nc{'sclm_time'}.long_name = ncchar('time for salinity climatology');
child_clim_nc{'sclm_time'}.long_name = 'time for salinity climatology';
child_clim_nc{'sclm_time'}.units = ncchar('day');
child_clim_nc{'sclm_time'}.units = 'day';
child_clim_nc{'sclm_time'}.cycle_length = scycle;
child_clim_nc{'sclm_time'}.field = ncchar('sclm_time, scalar, serie');
child_clim_nc{'sclm_time'}.field = 'sclm_time, scalar, serie';
%
child_clim_nc{'uclm_time'}.long_name = ncchar('time climatological u');
child_clim_nc{'uclm_time'}.long_name = 'time climatological u';
child_clim_nc{'uclm_time'}.units = ncchar('day');
child_clim_nc{'uclm_time'}.units = 'day';
child_clim_nc{'uclm_time'}.cycle_length = ucycle;
child_clim_nc{'uclm_time'}.field = ncchar('uclm_time, scalar, serie');
child_clim_nc{'uclm_time'}.field = 'uclm_time, scalar, serie';
%
child_clim_nc{'vclm_time'}.long_name = ncchar('time climatological v');
child_clim_nc{'vclm_time'}.long_name = 'time climatological v';
child_clim_nc{'vclm_time'}.units = ncchar('day');
child_clim_nc{'vclm_time'}.units = 'day';
child_clim_nc{'vclm_time'}.cycle_length = vcycle;
child_clim_nc{'vclm_time'}.field = ncchar('vclm_time, scalar, serie');
child_clim_nc{'vclm_time'}.field = 'vclm_time, scalar, serie';
%
child_clim_nc{'ssh_time'}.long_name = ncchar('time for sea surface height');
child_clim_nc{'ssh_time'}.long_name = 'time for sea surface height';
child_clim_nc{'ssh_time'}.units = ncchar('day');
child_clim_nc{'ssh_time'}.units = 'day';
child_clim_nc{'ssh_time'}.cycle_length = sshcycle;
child_clim_nc{'ssh_time'}.field = ncchar('ssh_time, scalar, serie');
child_clim_nc{'ssh_time'}.field = 'ssh_time, scalar, serie';
%
% Create global attributes
%
child_clim_nc.title = ncchar(title);
child_clim_nc.title = title;
child_clim_nc.date = ncchar(date);
child_clim_nc.date = date;
child_clim_nc.clim_file = ncchar(climfile);
child_clim_nc.clim_file = climfile;
child_clim_nc.grd_file = ncchar(gridfile);
child_clim_nc.grd_file = gridfile;
child_clim_nc.parent_file = ncchar(parentfile);
child_clim_nc.parent_file = parentfile;
child_clim_nc.type = ncchar(type);
child_clim_nc.type = type;
child_clim_nc.history = ncchar(history);
child_clim_nc.history = history;
%
% Leave define mode
%
endef(child_clim_nc);
%
[s_rho,Cs_rho,s_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%
%
% Write variables
%
child_clim_nc{'spherical'}(:)='T';
child_clim_nc{'Vtransform'}(:)=vtransform;
child_clim_nc{'tstart'}(:) =  min([min(ttime) min(stime) min(utime)]);
child_clim_nc{'tend'}(:) =  max([max(ttime) max(stime) max(utime)]);
child_clim_nc{'theta_s'}(:) =  theta_s;
child_clim_nc{'theta_b'}(:) =  theta_b;
child_clim_nc{'Tcline'}(:) =  Tcline;
child_clim_nc{'hc'}(:) =  hc;
child_clim_nc{'s_rho'}(:) =  s_rho;
child_clim_nc{'Cs_rho'}(:) =  Cs_rho;
child_clim_nc{'sclm_time'}(:) =  stime;
child_clim_nc{'tclm_time'}(:) =  ttime;
child_clim_nc{'uclm_time'}(:) = utime ;
child_clim_nc{'vclm_time'}(:) = vtime ;
child_clim_nc{'ssh_time'}(:) = sshtime;
child_clim_nc{'u'}(:) =  0;
child_clim_nc{'v'}(:) =  0;
child_clim_nc{'ubar'}(:) =  0;
child_clim_nc{'vbar'}(:) =  0;
child_clim_nc{'SSH'}(:) =  0;
child_clim_nc{'temp'}(:) =  0;
child_clim_nc{'salt'}(:) =  0;
%
sync(child_clim_nc);

%
% parent indices
%
[igrd_r,jgrd_r]=meshgrid((1:1:Lp),(1:1:Mp));
[igrd_p,jgrd_p]=meshgrid((1:1:Lp-1),(1:1:Mp-1));
[igrd_u,jgrd_u]=meshgrid((1:1:Lp-1),(1:1:Mp));
[igrd_v,jgrd_v]=meshgrid((1:1:Lp),(1:1:Mp-1));
%
% the children indices
%
[ichildgrd_p,jchildgrd_p]=meshgrid((1:1:Lpc-1),(1:1:Mpc-1));
[ichildgrd_r,jchildgrd_r]=meshgrid((1:1:Lpc),(1:1:Mpc));
[ichildgrd_u,jchildgrd_u]=meshgrid((1:1:Lpc-1),(1:1:Mpc));
[ichildgrd_v,jchildgrd_v]=meshgrid((1:1:Lpc),(1:1:Mpc-1));

%
% interpolations
%
disp(' ')
disp(' Do the interpolations...')



close(child_clim_nc);
