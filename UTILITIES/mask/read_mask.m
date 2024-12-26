function [spherical,x,y,bath,rmask]=read_mask(Gname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hernan G. Arango %%%%
% Copyright (c) 2001 Rutgers University.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% function [spherical,x,y,bath,rmask]=read_mask(Gname)                 %
%                                                                      %
% This routine reads in domain grid, bathymetry, and Land/Sea mask at  %
% from GRID NetCDF file.                                               %
%                                                                      %
% On Input:                                                            %
%                                                                      %
%    Gname       GRID NetCDF file name (character string).             %
%                                                                      %
% On Output:                                                           %
%                                                                      %
%    spherical   grid type switch (logical):                           %
%                spherical=1, spherical grid set-up.                   %
%                spherical=0, Cartesian grid set-up.                   %
%    x           X-location of RHO-points (real matrix).               %
%    y           Y-location of RHO-points (real matrix).               %
%    bath        raw bathymetry at RHO-points (real matrix; meters).   %
%    rmask       Land/Sea mask on RHO-points (real matrix):            %
%                rmask=0 land, rmask=1 Sea.                            %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------
% Inquire about spatial dimensions.
%-----------------------------------------------------------------------

sourcenc = netcdf(Gname);
Im = length( dim(sourcenc, 'xi_rho'));
Jm = length( dim(sourcenc, 'eta_rho'));
got.spher=0;  Vname.spher='spherical';
got.h    =0;  Vname.h    ='h';
got.hraw =0;  Vname.hraw ='hraw';
got.rmask=0;  Vname.rmask='mask_rho';
got.rlon =0;  Vname.rlon ='lon_rho';
got.rlat =0;  Vname.rlat ='lat_rho';
got.xr   =0;  Vname.xr   ='x_rho';
got.yr   =0;  Vname.yr   ='y_rho';


%-----------------------------------------------------------------------
% Read in relevant Land/Sea mask variables.
%-----------------------------------------------------------------------

% Spherical switch.
spherical = 0;
if sourcenc{'spherical'}(1)=='T' | sourcenc{'spherical'}(1)=='t'
    spherical=1;
end

% Grid positions at RHO-points.

if (spherical)
    try
        x = sourcenc{Vname.rlon}(:);
        y = sourcenc{Vname.rlat}(:);
    catch
        [y,x] = meshgrid(1:Jm, 1:Im);
    end
else
    try
        x = sourcenc{Vname.xr}(:);
        y = sourcenc{Vname.yr}(:);
    catch
        [y,x] = meshgrid(1:Jm, 1:Im);
    end
end

% Mask on RHO-points.
try
    rmask = sourcenc{Vname.rmask}(:);
catch
    rmask = ones([Im Jm]);
end

% Bathymetry.
try
    bath = sourcenc{Vname.hraw}(1,:);
catch
    try
        bath = sourcenc{Vname.h};
    catch
        bath = zeros([Im Jm]);
    end
end

close(sourcenc);
return
