function [imin,imax,jmin,jmax,refinecoeff] = detect_grid(parent_grid_nc, child_grid_nc)
% Create vertical and horizontal that places child on 
% parent's 2Dgrid
    
    parent_lat = sort(parent_grid_nc{'lat_rho'}(:,1))';
    parent_lon = sort(parent_grid_nc{'lon_rho'}(1,:));
    child_lat = sort(child_grid_nc{'lat_rho'}(:,1))';
    child_lon = sort(child_grid_nc{'lon_rho'}(1,:));
    xrho_p = parent_grid_nc{'x_rho'}(1,:);
    yrho_p = parent_grid_nc{'y_rho'}(:,1)';
    xrho_c = child_grid_nc{'x_rho'}(1,:);
    yrho_c = child_grid_nc{'y_rho'}(:,1)';
    scale_x = (xrho_c(2)-xrho_c(1))/(xrho_p(2)-xrho_p(1));
    scale_y = (yrho_c(2)-yrho_c(1))/(yrho_p(2)-yrho_p(1));

    % Notes:
    % L corresponds to X_rho and longitude
    % M corresponds to Y_rho and latitude
    Lp=length(parent_grid_nc('xi_rho'));
    Mp=length(parent_grid_nc('eta_rho'));
    Lpc=length(child_grid_nc('xi_rho'));
    Mpc=length(child_grid_nc('eta_rho'));

    % Sanity check: ensure latitude&longitude of the child are within
    % parent
    latmin_c = child_lat(1);
    latmax_c = child_lat(end);
    lonmin_c = child_lon(1);
    lonmax_c = child_lon(end);
    latmin_p = parent_lat(1);
    latmax_p = parent_lat(end);
    lonmin_p = parent_lon(1);
    lonmax_p = parent_lon(end);
    if (latmin_c < latmin_p || latmax_c > latmax_p)
        error('Child latitude is outside of parent')
    else if (lonmin_c < lonmin_p || lonmax_c > lonmax_p)
            error('Child longitude is outside of parent')
    end
    if abs(scale_x-scale_y)>1e-4
        error('X and Y scaling are not the same')
    end

    imin = find(parent_lon>lonmin_c, 1)-1;
    imax = find(parent_lon>lonmax_c, 1)-1;
    jmin = find(parent_lat>latmin_c, 1)-1;
    jmax = find(parent_lat>latmax_c, 1)-1;
    refinecoeff = 1/scale_x;

end

