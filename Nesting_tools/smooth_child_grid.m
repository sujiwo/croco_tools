function smooth_child_grid(parent_grid, child_grid)

    %crocotools_param;

    parent_nc = netcdf(parent_grid);
    child_nc  = netcdf(child_grid, 'write');
    h_parent  = parent_nc{'h'}(:);
    [Mp,Lp]   = size(h_parent);
    Lpc=length(child_nc('xi_rho'));
    Mpc=length(child_nc('eta_rho'));
    [imin,imax,jmin,jmax,refinecoeff] = detect_grid(parent_nc,child_nc);

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
    [ichildgrd_p,jchildgrd_p,...
    ichildgrd_r,jchildgrd_r,...
    ichildgrd_u,jchildgrd_u,...
    ichildgrd_v,jchildgrd_v] = create_query_grid(imin,imax,jmin,jmax,refinecoeff,Lpc,Mpc);

    % Interpolation
    hrawchild =interp2(igrd_r,jgrd_r,h_parent,ichildgrd_r,jchildgrd_r,'cubic');
    child_nc{'hraw'}(:) = hrawchild;
    close(child_nc);
    close(parent_nc);

end
