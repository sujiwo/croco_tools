function ...
[ichildgrd_p,jchildgrd_p,...
ichildgrd_r,jchildgrd_r,...
ichildgrd_u,jchildgrd_u,...
ichildgrd_v,jchildgrd_v] = create_query_grid(imin,imax,jmin,jmax,refinecoeff,Lpc,Mpc)

    ipchild=linspace(imin,imax,Lpc-1);
    jpchild=linspace(jmin,jmax,Mpc-1);
    irchild = linspace(imin+0.5-0.5/refinecoeff,imax+0.5+0.5/refinecoeff,Lpc);
    jrchild = linspace(jmin+0.5-0.5/refinecoeff,jmax+0.5+0.5/refinecoeff,Mpc);
    [ichildgrd_p,jchildgrd_p]=meshgrid(ipchild,jpchild);
    [ichildgrd_r,jchildgrd_r]=meshgrid(irchild,jrchild);
    [ichildgrd_u,jchildgrd_u]=meshgrid(ipchild,jrchild);
    [ichildgrd_v,jchildgrd_v]=meshgrid(irchild,jpchild);

