tic
nc = netcdf('/home/sujiwo/Data/LombokProj/Lombok3/CROCO_FILES/croco_grd.nc');
v = nc{'angle'}(100:110,50:59);
toc

tic
nc2 = nc4.netcdf('/home/sujiwo/Data/LombokProj/Lombok3/CROCO_FILES/croco_grd.nc');
v1 = nc2.var('angle').get(100:110,50:59);
toc

nc.close();
nc2.close();