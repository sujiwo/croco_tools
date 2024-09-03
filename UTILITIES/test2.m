nc=nc4.netcdf.create('/tmp/randomval.ncf');
nc.redef();

nc('lonT')=10;
nc('latT')=10;
nc('lonU')=10;
nc('latU')=10;
nc('lonV')=10;
nc('latV')=10;
        
nc{'temp'}=nc4.ncdouble('lonT','latT');
nc{'temp'}(:)=0;
nc{'temp'}.createAttribute('long_name', 'TEMPERATURE');

nc.endef();
nc.close();