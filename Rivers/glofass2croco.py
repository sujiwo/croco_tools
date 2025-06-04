import warnings
import netCDF4 as nc4
import numpy as np
import rasterio as rio
import rasterio.plot
import matplotlib.patches as mpatches
from matplotlib.pyplot import get_cmap, legend



class CrocoRunoff:
    def __init__(self, roPath, mode='r'):
        self.fd = nc4.Dataset(roPath, mode)
        
    def __del__(self):
        self.fd.close()
    
    @staticmethod
    def create(runoffPath):
        fd = nc4.Dataset(runoffPath, mode='w')
        fd.type = 'CROCO runoff file'
        
        fd.createDimension('qbar_time')
        fd.createDimension('n_qbar')
        fd.createDimension('one', 1)
        fd.createDimension('two', 2)
        
        __=fd.createVariable('qbar_time', 'f8', dimensions=('qbar_time'))
        __.long_name = 'runoff time'
        __.units = 'days'
        __.cycle_length = 365.25
        
        __=fd.createVariable('Qbar', 'f8', dimensions=('n_qbar', 'qbar_time'))
        __.long_name = 'runoff discharge'
        __.units = 'm3.s-1'
        
        fd.close()
        return CrocoRunoff(runoffPath, mode='a')


class CGrid:
    def __init__(self, grPath):
        self.fd = nc4.Dataset(grPath)
        self.latitudes_v = np.array( self.fd['lat_v'][:,0] )
        self.longitudes_v = np.array( self.fd['lon_v'][0,:] )
        self.latitudes_u = np.array( self.fd['lat_u'][:,0] )
        self.longitudes_u = np.array( self.fd['lon_u'][0,:] )
        self.latitudes_r = np.array( self.fd['lat_rho'][:,0] )
        self.longitudes_r = np.array( self.fd['lon_rho'][0,:] )
        self.latitudes_p = np.array( self.fd['lat_psi'][:,0] )
        self.longitudes_p = np.array( self.fd['lon_psi'][0,:] )
        
        # build rio data for mask_rho
        self.transform_r = rio.transform.from_origin(
            west=self.longitudes_r[0],
            north=self.latitudes_r[-1],
            xsize=self.longitudes_r[1]-self.longitudes_r[0],
            ysize=self.latitudes_r[1]-self.latitudes_r[0]
            )
        
        # Mask Notes: False==Land, True==Sea
        self.mask_r = np.array(self.fd['mask_rho'][:], dtype=bool)
        self.mask_u = np.array(self.fd['mask_u'][:], dtype=bool)
        self.mask_v = np.array(self.fd['mask_v'][:], dtype=bool)
        self.mask_p = np.array(self.fd['mask_psi'][:], dtype=bool)
        self.bathy = np.array(self.fd['h'])
        
    def get_mask(self, kind='rho'):
        match kind:
            case 'rho':
                table = self.mask_r
            case 'psi':
                table = self.mask_psi
            case 'u':
                table = self.mask_u
            case 'v':
                table = self.mask_v
            case _:
                raise ValueError("Unknown kind")
        return np.flipud(table)

    def get_lon_limit(self):
        lonmin = max(self.longitudes_p[0], self.longitudes_r[0],
                     self.longitudes_u[0], self.longitudes_v[0])
        lonmax = min(self.longitudes_p[-1], self.longitudes_r[-1],
                     self.longitudes_u[-1], self.longitudes_v[-1])
        return lonmin, lonmax
    
    def get_lat_limit(self):
        latmin = max(self.latitudes_p[0], self.latitudes_r[0],
                     self.latitudes_u[0], self.latitudes_v[0])
        latmax = min(self.latitudes_p[-1], self.latitudes_r[-1],
                     self.latitudes_u[-1], self.latitudes_v[-1])
        return latmin, latmax
        
    def plot_mask_rho(self):
        rasterio.plot.show(np.flipud(self.mask_r), transform=self.transform_r)
        cmap = get_cmap()
        patches = [mpatches.Patch(color=cmap(0.0), label="Land"),
                   mpatches.Patch(color=cmap(1.0), label="Sea")]
        legend(handles=patches)
        
    def _select_table(self, kind):
        match kind:
            case 'rho':
                lon_tg = self.longitudes_r
                lat_tg = self.latitudes_r
                table = self.mask_r
            case 'psi':
                lon_tg = self.longitudes_p
                lat_tg = self.latitudes_p
                table = self.mask_psi
            case 'u':
                lon_tg = self.longitudes_u
                lat_tg = self.latitudes_u
                table = self.mask_u
            case 'v':
                lon_tg = self.longitudes_v
                lat_tg = self.latitudes_v
                table = self.mask_v
            case _:
                raise ValueError("Unknown kind")
        return (lat_tg, lon_tg, table)

    def where(self, lat, lon, kind='rho'):
        lat_tg, lon_tg, _ = self._select_table(kind)
        lon_num = np.searchsorted(lon_tg, lon)
        lat_num = np.searchsorted(lat_tg, lat)
        return (lon_num,lat_num)
    
    def query(self, lat, lon, kind='rho'):
        lat_tg, lon_tg, table = self._select_table(kind)
        lonnum = np.searchsorted(lon_tg, lon)
        latnum = np.searchsorted(lat_tg, lat)
        return table[latnum, lonnum]


class Glofass:
    def __init__(self, glPath):
        self.fd = nc4.Dataset(glPath)
        # Filter all valid sources
        discharges = self.fd['dis24'][0,:,:]
        coordshapes = discharges.shape
        self.truerivers = []
        self.valid_times = np.array(self.fd['valid_time'][:], dtype=self.fd['valid_time'].dtype)
        self.latitude = np.array(self.fd['latitude'])
        self.longitude = np.array(self.fd['longitude'])
        for lat_n in range(coordshapes[0]):
            for lon_n in range(coordshapes[1]):
                if discharges.mask[lat_n, lon_n]==True:
                    continue
                latitude = float(self.fd['latitude'][lat_n])
                longitude = float(self.fd['longitude'][lon_n])
                self.truerivers.append({'lat':latitude, 
                    'lon':longitude, 
                    'gindex': (lat_n,lon_n)}
                )

    def filter_rivers(self, lat_min, lat_max, lon_min, lon_max):
        self.truerivers = [R for R in self.truerivers 
                           if R['lat']>=lat_min and R['lat']<=lat_max
                           and R['lon']>=lon_min and R['lon']<=lon_max
                          ]
        self.latitude = self.latitude[
            (self.latitude>=lat_min) & (self.latitude<=lat_max)
            ]
        self.longitude = self.longitude[
            (self.longitude>=lon_min) & (self.longitude<=lon_max)
            ]
        
    def convert_croco(self, cgrid: CGrid, crocoRunoffPath: str):
        lonmin,lonmax = cgrid.get_lon_limit()
        latmin,latmax = cgrid.get_lat_limit()
        self.filter_rivers(latmin, latmax, lonmin, lonmax)

        P = []
        for R in self.truerivers:
            lon_num, lat_num = cgrid.where(R['lat'], R['lon'], kind='rho')
            if lat_num>=cgrid.mask_u.shape[0]-1 or \
                lat_num>=cgrid.mask_v.shape[0]-1:
                    continue
            if lon_num>=cgrid.mask_u.shape[1]-1 or \
                lon_num>=cgrid.mask_v.shape[1]-1:
                    continue
            if cgrid.mask_r[lat_num,lon_num]==False and \
                cgrid.mask_u[lat_num,lon_num]==False and \
                cgrid.mask_u[lat_num,lon_num+1]==False and \
                cgrid.mask_v[lat_num,lon_num]==False and \
                cgrid.mask_v[lat_num+1,lon_num]==False:
                continue
            if cgrid.mask_r[lat_num,lon_num]==True and \
                cgrid.mask_u[lat_num,lon_num]==True and \
                cgrid.mask_u[lat_num,lon_num+1]==True and \
                cgrid.mask_v[lat_num,lon_num]==True and \
                cgrid.mask_v[lat_num+1,lon_num]==True:
                continue

            # Find position and direction in grid
            # if mask change is in U axis => zonal
            # if mask change is in V axis => meridional
            if cgrid.mask_u[lat_num,lon_num]!= \
                cgrid.mask_u[lat_num,lon_num+1]:
                    R['type']='U'
                    if cgrid.mask_u[lat_num,lon_num]==False:
                        # Left-to-right
                        R['dir']=1
                        R['lat_num']=lat_num+1
                        R['lon_num']=lon_num+1
                    else:
                        # Right-to-left
                        R['dir']=-1
                        R['lat_num']=lat_num+1
                        R['lon_num']=lon_num+1+1
            elif cgrid.mask_v[lat_num,lon_num]!= \
                cgrid.mask_v[lat_num+1,lon_num]:
                    R['type']='V'
                    if cgrid.mask_v[lat_num,lon_num]==False:
                        # Bottom-to-top
                        R['dir']=1
                        R['lat_num']=lat_num+1
                        R['lon_num']=lon_num+1
                    else:
                        # Top-to-bottom
                        R['dir']=-1
                        R['lat_num']=lat_num+1+1
                        R['lon_num']=lon_num+1
            else:
                continue
            
            # Replace latitude using rho points
            R['lat'] = cgrid.latitudes_r[lat_num]
            R['lon'] = cgrid.longitudes_r[lon_num]
            glat = R['gindex'][0]
            glon = R['gindex'][1]
            # Collect flows, filter those with average discharge below cutoff
            flow = glo.fd['dis24'][:,glat,glon]
            flowAvg = np.average(flow)
            R['disAvg'] = flowAvg
            if flowAvg < dischargeCutOff:
                continue
            
            P.append(R)
        self.truerivers = P
        
        # Create croco runoff file
        croco_runoff = CrocoRunoff.create(crocoRunoffPath)
        Qbar = np.zeros((len(self.truerivers), len(self.valid_times)))
        for i in range(len(self.truerivers)):
            R = self.truerivers[i]
            cur_flows = self.fd['dis24'][:,R['gindex'][0],R['gindex'][1]]
            Qbar[i] = cur_flows
        croco_runoff.fd['Qbar'][:] = Qbar

        
        # Filter
        
        # Write times
        
        # ...
        pass


if __name__=='__main__':
    glosrcname = '/home/sujiwo/Data/Natuna/DATA/glofas_merged_NATUNA.nc'
    crocdstname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_runoff_test.nc'
    gridname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_grd.nc'
    dischargeCutOff = 1.0

    cgrid = CGrid(gridname)
    glo = Glofass((glosrcname))
    glo.convert_croco(cgrid, crocdstname)

    pass
