import netCDF4 as nc4
import numpy as np
import rasterio as rio
import rasterio.plot
from copy import copy
import warnings
import matplotlib.patches as mpatches
from matplotlib.pyplot import get_cmap, legend

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


class CrocoRunoff:
    def __init__(self, roPath):
        pass
    
    @staticmethod
    def create(runoffPath):
        fd = nc4.Dataset(runoffPath, mode='w')
        fd.createDimension('qbar_time')
        fd.createDimension('n_qbar')
        fd.createDimension('one', 1)
        fd.createDimension('two', 2)
        
        qb=fd.createVariable('qbar_time', 'f8', dimensions=('qbar_time'))
        qb.long_name = 'runoff time'
        qb.units = 'days'
        qb.cycle_length = 365.25
        
        fd.close()
        


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

def convert_glofass2croco(glosrc_, crocogrd_):
    glofd = nc4.Dataset(glosrc_)
    cgrid = CGrid(crocogrd_)
    
    if glofd['latitude'][0] < cgrid.latitudes_r[0] or \
        glofd['latitude'][-1] > cgrid.latitudes_r[-1]:
            warnings.warn("Glofass latitude is outside grid range")
    if glofd['longitude'][0] < cgrid.longitudes_r[0] or \
        glofd['longitude'][-1] > cgrid.longitudes_r[-1]:
            warnings.warn("Glofass longitude is outside grid range")

    glofd.filter_rivers(cgrid.latitudes_r[0], cgrid.latitudes_r[-1],
                        cgrid.longitudes_r[0], cgrid.longitudes_r[-1])
    # P = []
    # for R in glofd.truerivers:
    #     xr, yr = cgrid.where(R['lat'], R['lon'], kind='rho')
    #     if cgrid.mask_r[yr,xr]==False and \
    #         cgrid.mask_u[yr,xr]==False:
    #         pass


if __name__=='__main__':
    glosrcname = '/home/sujiwo/Data/Natuna/DATA/glofas_merged_NATUNA.nc'
    crocdstname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_runoff_test.nc'
    gridname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_grd.nc'
    dischargeCutOff = 7.5

    cgrid = CGrid(gridname)
    glo = Glofass((glosrcname))
    lonmin,lonmax = cgrid.get_lon_limit()
    latmin,latmax = cgrid.get_lat_limit()
    glo.filter_rivers(latmin, latmax, lonmin, lonmax)

    P = []
    for R in glo.truerivers:
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
        
        # Find position in grid and 
        
        P.append(R)

    pass
