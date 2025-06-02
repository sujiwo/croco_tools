import netCDF4 as nc4
import numpy as np
import rasterio as rio
import rasterio.plot
from copy import copy
import warnings


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
                    'index': (lat_n,lon_n)}
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
        self.mask_r = np.array(np.flipud(self.fd['mask_rho']), dtype=bool)
        self.mask_u = np.array(np.flipud(self.fd['mask_u']), dtype=bool)
        self.mask_v = np.array(np.flipud(self.fd['mask_v']), dtype=bool)
        self.mask_p = np.array(np.flipud(self.fd['mask_psi']), dtype=bool)
        self.bathy = np.flipud(self.fd['h'])
        
    def plot_mask_rho(self):
        rasterio.plot.show(self.mask_r, transform=self.transform_r)
        
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
        xd = np.searchsorted(lon_tg, lon)
        yd = np.searchsorted(lat_tg, lat, side='right')
        return (xd,yd)
    
    def query(self, lat, lon, kind='rho'):
        lat_tg, lon_tg, table = self._select_table(kind)
        xd = np.searchsorted(lon_tg, lon)
        yd = np.searchsorted(lat_tg, lat, side='right')
        return table[yd, xd]

def convert_glofass2croco(glosrc_, crocogrd_):
    glofd = nc4.Dataset(glosrc_)
    cgrid = CGrid(crocogrd_)
    
    if glofd['latitude'][0] < cgrid.latitudes_r[0] or \
        glofd['latitude'][-1] > cgrid.latitudes_r[-1]:
            warnings.warn("Glofass latitude is outside grid range")
    if glofd['longitude'][0] < cgrid.longitudes_r[0] or \
        glofd['longitude'][-1] > cgrid.longitudes_r[-1]:
            warnings.warn("Glofass longitude is outside grid range")

    discharges = glofd['dis24'][0,:,:]
    coordshapes = discharges.shape
    truerivers = []
    for lat_n in range(coordshapes[0]):
        for lon_n in range(coordshapes[1]):
            if discharges.mask[lat_n, lon_n]==True:
                continue
            latitude = float(glofd['latitude'][lat_n])
            longitude = float(glofd['longitude'][lon_n])
            # Find position in domain grid
            lon_rg, lat_rg = cgrid.where(latitude, longitude)
            # if cgrid.mask_r[lat_rg, lon_rg]==False:
            pass


if __name__=='__main__':
    glosrcname = '/home/sujiwo/Data/Natuna/DATA/glofas_merged_NATUNA.nc'
    crocdstname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_runoff_test.nc'
    gridname = '/home/sujiwo/Data/Natuna/CROCO_FILES/croco_grd.nc'

    cgrid = CGrid(gridname)
    glo = Glofass((glosrcname))
    glo.filter_rivers(cgrid.latitudes_r[0], cgrid.latitudes_r[-1], 
                      cgrid.longitudes_r[0], cgrid.longitudes_r[-1])

    # convert_glofass2croco(glosrcname, gridname)
    pass
# glosrc = nc4.Dataset(glosrcname)
# grid = nc4.Dataset(gridname)


# discharges = glosrc['dis24'][0,:,:]
# print("There are {} suspected river(s)".format(discharges.count()))
# coordshapes = discharges.shape
# truerivers = []

# xc = 0
# latitudes_v = grid['lat_v'][:,0]
# longitudes_v = grid['lon_v'][0,:]
# latitudes_u = grid['lat_u'][:,0]
# longitudes_u = grid['lon_u'][0,:]
# latitudes_r = grid['lat_rho'][:,0]
# longitudes_r = grid['lon_rho'][0,:]

# for j in range(coordshapes[0]):
#     for i in range(coordshapes[1]):
#         if discharges.mask[j,i]==True:
#             continue
#         latitude = float(glosrc['latitude'][j])
#         longitude = float(glosrc['longitude'][i])
#         # Find position in domain grid
#         ri = np.searchsorted(longitudes_r, longitude)
#         rj = np.searchsorted(latitudes_r, latitude)
#         if grid['mask_rho'][rj,ri]==0 and \
#            grid['mask_v'][rj,ri]==0 and \
#            grid['mask_v'][rj+1,ri]==0 and \
#            grid['mask_u'][rj,ri]==0 and \
#            grid['mask_u'][rj,ri+1]==0:
#                print('In sea: ({},{})'.format(i,j))
        

#glosrc.close()
#grid.close()