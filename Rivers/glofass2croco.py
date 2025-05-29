import netCDF4 as nc4
import numpy as np
import rasterio as rio
import rasterio.plot
from copy import copy


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
        yd = np.searchsorted(lat_tg, lat)
        return (xd,yd)
    
    def query(self, lat, lon, kind='rho'):
        lat_tg, lon_tg, table = self._select_table(kind)
        xd = np.searchsorted(lon_tg, lon)
        yd = np.searchsorted(lat_tg, lat)
        return table[yd, xd]

# glosrcname = 'DATA/glofas_merged_NATUNA.nc'
# crocdstname = 'CROCO_FILES/croco_runoff_test.nc'
# gridname = 'CROCO_FILES/croco_grd.nc'

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