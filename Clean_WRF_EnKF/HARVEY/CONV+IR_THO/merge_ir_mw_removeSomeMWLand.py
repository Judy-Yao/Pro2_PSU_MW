#!/work2/06191/tg854905/App/anaconda3/bin/python
import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from itertools import compress
import sys
# Do not filter MW
# Do not use IR



SFC_CHANNELS={
        'gmi_gpm_lf': [1,2,3,4,5,6,7,8,9],
        'ssmi_f15': [1,2,3,4,5,6,7],
        'ssmis_f16': [1,2,12,13,14,15,16,17,18],
        'ssmis_f17': [1,2,12,13,14,15,16,17,18],
        'ssmis_f18': [1,2,12,13,14,15,16,17,18],
        'atms_npp': [1,2,3,4,5,16],
        'amsr2_gcom-w1': [1,2,3,4,5,6,7,8,9,10,11,12,13,14],
        'saphir_meghat': [6],
        }

def read_IR(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    nobs = len(lines)
    ir = {}
    ir['lat'] = np.zeros(nobs)
    ir['lon'] = np.zeros(nobs)
    ir['tb'] = np.zeros(nobs)
    
    for i in range(nobs):
        data = lines[i].split()
        ir['lat'][i] = float(data[3])
        ir['lon'][i] = float(data[4])
        ir['tb'][i] = float(data[5])
    return ir

def read_MW(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    nobs = len(lines)
    mw = {}
    mw['lat'] = np.zeros(nobs)
    mw['lon'] = np.zeros(nobs)
    mw['tb'] = np.zeros(nobs)
    mw['sfc'] = np.zeros(nobs,dtype=bool)
    mw['fov'] = np.zeros(nobs)
    
    for i in range(nobs):
        data = lines[i].split()
        mw['lat'][i] = float(data[3])
        mw['lon'][i] = float(data[4])
        mw['tb'][i] = float(data[5])
        if data[1] in SFC_CHANNELS:
            #if int(data[2]) == SFC_CHANNELS[ data[1]]:
            if int(data[2]) in SFC_CHANNELS[ data[1]]:
                mw['sfc'][i] = True
        mw['fov'][i] = max(float(data[8]), float(data[9]) ) 
    return mw


def read_WRF(filename):
    print('entered read_WRF')
    wrf = {}
    with Dataset(filename, 'r') as f:
        print(f['XLAT'][0,10,0])
        wrf['DX'] = f.DX / 1000
        wrf['lat'] = f['XLAT'][0,:,0]
        wrf['lon'] = f['XLONG'][0,0,:]
        wrf['ocean'] = np.logical_and(f['LANDMASK'][0,:,:]==0, f['LAKEMASK'][0,:,:]==0)
    print(wrf['DX'])
    return wrf




#def find_valid_idx(ir_obs, mw, wrf, radius = 2):
def find_valid_idx(mw, wrf, radius = 2):
    nlat = len( wrf['lat'] )
    nlon = len( wrf['lon'] )
    latmin = np.min(wrf['lat'])
    latmax = np.max(wrf['lat'])
    lonmin = np.min(wrf['lon'])
    lonmax = np.max(wrf['lon'])
    
    n_mw = len(mw['tb'])
    use_mw = np.ones( shape=mw['tb'].shape, dtype=bool )
    for i in range(n_mw):
        if mw['lat'][i] <= latmin or  mw['lat'][i] >= latmax or  mw['lon'][i] <= lonmin or  mw['lon'][i] >= lonmax:
            use_mw[i] = False
            continue
        ilat = np.searchsorted( wrf['lat'], mw['lat'][i] )
        ilon = np.searchsorted( wrf['lon'], mw['lon'][i] )
        ilat0 = max( 0, ilat-radius )
        ilat1 = min( nlat, ilat+radius+1)
        ilon0 = max( 0, ilon-radius )
        ilon1 = min( nlon, ilon+radius+1)
        #if np.any( ir_obs[ilat0:ilat1,ilon0:ilon1] ) or  not np.all( wrf['ocean'][ilat0:ilat1,ilon0:ilon1] ):
        if (np.mean( wrf['ocean'][ilat0:ilat1,ilon0:ilon1]) < 0.1):
            use_mw[i] = False
            continue
    return use_mw

def filter_obs(filename_in, useobs, filename_out):
    with open(filename_in, 'r') as f:
        lines = f.readlines()
    validlines = list(compress(lines,useobs))
    with open(filename_out, 'w') as f:
        f.writelines(validlines)

def calc_clearmask_2d(wrf, ir, threshold):
    nlat = len( wrf['lat'] )
    nlon = len( wrf['lon'] )
    latmin = np.min(wrf['lat'])
    latmax = np.max(wrf['lat'])
    lonmin = np.min(wrf['lon'])
    lonmax = np.max(wrf['lon'])
    isclear_2d = np.zeros( (nlat, nlon), dtype=bool )
    
    isclear_1d = np.logical_and.reduce( (ir['tb']>threshold, ir['lat']>latmin, ir['lat']<latmax, ir['lon']>lonmin, ir['lon']<lonmax ) )
    idx_lat = np.searchsorted( wrf['lat'], ir['lat'][isclear_1d] )
    idx_lon = np.searchsorted( wrf['lon'], ir['lon'][isclear_1d] )
    for i in range(len(idx_lat)):
        isclear_2d[ idx_lat[i], idx_lon[i] ] = True
    return isclear_2d

def test():
    filename_IR = 'radiance_201708231300_so_orig'
    ir = read_IR(filename_IR)

    filename_WRF = 'fort.80011'
    wrf = read_WRF(filename_WRF)

    filename_MW = 'microwave_201708231300_so_orig'
    mw = read_MW(filename_MW)

    
    nlat = len( wrf['lat'] )
    nlon = len( wrf['lon'] )
    latmin = np.min(wrf['lat'])
    latmax = np.max(wrf['lat'])
    lonmin = np.min(wrf['lon'])
    lonmax = np.max(wrf['lon'])
    
    isclear_2d = calc_clearmask_2d(wrf, ir, 230)
    
    use_mw = find_valid_idx(isclear_2d, mw, wrf, radius=int(np.ceil(mw['fov'] / wrf['DX']) ))
    use_ir = find_valid_idx(np.logical_not(isclear_2d), ir, wrf, radius=0)
  
    filter_obs(filename_IR, use_ir.tolist(), filename_IR + '_mwir')
    filter_obs(filename_MW, use_mw.tolist(), filename_MW + '_mwir')

    ir_filtered = read_IR(filename_IR+'_mwir')
    mw_filtered = read_IR(filename_MW+'_mwir')

    f, ax=plt.subplots(3, 2, subplot_kw=dict(adjustable='box', aspect='equal'), sharex='all', sharey='all', figsize=(4,6), dpi=300)
    ax[0,0].scatter( ir['lon'], ir['lat'], s=0.3, c=ir['tb'], cmap='jet')
#    plt.colorbar()
    
    ax[0,1].pcolormesh( wrf['lon'], wrf['lat'], isclear_2d)
#    plt.colorbar()
    
    ax[1,0].scatter( ir['lon'][use_ir], ir['lat'][use_ir], s=0.3, c=ir['tb'][use_ir], cmap='jet')
#    plt.colorbar()
    
    ax[1,1].scatter( mw['lon'][use_mw], mw['lat'][use_mw], s=0.3, c=mw['tb'][use_mw], cmap='jet')
#    plt.colorbar()
    
    ax[2,0].scatter( ir_filtered['lon'], ir_filtered['lat'], s=0.3, c=ir_filtered['tb'], cmap='jet')
    ax[2,1].scatter( mw_filtered['lon'], mw_filtered['lat'], s=0.3, c=mw_filtered['tb'], cmap='jet')
    
    
    ax[0,0].set_xlim( (lonmin, lonmax) )
    ax[0,0].set_ylim( (latmin, latmax) )
    plt.show()

def driver(filename_ir_in, filename_mw_in, filename_wrf_in, filename_ir_out, filename_mw_out):
    print(filename_wrf_in)
    ir = read_IR(filename_ir_in)
    wrf = read_WRF(filename_wrf_in)
    mw = read_MW(filename_mw_in)
    
    #isclear_2d = calc_clearmask_2d(wrf, ir, 230)
    
    use_ir = np.ones( shape=ir['tb'].shape, dtype=bool )
    #use_mw = find_valid_idx(isclear_2d, mw, wrf)
    use_mw = find_valid_idx(mw, wrf)
    use_mw = np.logical_or( use_mw, np.logical_not(mw['sfc']) ) # DO not filter non-surface-affected MW obs
 
 
    filter_obs(filename_ir_in, use_ir.tolist(), filename_ir_out)
    filter_obs(filename_mw_in, use_mw.tolist(), filename_mw_out)

if __name__=='__main__':
    #test()
    driver( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] )

