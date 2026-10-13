"""Save bounded true projected positions for the three illustration systems."""
import socket
import numpy as np
from measure_halo import ROOT,load,sb,rv
if socket.gethostname().lower()!='lageunha':raise SystemExit('Lageunha only')
for step,nstar in [('00475',5557744),('00967',17494810),('00178',7062464)]:
    dest=ROOT/'results'/f'{step}_{nstar}'
    h=load(step,nstar);scale,n=rv.angular_scale(h['a'],h['z'])
    with np.load(dest/'images.npz') as d:sample=d['sample_indices']
    coords=[h[c][sample]*scale for c in ['dx','dy','dz']]
    arrays={}
    for los,(j,k) in zip(rv.LOS,sb.LOS_AXES):arrays['x_'+los]=coords[j];arrays['y_'+los]=coords[k]
    np.savez_compressed(dest/'drawing_positions.npz',**arrays)
    del h,coords,arrays
    print('SAVED_PROJECTED_POSITIONS',step,nstar,flush=True)
