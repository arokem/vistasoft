import argparse as arg

import numpy as np
import scipy.io as sio

from dipy.segment.quickbundles import QuickBundles

parser = arg.ArgumentParser(description='Segment a mrDiffusion FiberGroup
object into clusters using QuickBundles (Garyfallidise et al. 2012)')

parser.add_argument('in_file', action='store', metavar='File', 
                    help='FG file (.mat)')

parser.add_argument('--dist_th', action='store', metavar='Float',
                    help='Distance threshold (default: 30.)', default=30.)

parser.add_argument('--pts',  action='store', metavar='Int',
                    help='Points', default=18))

params = parser.parse_args()


if __name__ == "__main__":
    fg_mat = sio.loadmat(params.in_file, squeeze_me=True)
    fg = fg_mat['fg']
    name_orig = str(fg['name'])
    fibers = fg['fibers'].item()
    streamlines = [np.array(ff).T for ff in fibers]
    qb = QuickBundles(streamlines, dist_thr=params.dist_th, pts=params.pts)
    clusters = qb.clusters()

    for c in clusters:
        new_fg = 


