import argparse
import json
import matplotlib.pyplot as plt
import numpy as np

# +/- 3 hours, +/- 2.5 hours, +/- 2 hours, +/- 1.5 hours, +/- 1 hour and zenith
ha_vals = 15*np.array([-3.0, -2.5, -2.0, -1.5, -1.0, 0.0, 1.0, 1.5, 2.0, 2.5, 3.0])
dec_vals = [60.0, 70.0, 75.0]

def fvc_request():
    req = {"sequence": "FVC",
           "flavor": "science",
           "exptime": 1.0,
           "fiducials": "on",
           "leave_fiducials": "off",
           "program": "FVC astrometry fields"}

    return req

def one_image_request(reqha, reqdec, exptime=60.0):

    request = {"sequence": "CI",
             "flavor": "science",
             "exptime": exptime,
             "reqha" : reqha,
             "reqdec" : reqdec,
             "track": True,
             "useadc": False,
             "program": "astrometry fields"}

    return request

def ha_dec_grid(do_plot=False):

    _dec_vals = dec_vals

    ha_sequence = []
    dec_sequence = []
    for ha in ha_vals:
        for dec in _dec_vals:
            ha_sequence.append(ha)
            dec_sequence.append(dec)
        _dec_vals.reverse()

    if do_plot:
        plt.plot(ha_sequence, dec_sequence)
        plt.scatter(ha_sequence, dec_sequence)
        plt.title('lines connect time-adjacent exposures')
        plt.xlabel('HA (deg)')
        plt.ylabel('Dec (deg)')
        plt.savefig('astrometry_fields.png', bbox_inches='tight')

    return ha_sequence, dec_sequence

def all_observations(exptime=60.0, do_plot=False):
    ha, dec = ha_dec_grid(do_plot=do_plot)

    requests = []
    for t in zip(ha, dec):
        requests.append(one_image_request(t[0], t[1], exptime=exptime))
        requests.append(fvc_request())

    return requests

if __name__ == "__main__":
    descr = 'create script to take exposures at a high Dec grid of (HA, Dec)'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--exptime', default=60.0, type=float,
                        help='exposure time at each pointing, default 60')

    parser.add_argument('--outname', default='ha_dec_grid.json',
                        help='output file name for observing script')

    parser.add_argument('--do_plot', default=False, action='store_true',
                        help='save plot showing sequence of pointings')

    args = parser.parse_args()

    requests = all_observations(exptime=args.exptime, do_plot=args.do_plot)

    outname = args.outname

    with open(outname, 'w') as outfile:
        json.dump(requests, outfile, indent=2)
