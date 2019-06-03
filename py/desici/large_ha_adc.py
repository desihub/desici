import argparse
import json
import matplotlib.pyplot as plt
import numpy as np

# +/- 4 hours
ha_vals = 15*np.array([-4.0, 4.0])
dec_vals = [-10, 0, 10, 20, 30, 40]

# probably better to start/end closer to the KPNO's latitude
dec_vals.reverse()

def fvc_request():
    req = {"sequence": "FVC",
           "flavor": "science",
           "exptime": 1.0,
           "fiducials": "on",
           "leave_fiducials": "off",
           "program": "FVC astrometry fields with ADC"}

    return req

def one_image_request(reqha, reqdec, exptime=60.0):

    request = {"sequence": "CI",
             "flavor": "science",
             "exptime": exptime,
             "reqha" : reqha,
             "reqdec" : reqdec,
             "track": True,
             "useadc": True,
             "correct_for_adc": True,
             "program": "astrometry fields with ADC"}

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
        plt.savefig('astrometry_fields_with_adc.png', bbox_inches='tight')

    return ha_sequence, dec_sequence

def all_observations(exptime=60.0, do_plot=False):
    ha, dec = ha_dec_grid(do_plot=do_plot)

    requests = []
    for t in zip(ha, dec):
        requests.append(one_image_request(t[0], t[1], exptime=exptime))
        requests.append(fvc_request())

    return requests

def total_time_estimate(requests):
    t = 0.0 # seconds

    # based on recent experience the full readout/OCS
    # overhead is fractionally a significant amount larger than the 
    # nominal 6.5 second readout time
    # so this isn't really the detector readout time, but is more 
    # realistic for calculating overheads

    # assume this value is also valid for FVC

    t_readout = 10.0 # seconds

    # just a guess

    t_settle = 5.0 # seconds

    for request in requests:
        t += (request['exptime'] + t_readout)

    n_slews = len(ha_vals)*len(dec_vals) - 1
    
    t += t_settle*n_slews

    slew_tot_ha_deg = np.max(ha_vals)-np.min(ha_vals)
    slew_tot_dec_deg = (len(np.unique(ha_vals)) - 1)*(np.max(dec_vals)-np.min(dec_vals))

    seconds_per_deg = 2.7

    t += seconds_per_deg*(slew_tot_ha_deg + slew_tot_dec_deg)

    t_minutes = t/60.0

    return t_minutes

if __name__ == "__main__":
    descr = 'create script to take exposures at large HA and moderate Dec correcting for ADC'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--exptime', default=60.0, type=float,
                        help='exposure time at each pointing, default 60')

    parser.add_argument('--outname', default='large_ha_adc.json',
                        help='output file name for observing script')

    parser.add_argument('--do_plot', default=False, action='store_true',
                        help='save plot showing sequence of pointings')

    args = parser.parse_args()

    requests = all_observations(exptime=args.exptime, do_plot=args.do_plot)

    outname = args.outname

    with open(outname, 'w') as outfile:
        json.dump(requests, outfile, indent=2)

    print('CI exposure sequence will take ' + "{:.1f}".format(total_time_estimate(requests)) + ' minutes')
