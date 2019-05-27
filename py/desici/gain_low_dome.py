import argparse
import json
import numpy as np

exptimes = np.arange(2, 34, 2).astype(float)

def all_observations(n_per_exptime):
    # can I specify track=False in this request, or is that only valid for slews?

    requests = []
    for exptime in exptimes:
        request = {"sequence": "CI",
                   "flavor": "science",
                   "exptime": exptime,
                   "program": "gain low dome"}
        requests = requests + [request]*n_per_exptime

    return requests

if __name__ == "__main__":
    descr = 'create script to gather closed-dome data that can be used to measure gain'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--n_per_exptime', default=4, type=int,
                        help='number of exposures per EXPTIME value, default 4')

    parser.add_argument('--outname', default='gain_low_dome.json',
                        help='output file name for observing script')

    args = parser.parse_args()

    requests = all_observations(args.n_per_exptime)

    outname = args.outname

    with open(outname, 'w') as outfile:
        json.dump(requests, outfile, indent=2)

