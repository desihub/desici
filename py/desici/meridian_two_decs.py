import argparse
import json

def fvc_request():
    req = {"sequence": "FVC",
           "flavor": "science",
           "exptime": 1.0,
           "fiducials": "on",
           "leave_fiducials": "off",
           "program": "FVC meridian field rotation"}
    return req

def one_observation(ra_dest, dec_dest, exptime, pointing_exptime, visit_num=None, do_fvc=False):

    slew =   {"sequence": "Action",
              "action": "slew",
              "reqra": ra_dest,
              "reqdec": dec_dest,
              "track": True}

    pointing = {"sequence": "CI",
                "flavor": "science",
                "exptime": pointing_exptime,
                "program": "meridian field rotation two declinations pointing"}

    _break = {"sequence": "Break"}

    image = {"sequence": "CI",
             "flavor": "science",
             "exptime": exptime,
             "program": "meridian field rotation two declinations"}

    if visit_num is not None:
        image["object"] = ("north" if (dec_dest > 30) else "south") + " visit " + str(visit_num)

    requests = [slew, pointing, _break, image]

    if do_fvc:
        requests.append(fvc_request())

    return requests

def all_observations(visits_per_field, lst_deg, dec_north_deg, dec_south_deg, exptime, pointing_exptime, do_fvc=False):

    # want to go 1.5 East in RA relative to LST at start of test
    # also pad by 0.02 hours to account for time taken by very first slew
    pad_hours = 0.02
    ra_deg = lst_deg + (1.5 + pad_hours)*15

    requests = []
    for i in range(visits_per_field):
        for dec in [dec_south_deg, dec_north_deg]:
            requests += one_observation(ra_deg, dec, exptime, pointing_exptime, visit_num=(i+1), do_fvc=do_fvc)

    return requests

if __name__ == "__main__":
    descr = 'create script to watch two fields on the meridian for three hours'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('lst_degrees', type=float, nargs=1, 
                        help='LST in decimal degrees at start of test')

    parser.add_argument('--dec_north_degrees', default=60.0,
                        help='northern Dec in decimal degrees for CI exposures, default 60')

    parser.add_argument('--dec_south_degrees', default=10.0,
                        help='southern Dec in decimal degrees for CI exposures, default 10')

    parser.add_argument('--exptime', default=60.0,
                        help='exposure time for CI sky image in seconds, default 60')

    parser.add_argument('--pointing_exptime', default=15.0,
                        help='exposure time for pointing check after slewing in seconds, default 15')

    parser.add_argument('--visits_per_field', default=25,
                        help='number of visits per field, default of 25 errs larger than necessary')

    parser.add_argument('--outname', default='meridian_two_decs.json',
                        help='output file name for observing script')

    parser.add_argument('--do_fvc', default=False, action='store_true',
                        help='take FVC image during each visit to each field')

    args = parser.parse_args()

    requests = all_observations(args.visits_per_field, args.lst_degrees[0], 
                                args.dec_north_degrees, args.dec_south_degrees,
                                args.exptime, args.pointing_exptime, do_fvc=args.do_fvc)

    outname = args.outname
    with open(outname, 'w') as outfile:
        json.dump(requests, outfile, indent=2)
