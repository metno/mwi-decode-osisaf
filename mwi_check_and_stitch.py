'''Control routine to figure out which files to call and call mwu_pdu_stitch on these'''

import argparse
import re
import os

from mwi_pdu_stitch import mwi_pdu_stitch

orbpatt = '^mwi_sgb1_n{}_.*\.nc$'
patt = re.compile(orbpatt.replace('{}', '(\d{1,})'))

try:
    c_rows, c_cols = os.popen('stty size', 'r').read().split()
except ValueError as e:
    c_cols = 100


def parse_args():

    parser = argparse.ArgumentParser(description='mwi_check_and_stitch',
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=40, width=int(c_cols)))
    parser.add_argument('-i', '--indir', required=True,
                        default=None,
                        help='Input dir')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Output dir')
    parser.add_argument('-O', '--overwrite', dest='owrite', action='store_true',
                        default=False,
                        help='Overwrite existing outputfile')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        default=False,
                        help='Print more stuff to stdout')
    args = parser.parse_args()

    return args


def sort_files_by_timestamp(tlist):

    tstamps = [os.stat(t).st_mtime for t in tlist]
    zlist = sorted(zip(tstamps, tlist))

    return zlist


def mwi_check_and_stitch():

    args = parse_args()
    indir = args.indir
    outdir = args.outdir

    # Checking the input and output directories
    if not os.path.isdir(indir):
        print('\nInvalid input dir: {}'.format(indir))
        sys.exit(1)

    if not os.path.isdir(outdir):
        print('\nInvalid output dir: {}'.format(outdir))
        sys.exit(1)

    if os.path.dirname(indir) == os.path.abspath(outdir):
        print('\nCannot have same input and output dir\n')
        sys.exit(1)

    # Finding a file list
    flist = os.listdir(indir)
    for f in flist:
        if not patt.match(f):
            print('Skipping {}',format(f))
            continue
    flist = [os.path.join(indir, x) for x in flist]
    zflist = sort_files_by_timestamp(flist)

    # Check the last modified time of the output files
    outmodtime = None
    olist = os.listdir(outdir)
    for o in olist:
        if not o.startswith('mwi') or o[-3:] != '.nc':
            continue
    # There may not be output files
    if olist:
        olist = [os.path.join(indir, x) for x in olist]
        zolist = sort_files_by_timestamp(olist)
        outmodtime = zolist[-1][0]

    # Select the decoded files that are new from the zipped list, and
    # return to a single list of files
    if outmodtime:
        slist = []
        for z in zflist:
            # TODO: Should probably add a small time difference so files which
            # come in simultaneously with the last stitching get included
            if z[0] > outmodtime:
                slist.append(z[1])
    else:
        slist = [s[1] for s in zflist]

    # If the list of new files is empty, should exit at this point
    if not slist:
        print("No new files found")
        sys.exit(1)

    # Logic to select based on orbit number. Find the number of unique
    # orbit numbers that have arrived since the last created output file
    nlist = [patt.match(os.path.basename(s))[1] for s in slist]
    orbnum_to_rerun = set(nlist)

    # Find all files with each orbit number
    for orbnum in orbnum_to_rerun:
        opatt = re.compile(orbpatt.replace('{}', orbnum))
        rlist = []
        for f in flist:
            if opatt.match(os.path.basename(f)):
                rlist.append(f)

        # Calling the stitching routine
        mwi_pdu_stitch(sorted(rlist), outdir, overwrite=args.owrite,
                       verbose=args.verbose)


if __name__ == '__main__':

    mwi_check_and_stitch()
