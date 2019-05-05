import argparse
import json
from subprocess import Popen, PIPE


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Returns json with chromsizes of given genome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-genome",
                        type=str,
                        nargs='?',
                        dest='genome',
                        help='genome file')
    parser.add_argument("-outfile",
                        type=str,
                        nargs='?',
                        dest='outfilename',
                        help='output json filename')
    args = parser.parse_args()
    awk_command = '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
    proc = Popen("cat {} | awk '{}'".format(args.genome, awk_command),
                 stdout=PIPE,
                 stderr=PIPE,
                 shell=True)
    stdout, stderr = proc.communicate()
    if stderr:
        print(stderr.decode())
    chromsizes = dict([(line.split()[0], line.split()[-1])
                       for line in stdout.decode().strip().split("\n")])
    with open(args.outfilename, 'w') as outfile:
        json.dump(chromsizes, outfile)
