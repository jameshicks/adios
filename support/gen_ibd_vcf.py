from __future__ import print_function
import argparse
import sys
import scipy.stats
import numpy as np


def make_size(v):
    suffix_size = {'kb': 1000, 'mb': 1000000, 'bp': 1}
    if v[-2:] in suffix:
        return int(float(v[:-2]) * suffix_size[v[-2:]])
    else:
        return int(v)


def make_vcf_allele(a, b): return '{}|{}'.format(a, b)
parser = argparse.ArgumentParser()
parser.add_argument('--size', type=make_size, help='Sequence length',
                    default=150000000)
parser.add_argument('--lambd', type=float, default=0.05,
                    help='Exponential dist parameter')
parser.add_argument('--dist', type=int, default=250,
                    help='Distance between markers')
parser.add_argument('--out', default='IBD1_pair.vcf')
args = parser.parse_args()

ibd_start = int(args.size*0.25)
ibd_size = int(2e7)
ibd_stop = ibd_start+ibd_size


print("Generating {} markers.".format(int(args.size / args.dist)))
print("IBD region is {}-{} ({}bp)".format(ibd_start, ibd_stop, ibd_size))
print("Writing output to {}...".format(args.out))
with open(args.out, 'w') as of:
    print("##fileformat=VCF4.1", file=of)
    print("##kind=simulation", file=of)
    print("##source={}".format(' '.join(sys.argv)), file=of)
    print("##ibd=<Ind1=A,Ind2=B,Chr=1,start={},stop={},size={}>".format(
        ibd_start, ibd_stop, ibd_size), file=of)
    print("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele MAF\">",
          file=of)

    fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
              "INFO", "FORMAT", "DUMMY_A", "DUMMY_B"]
    print("#{}".format('\t'.join(fields)), file=of)
    for i, pos in enumerate(range(1, args.size, args.dist)):
        freq = np.random.exponential(args.lambd)
        ibd_state = ibd_start < pos < ibd_stop

        info = 'AF={0:.4f};IBD={1}'.format(freq, int(ibd_state))
        if ibd_state == 0:
            a1 = int(np.random.random() < freq)
            a2 = int(np.random.random() < freq)
            a = make_vcf_allele(a1, a2)

            b1 = int(np.random.random() < freq)
            b2 = int(np.random.random() < freq)
            b = make_vcf_allele(b1, b2)

        elif ibd_state == 1:
            # print('in state')
            shared_allele = int(np.random.random() < freq)
            a2 = int(np.random.random() < freq)
            b2 = int(np.random.random() < freq)

            a = make_vcf_allele(shared_allele, a2)
            b = make_vcf_allele(shared_allele, b2)

        record = ['1', str(pos), 'sim{}'.format(i),
                  'A', 'G', '255', 'PASS', info, 'GT',
                  a, b]
        print('\t'.join(str(x) for x in record), file=of)
    print('', file=of)
print("Done")
