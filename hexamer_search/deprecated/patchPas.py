import argparse, os, sys
from Kleat import *
import pysam

def revComp(seq):
    seq = seq.upper()
    ndic = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    revcomp = ''
    for nuc in seq:
        revcomp += ndic[nuc]
    return revcomp[::-1]

def findBindingSites(chrom, cleavage_site, strand, refseq, window=50):
    binding_sites = {'AATAAA':1,'ATTAAA':2,'AGTAAA':3,'TATAAA':4,
                     'CATAAA':5,'GATAAA':6,'AATATA':7,'AATACA':8,
                     'AATAGA':9,'AAAAAG':10,'ACTAAA':11,'AAGAAA':12,
                     'AATGAA':13,'TTTAAA':14,'AAAACA':15,'GGGGCT':16}
    results = []
    if strand == '+':
        try:
            #print 'Attempting to fetch {}:{}-{}'.format(chrom, cleavage_site-window, cleavage_site)
            seq = refseq.fetch(chrom, cleavage_site-window, cleavage_site).upper()
        except (ValueError, IndexError) as e:
            return None
        for i in xrange(len(seq)):
            hexamer = seq[i:i+6]
            if (hexamer in binding_sites):
                results.append([cleavage_site + i - window + 1, binding_sites[hexamer]])
    else:
        try:
            #print 'Attempting to fetch {}:{}-{}'.format(chrom, cleavage_site-window, cleavage_site)
            seq = refseq.fetch(chrom, cleavage_site, cleavage_site+window).upper()
        except (ValueError, IndexError) as e:
            return None
        for i in xrange(len(seq)):
            hexamer = revComp(seq[i:i+6])
            if (hexamer in binding_sites):
                results.append([cleavage_site + i + 6, binding_sites[hexamer]])
    if results:
        results = sorted(results, key=lambda(x):x[1])[0]
        return results
    return None

parser = argparse.ArgumentParser(description='Fix PAS signals for KLEAT runs')

parser.add_argument('kleat', nargs='+', help='KLEAT output file')
parser.add_argument('-r', '--reference_genome', help='Reference Genome in FASTA format. Default is "/home/dmacmillan/references/hg19/hg19.fa"', default='/home/dmacmillan/references/hg19/hg19.fa')
#parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
#parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

#if not os.path.isdir(args.outdir):
#    os.makedirs(args.outdir)

# Global variables
#outfile = os.path.join(args.outdir, args.name)
hg19 = pysam.FastaFile(args.reference_genome)
#kleat_header = 'gene\ttranscript\ttranscript_strand\tcoding\tcontig\tchromosome\tcleavage_site\twithin_UTR\tdistance_from_annotated_site\tESTs\tlength_of_tail_in_contig\tnumber_of_tail_reads\tnumber_of_bridge_reads\tmax_bridge_read_tail_length\tbridge_read_identities\ttail+bridge_reads\tnumber_of_link_pairs\tmax_link_pair_length\tlink_pair_identities\thexamer_loc+id\t3UTR_start_end'

#with open(outfile, 'w') as o:
#    o.write(kleat_header)
#    for kleat in Kleat.parseKleat(args.kleat):
#        kleat.pas = findBindingSites(kleat.chromosome, kleat.cleavage_site, kleat.transcript_strand, hg19)
#        if kleat.pas == []:
#            kleat.pas = '-'
#        o.write(str(kleat) + '\n')

print ('\t').join(Kleat.columns)
for _file in args.kleat:
    for kleat in Kleat.parseKleatFast(_file):
        #try:
        #    old = kleat.pas[:]
        #except TypeError:
        #    old = None
        kleat.pas = findBindingSites(kleat.chromosome, kleat.cleavage_site, kleat.transcript_strand, hg19)
        print kleat
        #old_dist = abs(kleat.cleavage_site - old[0])
        #new_dist = abs(kleat.cleavage_site - kleat.pas[0])
        #diff = abs(new_dist - old_dist)
        #old_str = old[1]
        #new_str = kleat.pas[1]
#        dold = dnew = sold = snew = '-'
#        if old:
#            dold = abs(kleat.cleavage_site - old[0])
#            sold = old[1]
#        if kleat.pas:
#            dnew = abs(kleat.cleavage_site - kleat.pas[0])
#            snew = kleat.pas[1]
#        if old != kleat.pas:
#            print '{}\t{}\t{}\t{}'.format(dold, sold, dnew, snew)
#        #print str(kleat)
