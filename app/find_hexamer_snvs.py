#!/usr/bin/env python
import argparse
import pysam
import sys
from pybedtools import BedTool
import os

def reverse_complement(seq):
    """Reverse complements sequence string"""
    from string import maketrans
    complement = maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)

def filter_vcf(vcf, min_dp, min_alt_fraction, out_file=None):
    """Extract SNVs based on minimum depth and minimum alternative allele fraction"""
    bcf_out = None
    if out_file is not None:
        bcf_out = pysam.VariantFile(out_file, 'w', header=vcf.header)

    filtered = []
    for variant in vcf.fetch():
        if 'INDEL' in variant.info.keys():
            continue
        
        if len(variant.alts) > 1:
            continue
        
        allele_depth = int(variant.info['DP4'][2]) + int(variant.info['DP4'][3])
        
        if int(variant.info['DP']) >= min_dp and float(allele_depth) / float(variant.info['DP']) >= min_alt_fraction:
            allele = variant.alts[0]
            if variant.samples[0]['GT'][0] == variant.samples[0]['GT'][1]:
                zygosity = 'hom'
            else:
                zygosity = 'het'
            filtered.append({'variant': variant,
                             'allele': allele,
                             'allele_depth': allele_depth,
                             'zygosity': zygosity,
                             })
            if bcf_out is not None:
                bcf_out.write(variant)
    
    return filtered

def report(events, outfile):
    header = ('chrom',
              'pos',
              'ref',
              'alt',
              'DP',
              'AD',
              'AD/DP',
              'zygosity',
              'gene',
              'transcript',
              'transcript_strand',
              '3UTR',
              'A_to_G',
              'hexamer_pos',
              'hexamer_effect',
              'hexamer_change',
              'hexamer_to_cleavage',
              'canonical_hexamer'
              )
    out = open(outfile, 'w')
    out.write('\t'.join(header) + '\n')
    for event in events:
        data = []
        data.append(event['variant'].chrom)
        data.append(event['variant'].pos)
        data.append(event['variant'].ref)
        data.append(event['allele'])
        data.append(event['variant'].info['DP'])
        data.append(event['allele_depth'])
        data.append(event['allele_depth']*100/event['variant'].info['DP'])
        data.append(event['zygosity'])
        data.append(event['gene'])
        data.append(event['transcript'])
        data.append(event['transcript_strand'])
        data.append(event['3UTR'])
        data.append(event['A_to_G'])
        data.append(event['hexamer_pos'])
        data.append(event['hexamer_effect'])
        data.append(event['hexamer_change'])
        data.append(event['hexamer_to_cleavage'])
        data.append(event['canonical_hexamer'])
        out.write('\t'.join(map(str, data)) + '\n')
    out.close()
    
def associate_utr(filtered_vcf, utr_gff, events):
    """associate filtered SNVs to UTR"""
    filtered_bed = BedTool(filtered_vcf)
    utr_bed = BedTool(utr_gff)
    
    utr_olapped = {}
    for olap in filtered_bed.intersect(utr_bed, wb=True).saveas('tmp.bed'):
        utr_start = int(olap[13])
        utr_end = int(olap[14])
        transcript_strand = olap[16]
        transcript, gene = olap[18].split('::')
        
        pos = '%s:%s' % (olap[0], olap[1])
        utr_olapped[pos] = gene, transcript, transcript_strand, utr_start, utr_end
        
    remove = []
    for i in range(len(events)):
        pos = '%s:%s' % (events[i]['variant'].chrom, events[i]['variant'].pos)
        if not utr_olapped.has_key(pos):
            remove.append(i)
        else:
            events[i]['gene'] = utr_olapped[pos][0]
            events[i]['transcript'] = utr_olapped[pos][1]
            events[i]['transcript_strand'] = utr_olapped[pos][2]
            events[i]['3UTR'] = '%s:%s-%s' % (events[i]['variant'].chrom, utr_olapped[pos][3], utr_olapped[pos][4])
            events[i]['utr_start'] = int(utr_olapped[pos][3])
            events[i]['utr_end'] = int(utr_olapped[pos][4])
            
            if (events[i]['transcript_strand'] == '+' and events[i]['variant'].ref == 'A' and events[i]['allele'] == 'G') or\
               (events[i]['transcript_strand'] == '-' and events[i]['variant'].ref == 'T' and events[i]['allele'] == 'C') :
                events[i]['A_to_G'] = True
            else:
                events[i]['A_to_G'] = False                
            
    for i in reversed(remove):
        del events[i]
        
def scan_for_hexamer(fasta, chrom, pos, ref, alt, motifs_plus, motifs_minus, ranks):
    """Look for possible hexamers at given chromosome position"""
    motifs = motifs_plus + motifs_minus
    hexamer_start = None
    effect = None
    hexamer_change = None
    
    result = None
    #chrom_ucsc = 'chr' + chrom
    for base in range(pos - 5, pos + 1):
        hexamer_ref = fasta.fetch(chrom, base - 1, base + 5)
        hexamer_idx = pos - base
        hexamer_alt_list = list(hexamer_ref)
        if not hexamer_ref:
            continue
        hexamer_alt_list[hexamer_idx] = alt
        hexamer_alt = ''.join(hexamer_alt_list)
        
        hexamer_change = (hexamer_ref, hexamer_alt)
        
        if (hexamer_ref.upper() in motifs_plus and\
            hexamer_alt.upper() not in motifs_plus) or\
           (hexamer_ref.upper() in motifs_minus and\
            hexamer_alt.upper() not in motifs_minus):
            if hexamer_ref.upper() in motifs_plus:
                strand = '+'
                hexamer_end = base + 5
            else:
                strand = '-'
                hexamer_end = base
            
            result = {'chrom': chrom,
                      'start': pos,
                      'ref': '%s(%s)' % (hexamer_ref.upper(), ranks[hexamer_ref.upper()]),
                      'alt':  hexamer_alt.upper(),
                      'strand': strand,
                      'effect': 'destroy',
                      'hexamer_end': hexamer_end
                      }
            break
        elif (hexamer_ref.upper() not in motifs_plus and\
              hexamer_alt.upper() in motifs_plus) or\
             (hexamer_ref.upper() not in motifs_minus and\
              hexamer_alt.upper() in motifs_minus):
            if hexamer_alt.upper() in motifs_plus:
                strand = '+'
                hexamer_end = base + 5
            else:
                strand = '-'
                hexamer_end = base
            result = {'chrom': chrom,
                      'start': pos,
                      'ref': hexamer_ref.upper(),
                      'alt': '%s(%s)' % (hexamer_alt, ranks[hexamer_alt.upper()]),
                      'strand': strand,
                      'effect': 'create',
                      'hexamer_end': hexamer_end
                      }
            break
        elif (hexamer_ref.upper() in motifs_plus and\
              hexamer_alt.upper() in motifs_plus) or\
             (hexamer_ref.upper() in motifs_minus and\
              hexamer_alt.upper() in motifs_minus):
            if hexamer_ref.upper() in motifs_plus:
                strand = '+'
                hexamer_end = base + 5
            else:
                strand = '-'
                hexamer_end = base
            
            if ranks[hexamer_ref.upper()] < ranks[hexamer_alt.upper()]:
                effect = 'shift-'
            else:
                effect = 'shift+'
            
            result = {'chrom': chrom,
                      'start': pos,
                      'ref': '%s(%s)' % (hexamer_ref.upper(), ranks[hexamer_ref.upper()]),
                      'alt': '%s(%s)' % (hexamer_alt.upper(), ranks[hexamer_alt.upper()]),
                      'strand': strand,
                      'effect': effect,
                      'hexamer_end': hexamer_end
                      }
            break
        
    return result

def check_hexamer(events, fasta, distance_range):
    """screen events for possible hexamer effects"""
    motifs_plus, motifs_minus, ranks = get_motifs()

    remove = []
    for i in range(len(events)):
        event = events[i]
        if event['allele'] is None:
            remove.append(i)
            continue

        hexamer = scan_for_hexamer(fasta,
                                   event['variant'].chrom,
                                   event['variant'].pos,
                                   event['variant'].ref,
                                   event['allele'],
                                   motifs_plus,
                                   motifs_minus,
                                   ranks,
                                   )
        if hexamer is None:
            remove.append(i)
            continue
        
        if hexamer['strand'] != event['transcript_strand']:
            remove.append(i)
            continue
        
        chrom = event['variant'].chrom
        if hexamer['strand'] == '+':
            hexamer_to_cleavage = int(event['utr_end']) - int(hexamer['hexamer_end'])
            hexamer_pos = '%s:%d-%d' % (chrom, hexamer['hexamer_end'] - 5, hexamer['hexamer_end'])
        else:
            hexamer_to_cleavage = int(hexamer['hexamer_end']) - int(event['utr_start'])
            hexamer_pos = '%s:%d-%d' % (chrom, hexamer['hexamer_end'], hexamer['hexamer_end'] + 5)
            
        if hexamer_to_cleavage < distance_range[0] or hexamer_to_cleavage > distance_range[1]:
            remove.append(i)
            continue
        
        event['hexamer_to_cleavage'] = hexamer_to_cleavage
        event['hexamer_pos'] = hexamer_pos
        event['hexamer_effect'] = hexamer['effect']
        event['hexamer_change'] = '%s->%s' % (hexamer['ref'], hexamer['alt'])
        event['canonical_hexamer'] = find_canonical_hexamer(chrom,
                                                            event['utr_end'] if event['transcript_strand'] == '+' else  event['utr_start'],
                                                            hexamer['strand'],
                                                            distance_range,
                                                            motifs_plus,
                                                            motifs_minus,
                                                            ranks,
                                                            fasta)
            
    for i in reversed(remove):
        del events[i]
        
def get_motifs():
    """put hexamer motifs and their ranks into data structure"""
    motifs_plus = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA',
                   'AATATA', 'AATACA', 'AATAGA', 'AAAAAG', 'ACTAAA']
    motifs_minus = []
    ranks = {}
    rank = 1
    for motif in motifs_plus:
        motifs_minus.append(reverse_complement(motif))
        ranks[motif] = rank
        rank += 1
    rank = 1
    for motif in motifs_minus:
        ranks[motif] = rank
        rank += 1

    return motifs_plus, motifs_minus, ranks
        
def find_canonical_hexamer(chrom, utr_end, strand, distance_range, motifs_plus, motifs_minus, ranks, fasta):
    def search(window, motifs, found):
        for pos in range(window[0] - 5, window[1] + 5 + 1):
            hexamer = fasta.fetch(ucsc_chrom, pos - 1, pos + 5)
            if hexamer.upper() in motifs:
                coord = '%s:%d-%d' % (chrom, pos, pos + 5)
                found[coord] = hexamer.upper()

    ucsc_chrom = chrom
    #if not 'chr' in ucsc_chrom:
        #ucsc_chrom = 'chr' + chrom
    found = {}
    if strand == '+':
        window = utr_end - distance_range[1], utr_end - distance_range[0]
        search(window, motifs_plus, found)
    elif strand == '-':
        window = utr_end + distance_range[0], utr_end + distance_range[1]
        search(window, motifs_minus, found)

    if found:
        motifs_sorted = sorted(found.keys(), key = lambda m: ranks[found[m]])
        return motifs_sorted[0]
        #return '%s(%s)' % (motifs_sorted[0], ranks[found[motifs_sorted[0]]])
    else:
        return '-'
    
    
def parse_args():
    parser = argparse.ArgumentParser(description='Find SNVs from samtools mpileup that potentially affect hexamers(PAS)')
    parser.add_argument("vcf", type=str, help="samtools mpileup vcf")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("utr3", type=str, help="3utr gff")
    parser.add_argument("outfile", type=str, help="output file")
    parser.add_argument("--min_dp", type=int, default=5, help="min depth Default:5")
    parser.add_argument("--min_alt_fraction", type=float, help="minimum percentage of alternatvie allele. Default=0.5", default=0.5)
    parser.add_argument("--distance", nargs=2, type=int, default=(10,40),\
                        help="distance range between hexamer and cleavage site. Default: 10,40")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    events = []
    if os.path.getsize(args.vcf) > 0:
        input_vcf = pysam.VariantFile(args.vcf)

        filtered_vcf = '%s/%s_filtered.vcf' % (os.path.dirname(os.path.abspath(args.vcf)),
                                               os.path.splitext(os.path.basename(args.vcf))[0])
        fasta = pysam.FastaFile(args.genome_fasta)
    
        events = filter_vcf(input_vcf, args.min_dp, args.min_alt_fraction, out_file=filtered_vcf)
        print 'filtered:%s' % len(events)

        associate_utr(filtered_vcf, args.utr3, events)
        print 'utr:%s' % len(events)
        
        check_hexamer(events, fasta, args.distance)
        print 'hexamer:%s' % len(events)
        
    report(events, args.outfile)
    
main()
