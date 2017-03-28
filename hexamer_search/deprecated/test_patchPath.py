def revComp(seq):
    seq = seq.upper()
    ndic = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    revcomp = ''
    for nuc in seq:
        revcomp += ndic[nuc]
    return revcomp[::-1]


binding_sites = {'AATAAA':1,'ATTAAA':2,'AGTAAA':3,'TATAAA':4,
                 'CATAAA':5,'GATAAA':6,'AATATA':7,'AATACA':8,
                 'AATAGA':9,'AAAAAG':10,'ACTAAA':11,'AAGAAA':12,
                 'AATGAA':13,'TTTAAA':14,'AAAACA':15,'GGGGCT':16}


def plus_find(seq, cleavage_site, window=9):
    results = []
    for i in xrange(len(seq)):
        hexamer = seq[i:i+6]
        if (hexamer in binding_sites):
            results.append([cleavage_site + i - window + 1, binding_sites[hexamer]])
    return results


def minus_find(seq, cleavage_site, window=9):
    results = []
    for i in xrange(len(seq)):
        hexamer = revComp(seq[i:i+6])
        if (hexamer in binding_sites):
            results.append([cleavage_site + i + 6, binding_sites[hexamer]])
    return results


print(plus_find('GGGAATAAA', 9))
print(minus_find('TTTATTCCC', 1))
