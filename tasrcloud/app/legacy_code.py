

@R.mkdir(sort_bam_by_name, R.formatter(), '{subpath[0][1]}/bam2fastq')
@R.transform(sort_bam_by_name, R.formatter(), [
    '{subpath[0][1]}/bam2fastq/a_1.fq',
    '{subpath[0][1]}/bam2fastq/a_2.fq',
    '{subpath[0][1]}/bam2fastq/bam2fastq.log',
    '{subpath[0][1]}/bam2fastq/bam2fastq.SUCCESS'
])
@U.timeit
def bam2fastq(inputs, outputs):
    input_bam, _, _ = inputs
    fq1, fq2, log, flag = outputs
    cmd = ('bedtools bamtofastq -i {input_bam} -fq {fq1} -fq2 {fq2} 2>&1 '
           '| tee {log}').format(**locals())
    U.execute(cmd, flag)

# instead of converting bam to fastq, which takes a lot of time, just remove
# the @CO lines from the bam, which causes problematic parsing in
# biobloomcategorizer, the same problem as described in
# https://groups.google.com/forum/#!msg/abyss-users/uDoJjgPWeu4/fJ-SYGN8XLsJ
@R.mkdir(sort_bam_by_name, R.formatter(), '{subpath[0][1]}/remove_CO_header')
@R.transform(sort_bam_by_name, R.formatter(), [
    '{subpath[0][1]}/remove_CO_header/a.bam',
    '{subpath[0][1]}/remove_CO_header/remove_CO_header.log',
    '{subpath[0][1]}/remove_CO_header/remove_CO_header.SUCCESS'
])
@U.timeit
def remove_CO_header(inputs, outputs):
    input_bam, _, _ = inputs
    output_bam, log, flag = outputs
    num_cpus = CONFIG['num_cpus']
    cmd = ('samtools view -h -@ {num_cpus} {input_bam} '
           '| grep -v "^@CO" '
           '| samtools view -Sb -@ {num_cpus} - > {output_bam} 2>&1'
           '| tee {log}').format(**locals())
    U.execute(cmd, flag)
