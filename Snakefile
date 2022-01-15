
import pandas
from Bio import SeqIO
try:
    from jme.jupy_tools import concatemers
except:
    from jme.jupy_tools.experimental import concatemers

default_template = "test/{sample}.fasta"
output_dir = config.get('output_dir', 'test/output')
reads_template = config.get('reads_template', default_template)
reads_format = config.get('reads_format', 'fasta')

len_cutoff = config.get('len_cutoff', 5000)
mlen_cutoff = config.get('mlen_cutoff', 500)
pctid_cutoff = config.get('pctid_cutoff', 0)
cl_linkage = config.get('cluster_linkage', 'average') # linkage shouldn't matter much
cl_cutoff = config.get('cluster_cutoff', 750)

methods = ['fft', 'cluster']
method_table_template = "{output_dir}/self.{sample}.cmers.{method}.tsv"

samples, = glob_wildcards(reads_template)
samples = [s for s in samples if re.search('/', s) is None] # no directory breaks in samples
print(f"processing {len(samples)} fasta files")

rule all:
    input:
        expand("{output_dir}/{sample}.cmers.merged.tsv", 
               sample=samples,
               output_dir=[output_dir,],
              ),

rule merge_methods:
    input: lambda w: expand(method_table_template, \
                            sample=[w.sample], \
                            output_dir=[w.output_dir], \
                            method=methods)
    output:
        table="{output_dir}/{sample}.cmers.merged.tsv"
    run:
        merged_table = None
        last_method = None
        for method in methods:
            method_table_file = method_table_template.format(sample=wildcards.sample,
                                                             output_dir=wildcards.output_dir,
                                                             method=method)
            method_table = pandas.read_csv(method_table_file, sep='\t', index_col=0)
            method_table.columns = [f"{c}_{method}" for c in method_table.columns]
            if merged_table is None:
                merged_table = method_table
            else:
                merged_table = merged_table.join(method_table)
                
        merged_table.to_csv(output.table, sep='\t')
            
        
rule filter_by_length:
    input: 
        reads=reads_template
    output: 
        fasta="{output_dir}/{sample}/reads.filtered.fasta"
    benchmark: '{output_dir}/benchmarks/filter_by_len.{sample}.time'
    run:
        kept, dropped = 0, 0
        with open(output.fasta, 'wt') as fasta_handle:
            for read in SeqIO.parse(input.reads, reads_format):
                if len(read) < len_cutoff:
                    dropped += 1
                    continue
                kept += 1
                fasta_handle.write(read.format('fasta'))
        with open(f"{output.fasta}.log", 'wt') as log_out:
            log_out.write(f"kept {kept} reads and dropped {dropped}")


# perl scriptlet to only keep self hits
PERL_SELF_FILTER_CMD = """
if ( $F[0] eq $F[5] ) {
    print;
}
"""

rule minimap:
    """ run minimap and pipe through above perl code to only keep self hits """
    input:
        fasta=rules.filter_by_length.output.fasta,
    output:
        m8="{output_dir}/{sample}/reads.v.self.paf"
    threads: 20
    benchmark: '{output_dir}/benchmarks/minimap.{sample}.time'
    shell:
        """
        minimap2 -t {threads} -x ava-ont {input.fasta} {input.fasta} \
        2> {output}.minimap.log \
          | perl -lane '{PERL_SELF_FILTER_CMD}' \
          > {output.m8}
        """

rule cmer_table:
    input:
        sample_m8="{output_dir}/{sample}/reads.v.self.paf"
    output: 
        sample_cmers=method_table_template
    benchmark: '{output_dir}/benchmarks/cmers.{sample}.{method}.time'
    run:
        concatemers.build_all_v_all_cmer_table(input.sample_m8, \
                                               output.sample_cmers, \
                                               add_full_hit=True, \
                                               qh_cutoff=mlen_cutoff, \
                                               pctid_cutoff=pctid_cutoff, \
                                               table_format="PAF", \
                                               method=wildcards.method)
