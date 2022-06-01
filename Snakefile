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

methods = ['fft', 'cluster']
searches = ['lastal','minimap']
method_table_template = "{output_dir}/self.{sample}.cmers.{search}.{method}.tsv"

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
                            search=searches, \
                            method=methods)
    output:
        table="{output_dir}/{sample}.cmers.merged.tsv"
    run:
        merged_table = None
        state_cols = []
        for method in methods:
            for search in searches:
                method_table_file = method_table_template.format(sample=wildcards.sample,
                                                                 output_dir=wildcards.output_dir,
                                                                 search=search,
                                                                 method=method)
                method_table = pandas.read_csv(method_table_file, sep='\t', index_col=0)
                other_cols = [c for c in method_table.columns if c != 'qlen']
                new_names = [f"{c}_{search}_{method}" for c in other_cols]
                state_cols.append(f"state_{search}_{method}")
                if merged_table is None:
                    # make length the first column
                    method_table = method_table[['qlen'] + other_cols]
                    # rename columns
                    method_table.columns = ['Length'] + new_names
                    merged_table = method_table
                else:
                    # drop the lenth column, it's redundant
                    method_table = method_table[other_cols]
                    # rename columns
                    method_table.columns = new_names
                    # merge into previous data
                    merged_table = merged_table.join(method_table)
                
        # make a master state col
        state_eval = "state = " + (" and ".join(state_cols))
        merged_table = merged_table.eval(state_eval)
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


# perl scriptlet to only keep self hits (minimap)
PERL_SELF_FILTER_CMD_MM = """
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
    threads: config.get('minimap_threads', 20)
    benchmark: '{output_dir}/benchmarks/minimap.{sample}.time'
    shell:
        """
        minimap2 -t {threads} -x ava-ont {input.fasta} {input.fasta} \
        2> {output.m8}.minimap.err \
          | perl -lane '{PERL_SELF_FILTER_CMD_MM}' \
          > {output.m8} 2> {output.m8}.perl.err
        """

# perl scriptlet to only keep self hits (lastal)
PERL_SELF_FILTER_CMD_LL = """
if ( $F[0] eq $F[1] ) {
    print;
}
"""

rule lastal:
    """ run lastal and pipe through above perl code to only keep self hits """
    input:
        fasta=rules.filter_by_length.output.fasta,
    output:
        m8="{output_dir}/{sample}/reads.v.self.lastn"
    threads: config.get('lastal_threads', 20)
    benchmark: '{output_dir}/benchmarks/minimap.{sample}.time'
    shell:
        """
        lastdb -P {threads} {input.fasta} {input.fasta} > {input.fasta}.lastdb.out 2>&1
        lastal -P {threads} -f blasttab+ {input.fasta} {input.fasta} \
         2> {output.m8}.lastal.err \
          | perl -lane '{PERL_SELF_FILTER_CMD_LL}' \
          > {output.m8} 2> {output.m8}.perl.err
        """

rule cmer_table:
    input:
        sample_m8=lambda w: "{output_dir}/{sample}/reads.v.self.paf"\
                            if w.search == "minimap" \
                            else "{output_dir}/{sample}/reads.v.self.lastn"
    output: 
        sample_cmers=method_table_template
    benchmark: '{output_dir}/benchmarks/cmers.{sample}.{search}.{method}.time'
    run:
        m8_format = "PAF" if wildcards.search == "minimap" else "BLASTTAB+" 
        concatemers.build_all_v_all_cmer_table(input.sample_m8, \
                                               output.sample_cmers, \
                                               add_full_hit=True, \
                                               qh_cutoff=mlen_cutoff, \
                                               pctid_cutoff=pctid_cutoff, \
                                               table_format=m8_format, \
                                               method=wildcards.method)
