# sketch signatures for each multifasta file
rule sketch:
    input:
        fa='../results/sourmash_sketches/split_fa/{sample}.fa',
    output:
        sig='../results/sourmash_sketches/split_sig/{sample}.sig.gz',
    conda:
        "branchwater"
    threads: 1
    shell:"""
    sourmash sketch dna \
    -p k=21,scaled=100 \
    {input.fa} --singleton -o {output.sig}
    """

rule rocksdb:
    input:
        sig='../results/sourmash_sketches/split_sig/{sample}.sig.gz',
    output:
        siglist = '../results/sourmash_sketches/split_sig/{sample}.txt',
    conda:
        "branchwater"
    threads: 1
    shell:"""
    readlink -f {input.sig} > {output.siglist} && \
    sourmash scripts index \
    {output.siglist} -m DNA -k 21 --scaled 100 \
    -o ../results/sourmash_sketches/rocksdb/{wildcards.sample}.rocksdb
    """