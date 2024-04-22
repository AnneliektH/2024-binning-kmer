rule fmg:
    input:
        contigsketch='../results/sourmash_sketches/rocksdb/{csketch}.rocksdb',
        readsketch= '../results/sourmash_sketches/sketch_reads/{rsketch}.k21.zip'
    output:
        fmg_csv='../results/sourmash/fastmultigather/{rsketch}x{csketch}.csv',
    conda:
        "branchwater"
    threads: 12
    shell:"""
    sourmash scripts fastmultigather \
    {input.readsketch} {input.contigsketch} \
    --output {output.fmg_csv} -c {threads} -k 21 -t 1000 -s 100
    """