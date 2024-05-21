# rule vamb
rule vamb:
    input:
        depthf='../results/depth_files/{manysearch}.h{intersect_hash}.depth.csv',
        fasta = '../resources/ERR11351.fasta'
    output:
        check='../results/vamb/check/{manysearch}.h{intersect_hash}.done.txt',
    conda:
        "vamb"
    threads: 1
    shell:"""
    vamb --fasta {input.fasta} --jgi {input.depthf} \
    --outdir ../results/vamb/{wildcards.manysearch}.h{wildcards.intersect_hash}/ && \
    touch {output.check} --minfasta 50000 -p {threads} 
    """

# ule vamb get the actual bins
rule vamb_get_bins:
    input:
        #check='../results/depth_files/{manysearch}.h{intersect_hash}.depth.csv',
        fasta = '../resources/ERR11351.fasta',
        clusterfile = '../results/vamb/{manysearch}.h{intersect_hash}/clusters.tsv'
    output:
        check='../results/vamb/check/{manysearch}.h{intersect_hash}.getbins.txt',
    conda:
        "vamb_cuda"
    threads: 1
    shell:"""
    python ./scripts/vamb_create_fasta.py \
    {input.fasta} {input.clusterfile} 50000 \
    ../results/vamb/{wildcards.manysearch}.h{wildcards.intersect_hash}/vamb_bins && \
    touch {output.check}
    """

# rule checkm for vamb
rule checkm2_vamb:
    input:
        vamb='../results/vamb/check/{manysearch}.h{intersect_hash}.getbins.txt',
    output:
        check='../results/checkm/check/checkm_{manysearch}.h{intersect_hash}.vamb.done.txt',
    conda:
        "checkm2"
    threads: 8
    shell:"""
    checkm2 predict --threads {threads} \
    --input ../results/vamb/{wildcards.manysearch}.h{wildcards.intersect_hash}/vamb_bins/ \
    -x .fna --output-directory ../results/checkm/vamb/{wildcards.manysearch}.h{wildcards.intersect_hash}/ && \
    touch {output.check}
    """