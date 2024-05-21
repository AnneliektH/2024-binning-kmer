# rule metabat
rule metabat:
    input:
        depthf='../results/depth_files/{manysearch}.h{intersect_hash}.95.depth.csv',
        fasta = '../resources/ERR11351.fasta'
    output:
        check='../results/metabat2/check/{manysearch}.h{intersect_hash}.95.done.txt',
    conda:
        "metabat2"
    threads: 8
    shell:"""
    metabat2 -m 1500 \
    -i {input.fasta} -a {input.depthf} -t {threads} \
    -o ../results/metabat2/{wildcards.manysearch}.h{wildcards.intersect_hash}.95/{wildcards.manysearch}.h{wildcards.intersect_hash}.95 && \
    touch {output.check}
    """

# rule checkm
rule checkm2_mb2:
    input:
        metabat ='../results/metabat2/check/{manysearch}.h{intersect_hash}.95.done.txt'
    output:
        check='../results/checkm/checkm_{manysearch}.h{intersect_hash}.mb2.95.done.txt',
    conda:
        "checkm2"
    threads: 8
    shell:"""
    checkm2 predict --threads {threads} \
    --input ../results/metabat2/{wildcards.manysearch}.h{wildcards.intersect_hash}.95/ \
    -x .fa --output-directory ../results/checkm/metabat2/{wildcards.manysearch}.h{wildcards.intersect_hash}.95/ && \
    touch {output.check}
    """