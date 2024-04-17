rule manifest_smashm:
    input:
        fmg='mgmanysearch/{sample}.gather.csv',
    output:
        mf='mgmanysearch/{sample}.mf.csv',
    shell:"""
    sourmash sig check -k 21 {DBCONTIG} --picklist {input.fmg}:match_md5:md5 -m {output.mf}
    """

# Do mgmanysearch
rule do_mgsearch_m:
    input:
        sample='mgmanysearch/sketch_reads/{sample}.sig.gz',
        mf='mgmanysearch/{sample}.mf.csv',
    output:
        csv='mgmanysearch/{sample}.mgm.csv'
    shell: """
    sourmash scripts mgmanysearch --queries {input.mf} \
    --against {input.sample} -k 21 --scaled 1000 -o {output.csv} 
    """