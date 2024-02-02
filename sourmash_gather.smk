# fastgather MAGS : MAKING A PICKLIST WORKS BUT YOU NEED A MANIFEST 
# FOR THE INPUT DB

## Maybe chance treshold bp --threshold-bp
rule sketch:
    input:
       fasta='{sample}_contigs.fasta',
       read1='atlas_{sample}/{sample}/sequence_quality_control/{sample}_QC_R1.fastq.gz',
    output:
        fasta_sketch='sourmash/sketch/{sample}.contigs.sig.gz',
        read_sketch='sourmash/sketch/{sample}.read1.sig.gz',
        fasta_zip = 'sourmash/sketch/{sample}.zip',
        fasta_sql = 'sourmash/sketch/{sample}.sqlmf'
    conda: 
        "sourmash"
    shell: """
        sourmash sketch dna \
        -p k=21,scaled=1000,k=31,scaled=1000,k=51,scaled=1000 \
        {input.fasta} --singleton -o {output.fasta_sketch} && \
        sourmash sketch dna \
        -p k=21,scaled=1000,k=31,scaled=1000,k=51,scaled=1000 \
        {input.read1} --singleton -o {output.read_sketch} && \
        sourmash sig cat {output.fasta_sketch} -o {output.fasta_zip} && \
        sourmash sig collect {output.fasta_sketch} -o {output.fasta_sql}
    """

rule fastgather:
    input: 
        fasta_zip='sourmash/sketch/{sample}.zip',
        read_sketch='sourmash/sketch/{sample}.read1.sig.gz'
    output: 
        csv = "sourmash/fgather/{sample}.{ksize}.csv"
    resources:
        # limit to one fastgather with --resources rayon_exclude=1
        rayon_exclude=1,
        mem_mb=15000
    conda: 
        "branchwater"
    shell:
        """
        sourmash scripts fastgather {input.read_sketch} \
        {input.fasta_zip} -k {wildcards.ksize} \
        --scaled 1000 -o {output.csv} -c 1 --threshold-bp=0
        """
rule gather:
    input: 
        fasta_zip='sourmash/sketch/{sample}.zip',
        picklist='sourmash/fgather/{sample}.{ksize}.csv',
        read_sketch='sourmash/sketch/{sample}.read1.sig.gz'
    output: 
        csv = "sourmash/gather/{sample}.{ksize}.csv"
    resources:
        # limit to one fastgather with --resources rayon_exclude=1
        mem_mb=15000
    conda: 
        "sourmash"
    shell:
        """
        sourmash gather {input.read_sketch} {input.fasta_zip}
        -k 21 --scaled 1000 -o {output.csv} \
        --picklist {input.picklist}:match_md5:md5 --threshold-bp=0
        """
