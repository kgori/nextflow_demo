params.outputDir = ""
params.from_reference = ""
params.to_reference = ""
params.inputDir = ""

process index_reference {
    input:
    path(reference)

    output:
    path("${reference}.*")

    script:
    """
    bwa index ${reference}
    samtools faidx ${reference}
    """
}

process cram_to_fastq {
    input:
    tuple val(id), path(cramfile)

    output:
    tuple val(id), path("${id}_1.fastq.gz"), path("${id}_2.fastq.gz")

    script:
    """
    samtools fastq -1 ${id}_1.fastq.gz -2 ${id}_2.fastq.gz -0 /dev/null -s /dev/null -n ${cramfile[0]}
    """
}

process align_fastq {
    input:
    tuple val(id), path(fastq1), path(fastq2), path(reference)

    output:
    path("${id}.bam")

    script:
    """
    bwa mem -Y ${reference} ${fastq1} ${fastq2} |\
      samtools sort |\
      samtools view -O CRAM -T ${reference} -o ${id}.bam
    """
}

process merge_alignments {
    input:
    tuple val(name), path(bamfiles)

    output:
    path("${name}.bam*")

    script:
    """
    samtools merge ${name}.bam ${bamfiles}
    samtools index ${name}.bam
    """
}

workflow {
    println("Running workflow")
    // Make new channels from input data
    ref_from = Channel.fromPath("${params.from_reference}", checkIfExists: true)
    ref_to = Channel.fromPath("${params.to_reference}", checkIfExists: true)
    crams = Channel.fromFilePairs("${params.inputDir}/*.cram{,.crai}", checkIfExists: true)

    // Index the reference
    ref_index = index_reference(ref_from.concat(ref_to))

    ref_index | view
}

