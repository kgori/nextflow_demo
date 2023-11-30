process depth {
    input:
    path ref
    path bam

    output:
    path "depth.txt"

    publishDir "${params.outputDir}/depth", mode: 'copy'

    script:
    """
    depth.sh ${ref} ${bam[0]} > depth.txt
    """
}
