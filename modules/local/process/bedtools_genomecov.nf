process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam), path(flagstat)
    val opts

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt"), emit: scale_factor
    path "*.version.txt", emit: version

    script:
    prefix = opts.suffix ? "${meta.id}${opts.suffix}" : "${meta.id}"
    pe = meta.single_end ? '' : '-pc'
    extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
    """
    SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        $pe \\
        $extend \\
        | sort -T '.' -k1,1 -k2,2n > ${prefix}.bedGraph

    bedtools --version > bedtools.version.txt
    """
}