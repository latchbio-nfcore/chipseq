/*
 * Reformat design file, check validitiy and create IP vs control mappings
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path samplesheet
    val opts

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row, String seq_center) {
    def meta = [:]
    meta.id = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.antibody = row.antibody
    meta.control = row.control

    def rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = rg

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
    }
    return array
}