// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_CONCAT {
    tag "$patient"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bcftools=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
    }

    input:
    tuple val(patient), val(control), path(vcfs), path(tbi)

    output:
    tuple val(patient), val(control), path("*.platypus_unfiltered.vcf.gz"), path("*.platypus_unfiltered.vcf.gz.tbi"), emit: vcf
    path  "*.version.txt"        , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${patient}${options.suffix}" : "${patient}.platypus_unfiltered"
    """
    bcftools concat \\
        --output ${prefix}.vcf \\
        $options.args \\
        --threads $task.cpus \\
        ${vcfs}
    bgzip ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """

    stub:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${patient}${options.suffix}" : "${patient}.platypus_unfiltered"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
