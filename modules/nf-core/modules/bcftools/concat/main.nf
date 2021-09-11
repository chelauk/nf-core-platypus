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
    tuple val(patient), val(control), path(vcf), path(tbi)

    output:
    tuple val(patient), val(control), path("*concat.vcf.gz"), emit: vcf
    path  "*.version.txt"        , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${patient}${options.suffix}" : "${patient}_concat"
    """
    bcftools concat \\
        --output ${prefix}.vcf.gz \\
        $options.args \\
        --threads $task.cpus \\
        ${vcfs}

    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """

    stub:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${patient}${options.suffix}" : "${patient}_concat"
    """
    touch ${prefix}.vcf
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
