// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BGZIP {
    tag "$patient"
    label 'process_low'
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
    tuple val(patient), path(filtered_vcf), path(removed_vcf)

    output:
    tuple val(patient), path("*.platypus.filtered.vcf.gz"), path("*.platypus.filtered.vcf.gz.tbi"), path("*.platypus.removed.vcf.gz"), path("*.platypus.removed.vcf.gz.tbi"),emit: vcf
    path  "*.version.txt"        , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    bgzip $filtered_vcf
    tabix -p vcf ${filtered_vcf}.gz
    bgzip $removed_vcf
    tabix -p vcf ${removed_vcf}.gz
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """

    stub:
    def software = getSoftwareName(task.process)
    """
    touch ${filtered_vcf}.gz
    touch ${filtered_vcf}.gz.tbi
    touch ${removed_vcf}.gz
    touch ${removed_vcf}.gz.tbi
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
