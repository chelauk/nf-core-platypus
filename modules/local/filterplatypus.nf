// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTERPLATYPUS {
    tag "$patient"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(patient), val(id_sample_norm), file(vcf), file(tbi)

    output:
    tuple val(patient), file("*.filtered.vcf"), file("*.removed.vcf"), emit: filtered_vcf

    script:
    """
    filter_platypus.py $vcf ${id_sample_norm}
    """

    stub:
    """
    echo $id_sample_norm > norm.txt
    touch ${patient}.platypus.filtered.vcf
    touch ${patient}.platypus.removed.vcf
    """

}
