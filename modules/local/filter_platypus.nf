// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process FILTER_PLATYPUS {
    tag "$patient"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Platypus', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(patient), val(id_sample_norm), file(vcf)

    output:
    path '*.filtered.vcf.gz'

    script:
    prefix       = options.suffix ? "${patient}${options.suffix}" : "platypus_${patient}"
    """
    filter_platypus.py ${prefix}.vcf ${id_sample_norm}
    bgzip  ${prefix}_filtered.vcf
    tabix -p vcf ${prefix}.filtered.vcf.gz
    """

    stub:
    prefix       = options.suffix ? "${patient}${options.suffix}" : "platypus_${patient}"
    """
    echo $id_sample_norm > norm.txt
    touch ${prefix}.filtered.vcf.gz
    """
}
