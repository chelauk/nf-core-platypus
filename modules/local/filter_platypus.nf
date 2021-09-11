// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process FILTER_PLATYPUS {
    label 'process_low'
    tag "$patient"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Platypus', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/platypus-variant%3A0.8.1.2--py27hb763d49_0"
    } else {
        container "quay.io/biocontainers/platypus-variant%3A0.8.1.2--py27hb763d49_0"
    }

    input:
    tuple val(patient), val(id_sample_norm), file(vcf)

    output:
    tuple path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"), emit: filter_platypus_out

    script:
    """
    filter_platypus.py $vcf ${id_sample_norm}
    mv ${patient}_concat_filtered.vcf Platypus_${patient}.filtered.vcf
    bgzip  Platypus_${patient}.filtered.vcf
    tabix -p vcf Platypus_${patient}.filtered.vcf.gz
    """

    stub:
    """
    echo $id_sample_norm > norm.txt
    touch ${prefix}.filtered.vcf.gz
    """
}
