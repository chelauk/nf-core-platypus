// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLATYPUSVARIANT {

    tag "$patient"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    conda (params.enable_conda ? "bioconda::platypus-variant=0.8.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(patient), val(control), file(vcf), file(bam), val(chr)
    path fasta

    output:
    tuple val(patient), val(control), path("*.platypus.vcf"), emit: platypus_vcf
    path "*.version.txt"                                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${patient}${options.suffix}" : "${patient}"
    """
    platypus callVariants \\
    --bamFiles=${bam.join(',')} \\
    --refFile=$fasta \\
    --regions=$chr \\
    --output=${chr}_${patient}_platypus.vcf
    --source=${vcf.join(',')} \\
    --filterReadPairsWithSmallInserts=0 \\
    --maxReads=100000000 \\
    --minPosterior=0 \\
    --nCPU=${task.cpus}
    """

    stub:
    """
    touch ${patient}_${chr}.platypus.vcf
    echo ${bam.join(',')} > vcf_files
    echo ${vcf.join(',')} > vcf_files
    touch platypus.version.txt
    """
}
