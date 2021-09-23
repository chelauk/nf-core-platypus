// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLATYPUSVARIANT {

    tag "$patient"
    label 'eight_cpus'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    conda (params.enable_conda ? "bioconda::platypus-variant=0.8.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/platypus-variant%3A0.8.1.2--py27hb763d49_0"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(patient), val(control), file(vcf), file(bam), val(chr)
    path fasta
	path fasta_fai

    output:
    tuple val(patient), val(control), path("*_platypus.vcf.gz"), path("*_platypus.vcf.gz.tbi"), emit: platypus_vcf
    path "*.version.txt"                                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${patient}${options.suffix}" : "${patient}"
    """
    bams=(*bam)
    bams=\${bams[@]}
    bams=\${bams// /,}
    vcfs=(*vcf.gz)
    vcfs=\${vcfs[@]}
    vcfs=\${vcfs// /,}
    platypus callVariants \\
    --bamFiles=\$bams \\
    --refFile=$fasta \\
    --regions=$chr \\
    --output=${chr}_${patient}_platypus.vcf \\
    --source=\$vcfs \\
    --filterReadPairsWithSmallInserts=0 \\
    --maxReads=100000000 \\
    --minPosterior=0 \\
    --nCPU=${task.cpus} \\
    --getVariantsFromBAMs=0

    bgzip ${chr}_${patient}_platypus.vcf
    tabix -p vcf ${chr}_${patient}_platypus.vcf.gz
    echo "platypus.0.2.8.2" > platypus.version.txt
    """

    stub:
    """
    bams=(*bam)
    bams=\${bams[@]}
    bams=\${bams// /,}
    vcfs=(*vcf.gz)
    vcfs=\${vcfs[@]}
    vcfs=\${vcfs// /,}
    touch ${chr}_${patient}_platypus.vcf.gz
    touch ${chr}_${patient}_platypus.vcf.gz.tbi
    echo \$bams > bam_files
    echo \$vcfs > vcf_files
    touch platypus.version.txt
    """
}
