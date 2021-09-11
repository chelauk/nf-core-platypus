/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPlatypus.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fasta_fai
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options )

//
// MODULE: platypus module
//
def platypusvariant_options   = modules['platypusvariant']
include { PLATYPUSVARIANT } from '../modules/local/platypusvariant' addParams( options: platypusvariant_options )

//
// MODULE: bcftools module
//

def bcftools_options          =modules['bcftools']
include { BCFTOOLS_CONCAT } from '../modules/nf-core/modules/bcftools/concat/main' addParams( options: bcftools_options)

//
// MODULE: bcftools module
//

def filter_platypus_options          =modules['filter_platypus']
include { FILTER_PLATYPUS } from '../modules/local/filter_platypus' addParams( options: filter_platypus_options)

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope

fasta             = params.fasta             ? Channel.fromPath(params.fasta).collect()             : ch_dummy_file
fasta_fai         = params.fasta_fai         ? Channel.fromPath(params.fasta_fai).collect()         : ch_dummy_file

// Initialise input sample
csv_file = file(params.input)
input_samples  = extract_csv(csv_file)
ch_fasta_fai = Channel.from("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                            "chr20","chr21","chr22","chrX")
platypus_input = make_platypus_input(input_samples)
platypus_input = platypus_input.combine(ch_fasta_fai)
//platypus_input.view()

//
// input channel functions
//

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        // create control id necessary for filter platypus
        if (row.bam_c) {
            meta.control = file(row.bam_c).getName()
            meta.control = meta.control.toString().minus('.recal.bam')
        }

        // mapping with platypus
        if (row.vcf && row.bam_t) {
            def vcf         = file(row.vcf, checkIfExists: true)
            def vcf_tbi         = file(row.vcf_tbi, checkIfExists: true)
            def bam_t       = file(row.bam_t, checkIfExists: true)
            def bam_t_bai   = file(row.bam_t_bai, checkIfExists: true)
            def bam_c       = file(row.bam_c, checkIfExists: true)
            def bam_c_bai   = file(row.bam_c_bai, checkIfExists: true)
            return [meta, [vcf, vcf_tbi, bam_t, bam_t_bai, bam_c, bam_c_bai]]
        // recalibration
        }
    }
}

def make_platypus_input(input) {
    return input
        .map { meta, files -> [ meta.patient, meta.control,
        [files[0],files[1]],[files[2],files[3],files[4],files[5]]]}
        .groupTuple()
        .map { patient, control, vcfs, bams  -> [ patient,control,vcfs.flatten(),bams.flatten()] }
        .map { patient, control, vcfs, bams  -> [ patient,control.unique().join(""),vcfs,bams.unique()] }
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PLATYPUS {

    ch_software_versions = Channel.empty()

    //
    // MODULE: Run platypus
    //

    PLATYPUSVARIANT(platypus_input, fasta, fasta_fai)
    BCFTOOLS_CONCAT(PLATYPUSVARIANT.out.platypus_vcf.groupTuple())
    filter_vcf_in = BCFTOOLS_CONCAT.out.vcf
    filter_vcf_in = filter_vcf_in
                        .map{patient, control, vcf -> [patient,control.unique().join(""),vcf]}
    FILTER_PLATYPUS(filter_vcf_in)


    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPlatypus.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
