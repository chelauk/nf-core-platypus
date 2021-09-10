#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/platypus
========================================================================================
    Github : https://github.com/nf-core/platypus
    Website: https://nf-co.re/platypus
    Slack  : https://nfcore.slack.com/channels/platypus
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai     = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PLATYPUS } from './workflows/platypus'

//
// WORKFLOW: Run main nf-core/platypus analysis pipeline
//
workflow NFCORE_PLATYPUS {
    PLATYPUS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_PLATYPUS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
