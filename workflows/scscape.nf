/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'
include { validateParameters; paramsHelp; samplesheetToList } from 'plugin/nf-schema'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowScscape.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK          } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QC } from '../subworkflows/local/qc/'
include { POST_QC } from '../subworkflows/local/post_qc'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SUBSET             } from '../modules/local/subset'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
def validation = []

workflow SCSCAPE {

    ch_versions = Channel.empty()
    ch_validation_log = Channel.empty()

    def lhm = new LinkedHashMap()
    ch_contrasts_file = Channel.from(file(params.input))
    ch_contrasts_file.splitJson(path: 'meta')
    .map{ it -> if (it instanceof String){
        lhm[it] = ''
        return lhm
    } }
    .unique()
    .set{ ch_meta_init }


    ch_contrasts_file.splitJson(path: 'samples')
    .flatMap()
    .map { it -> it['value'] }
    .flatMap()
    .buffer( size: 4 )
    .set{ ch_samples_init }

    ch_meta_init.combine(ch_samples_init)
    .map { headers, id, data, groups, content ->

        def dataValues = data.toString().replaceAll("data=\\[|\\]", "").split(",")
        def groupsValues = groups.toString().replaceAll("groups=\\[|\\]", "").split(",")

        def metaStr = content.toString()
        def metaMatch = metaStr =~ /meta_content=\[(.*)\]/
        def metaContentValues = []
        if (metaMatch.find()) {
            metaContentValues = metaMatch[0][1].split(",").collect { it.trim() }
        }

        def meta = new LinkedHashMap()

        def cleanHeaders = headers.collect { it.toString().replace("=", "").trim() }

        cleanHeaders.eachWithIndex { header, i ->
            if (i < metaContentValues.size()) {
                meta[header] = metaContentValues[i].trim()
            }
        }
        groupsValues.eachWithIndex { group, i ->
            groupsValues[i] =  group.trim()
        }
        meta['groups'] = groupsValues.sort()

        return [meta.id, meta, dataValues]
    }.set{ ch_meta_samples }

    if (params.subset) {

        SUBSET(
            ch_meta_samples.map { id, meta, data -> [ meta, data[0], data[1] ] }
        )

        // ch_subset_samples map {meta, data -> [meta.id, data])
        //POST_QC(ch_subset)

    } else {

        QC(
            ch_meta_samples,
            ch_validation_log
        )

        POST_QC(
            QC.out.doublet_filtered_rds,
            QC.out.validation_log
        )

    }

    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    //
    // MODULE: MultiQC
    //
    //workflow_summary    = WorkflowScscape.paramsSummaryMultiqc(workflow, summary_params)
    //ch_workflow_summary = Channel.value(workflow_summary)

    //methods_description    = WorkflowScscape.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    //ch_methods_description = Channel.value(methods_description)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //MULTIQC (
        //ch_multiqc_files.collect(),
        //ch_multiqc_config.toList(),
        //ch_multiqc_custom_config.toList(),
        //ch_multiqc_logo.toList()
    //)
    //multiqc_report = MULTIQC.out.report.toList()


    // ch_validation_logs.map { it -> it.text }
    //    .collectFile(name: "agg_validation.log", newline: true)
    //    .set { ch_validation_log }

    validation = ch_validation_log.collectFile( name: "val.log" ).map{ it -> it.text }.toList()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log, validation)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
