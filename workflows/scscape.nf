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
include { MAKE_SEURAT          } from '../modules/local/makeseurat.nf'
include { NORMALIZE_QC         } from '../modules/local/normalize_qc.nf'
include { FIND_DOUBLETS        } from '../modules/local/doubletfinder.nf'
include { MERGE_SO             } from '../modules/local/merge.nf'
include { MERGE_SO as SCALE_SO } from '../modules/local/merge.nf'
include { RUN_PCA as PCA_MULT  } from '../modules/local/runpca.nf'
include { RUN_PCA as PCA_SING  } from '../modules/local/runpca.nf'
include { INTEGRATION          } from '../modules/local/integration.nf'
include { FIND_NN_CLUSTER      } from '../modules/local/find_NN_clusters.nf'
include { DISPLAY_REDUCTION    } from '../modules/local/plotting.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCSCAPE {

    ch_versions = Channel.empty()

    ch_samples = Channel.fromList(samplesheetToList(params.sample_sheet, "./assets/schema_input.json"))
    ch_contrasts_file = Channel.from(file(params.segmentation_sheet))
    ch_contrasts_file.splitCsv ( header:true, sep:(params.segmentation_sheet.endsWith('tsv') ? '\t' : ','))
                    .flatMap().filter { !(it.toString().toUpperCase().contains("FALSE")) }
                    .map { it ->
                        if (it.toString().contains("id")){
                            lhm = new LinkedHashMap()
                            lhm['id'] = it.toString().split("=")[1]
                            return lhm
                        } else {
                            return it.toString().split("=")[0]
                        }
                    }
                    .collect().map { it.reverse() }.flatMap()
                    .buffer { it instanceof LinkedHashMap }
                    .map { it.reverse() }
                    .set { ch_contrasts }

    ch_contrasts.join(ch_samples).flatMap()
                .map { it ->
                if ( it instanceof LinkedHashMap ){
                    group_ls = new ArrayList()
                    return it + [ groups: group_ls ]
                } else if (it instanceof java.lang.String ){
                    group_ls.add(it)
                    return
                } else if (it instanceof java.util.ArrayList ){
                    return it
                } else {
                    return it
                }
                }
                .map { it ->
                    if (it instanceof LinkedHashMap){
                        return [id: it.id, groups: it.groups.sort()]
                    } else {
                        return it
                    }
                }.collect(flat: false).map { it.reverse() }.flatMap()
                .buffer { it instanceof LinkedHashMap }
                .map { it.reverse() }
                .set { ch_updated_meta }

    //ch_updated_meta.view()
    ch_init_rds = MAKE_SEURAT (
        ch_updated_meta.map { [it[0], it[1]] },
        ch_updated_meta.map { [it[0], it[2]] },
        params.min_cells,
        params.min_features,
        params.gene_identifier
    )
    ch_init_rds.rds.join(ch_updated_meta).set { ch_init_rds }

    ch_normalized_qc = NORMALIZE_QC (
        ch_init_rds.map { [it[0], it[1]] },
        ch_init_rds.map { [it[0], it[3]] },
        params.nfeature_lower,
        params.nfeature_upper,
        params.ncount_lower,
        params.ncount_upper,
        params.max_mito_pct,
        params.vars_2_regress
    )
    ch_normalized_qc.rds.join(ch_updated_meta).set { ch_normalized_qc }

    ch_doublet_filtered_rds = FIND_DOUBLETS (
        ch_normalized_qc.map { [it[0], it[1]] },
        ch_normalized_qc.map { [it[0], it[2]] },
        params.vars_2_regress
    )

    ch_doublet_filtered_rds.rds.map { meta, rds -> [ rds, meta.groups ]}
                    .transpose()
                    .groupTuple(by: 1)
                    .map { data, group -> [ group, data ] }
                    .set { ch_merged_groups }

    ch_merged_groups.branch { group, data ->
                            single: data.size() == 1
                            multiple: data.size() > 1
                            none: data.size() < 1
                            }
                            .set { ch_merged_groups }

    ch_merged_so = MERGE_SO (
        ch_merged_groups.multiple,
        params.vars_2_regress
    )

    ch_scaled_so = SCALE_SO (
        ch_merged_groups.single,
        params.vars_2_regress
    )

    ch_pca_multiple = PCA_MULT (
        ch_merged_so.rds
    )
    ch_pca_single = PCA_SING (
        ch_scaled_so.rds
    )

    ch_pca_single.rds.join(ch_pca_single.log).set { ch_pca_single_updated }
    ch_pca_multiple.rds.join(ch_pca_multiple.log).set { ch_pca_multiple_updated }

    ch_pca_single_updated.map { meta, rds, logs ->
                            if (true) {
                                new_meta = new LinkedHashMap()
                                new_meta['group'] = meta
                                new_meta['integrated'] = false
                                return [ new_meta, rds, logs ]
                            } else {
                                return [ meta, rds, logs ]
                            }
                        }
                        .set { ch_pca_single_updated }

    if (params.integration_method){
        ch_integrated = INTEGRATION (
            ch_pca_multiple_updated.map { meta, rds, logs -> [ meta, rds ] },
            params.integration_method
        )
        ch_integrated.rds.join(ch_pca_multiple.log).set { ch_integrated_all }
        ch_integrated_all.map { meta, rds, logs ->
                                    if (true) {
                                        new_meta = new LinkedHashMap()
                                        new_meta['group'] = meta
                                        new_meta['integrated'] = true
                                        return [ new_meta, rds, logs ]
                                    } else {
                                        return [ meta, rds, logs ]
                                    }
                                }.set { ch_dimensions_def }
    } else {
        ch_pca_multiple.rds.map { meta, rds, logs ->
                                    if (true) {
                                        new_meta = new LinkedHashMap()
                                        new_meta['group'] = meta
                                        new_meta['integrated'] = false
                                        return [ new_meta, rds, logs ]
                                    } else {
                                        return [ meta, rds, logs ]
                                    }
                                }
                                .set { ch_dimensions_def }
    }

    ch_dim_def_all = Channel.empty()
    ch_dim_def_all.mix(ch_pca_single_updated).set { ch_dim_def_all }
    ch_dim_def_all.mix(ch_dimensions_def).set { ch_dim_def_all }
    //ch_dim_def_all.view()

    ch_nn_clusters = FIND_NN_CLUSTER (
        ch_dim_def_all.map {meta, rds, log -> [meta, rds]},
        ch_dim_def_all.map {meta, rds, log -> [meta, log]},
        params.resolutions,
        params.integration_method
    )


    ch_nn_clusters.rds.map { meta, rds -> [ meta.group, meta.integrated, rds ] }
                        .join( ch_pca_multiple.log, by: [0, 0], remainder: true)
                        .join( ch_pca_single.log, by: [0, 0], remainder: true)
                        .map { list -> list.findAll { it != null }}
                        .set { ch_nn_clusters_w_log }

    ch_nn_clusters_w_log.view()
    DISPLAY_REDUCTION (
        ch_nn_clusters_w_log.map { [ it[0], it[2] ] },
        ch_nn_clusters_w_log.map { [ it[0], it[3] ] },
        ch_nn_clusters_w_log.map { it[1] },
        params.resolutions,
        params.makeLoupe,
        params.integration_method
    )
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
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
