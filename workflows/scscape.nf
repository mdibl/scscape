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
include { MAKE_SEURAT          } from '../modules/local/makeseurat'
include { NORMALIZE_QC         } from '../modules/local/normalize_qc'
include { FIND_DOUBLETS        } from '../modules/local/doubletfinder'
include { MERGE_SO             } from '../modules/local/merge'
include { MERGE_SO as SCALE_SO } from '../modules/local/merge'
include { RUN_PCA as PCA_MULT  } from '../modules/local/runpca'
include { RUN_PCA as PCA_SING  } from '../modules/local/runpca'
include { INTEGRATION          } from '../modules/local/integration'
include { FIND_NN_CLUSTER      } from '../modules/local/find_NN_clusters'
include { DISPLAY_REDUCTION    } from '../modules/local/plotting'
include { GZIP                 } from '../modules/local/gzip'
include { FEATURE_NAMING       } from '../modules/local/feature_naming'

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

    ch_gzip = GZIP(ch_meta_samples.map { [ it[1], file(it[2][0].trim()) ] })
    ch_gzip.zip
            .map { meta, data ->
                meta.eachWithIndex { entry, i ->
                    if (entry.value instanceof ArrayList) {
                        meta[entry.key] = entry.value.sort()
                    }
                }
                return [ meta.id, meta, data ]
            }
            .join( ch_meta_samples, by: [0,0])
            .map { id, meta, gz, meta1, data -> [ meta, gz, file(data[1].trim()) ] }
            .set { ch_meta_samples }


    if (params.gene_identifier.toUpperCase() == "COMBINE"){
        ch_updated_features = FEATURE_NAMING(
            ch_meta_samples.map { [ it[0], it[1],  it[2] ] },
        )
        ch_init_rds = MAKE_SEURAT (
        ch_updated_features.data.map{ [ it[0], it[1], it[2] ] },
        params.min_cells,
        params.min_features,
        params.gene_identifier
        )
        ch_init_rds.rds.map { meta, rds -> [ meta.id, meta, rds ]}
        .join(ch_meta_samples)
        .map{ id, meta, rds, genes ->  [ meta, rds, genes ] }
        .set { ch_init_rds_meta }
        ch_validation_log.mix(ch_init_rds.log).set{ ch_validation_log }
    } else {

        ch_meta_samples.map { meta, data, genes ->
                meta.eachWithIndex { entry, i ->
                    if (entry.value instanceof ArrayList) {
                        meta[entry.key] = entry.value.sort()
                    }
                }
                return [ meta, data, genes ]
            }.set{ ch_sorted_meta_mk_seurat }

        ch_init_rds = MAKE_SEURAT (
        ch_sorted_meta_mk_seurat.map { [ it[0] , it[1], it[2] ] },
        params.min_cells,
        params.min_features,
        params.gene_identifier
        )

        ch_meta_samples.map { meta, data, genes ->
                        meta.eachWithIndex { entry, i ->
                            if (entry.value instanceof ArrayList) {
                                meta[entry.key] = entry.value.sort()
                            }
                        }
                        return [ meta.id, meta, data, genes ]
                    }.set { ch_meta_samples }

        ch_init_rds.rds.map { meta, rds -> [ meta.id, meta, rds ]}
        .join(ch_meta_samples)
        .map{ id, meta, rds, meta1, data, genes ->  [ meta, rds, genes ] }
        .set { ch_init_rds_meta }
        ch_validation_log.mix(ch_init_rds.log).set{ ch_validation_log }

    }

    ch_init_rds_meta.map { [ it[0] , it[1], it[2] ] }
                    .map { meta, rds, genes ->
                        meta.eachWithIndex { entry, i ->
                            if (entry.value instanceof ArrayList) {
                                meta[entry.key] = entry.value.sort()
                            }
                        }
                        return [ meta, rds, genes ]
                    }.set{ ch_sorted_meta_norm_qc }



    ch_normalized_qc = NORMALIZE_QC (
        ch_sorted_meta_norm_qc.map { [it[0], it[1], it[2]] },
        params.nfeature_lower,
        params.nfeature_upper,
        params.ncount_lower,
        params.ncount_upper,
        params.max_mito_pct,
        params.vars_2_regress
    )

    ch_normalized_qc.rds
    .map { meta, rds -> [ meta.id, meta, rds ]}
    .join(ch_meta_samples)
    .map { id, meta, rds, meta1, data, genes -> [ meta, rds, data ] }
    .set { ch_normalized_qc_meta }
    ch_validation_log.mix(ch_normalized_qc.log).set{ ch_validation_log }

    ch_normalized_qc_meta.map { [ it[0] , it[1], it[2] ] }
                        .map { meta, rds, data ->
                            meta.eachWithIndex { entry, i ->
                                if (entry.value instanceof ArrayList) {
                                    meta[entry.key] = entry.value.sort()
                                }
                            }
                            return [ meta, rds, data ]
                        }.set{ ch_sorted_meta_doublets }

    ch_doublet_filtered_rds = FIND_DOUBLETS (
        ch_sorted_meta_doublets.map { [it[0] , it[1], it[2]] },
        params.vars_2_regress
    )
    ch_validation_log.mix(ch_doublet_filtered_rds.log).set{ ch_validation_log }

    ch_doublet_filtered_rds.rds.map { meta, rds -> [ rds, meta.groups.toList() ]}
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
    ch_validation_log.mix(ch_merged_so.log).set{ ch_validation_log }

    ch_scaled_so = SCALE_SO (
        ch_merged_groups.single,
        params.vars_2_regress
    )
    ch_validation_log.mix(ch_scaled_so.log).set{ ch_validation_log }

    ch_pca_multiple = PCA_MULT (
        ch_merged_so.rds
    )
    ch_validation_log.mix(ch_pca_multiple.log).set{ ch_validation_log }

    ch_pca_single = PCA_SING (
        ch_scaled_so.rds
    )
    ch_validation_log.mix(ch_pca_single.log).set{ ch_validation_log }

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

        ch_validation_log.mix(ch_integrated.log).set{ ch_validation_log }
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

    ch_nn_clusters = FIND_NN_CLUSTER (
        ch_dim_def_all.map {meta, rds, log -> [meta, rds, log]},
        params.resolutions,
        params.integration_method
    )
    ch_validation_log.mix(ch_nn_clusters.log).set{ ch_validation_log }


    ch_nn_clusters.rds.map { meta, rds -> [ meta.group, meta.integrated, rds ] }
                        .join( ch_pca_multiple.log, by: [0, 0], remainder: true)
                        .join( ch_pca_single.log, by: [0, 0], remainder: true)
                        .map { list -> list.findAll { it != null }}
                        .set { ch_nn_clusters_w_log }

    DISPLAY_REDUCTION (
        ch_nn_clusters_w_log.map { [ it[0], it[2], it[3] ] },
        ch_nn_clusters_w_log.map { it[1] },
        params.resolutions,
        params.makeLoupe,
        params.integration_method
    )
    ch_validation_log.mix(DISPLAY_REDUCTION.out.log).set{ ch_validation_log }


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
