include { MERGE_SO             } from '../../../modules/local/merge'
include { MERGE_SO as SCALE_SO } from '../../../modules/local/merge'
include { RUN_PCA as PCA_MULT  } from '../../../modules/local/runpca'
include { RUN_PCA as PCA_SING  } from '../../../modules/local/runpca'
include { INTEGRATION          } from '../../../modules/local/integration'
include { FIND_NN_CLUSTER      } from '../../../modules/local/find_NN_clusters'
include { DISPLAY_REDUCTION    } from '../../../modules/local/plotting'

workflow POST_QC {
    take:
    ch_merged_groups
    ch_validation_log


    main:

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
                                def new_meta = new LinkedHashMap()
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
                                        def new_meta = new LinkedHashMap()
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
                                        def new_meta = new LinkedHashMap()
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

    emit:
    validation_log = ch_validation_log
    final_rds = DISPLAY_REDUCTION.out.rds

}
