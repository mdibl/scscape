include { MAKE_SEURAT          } from '../../../modules/local/makeseurat'
include { NORMALIZE_QC         } from '../../../modules/local/normalize_qc'
include { FIND_DOUBLETS        } from '../../../modules/local/doubletfinder'
include { GZIP                 } from '../../../modules/local/gzip'
include { FEATURE_NAMING       } from '../../../modules/local/feature_naming'


workflow QC {

    take:
    ch_meta_samples
    ch_validation_log

    main:

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



    emit:
    validation_log = ch_validation_log
    doublet_filtered_rds = ch_merged_groups

}
