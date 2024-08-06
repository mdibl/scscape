process NORMALIZE_QC {

    tag "${meta.id}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    tuple val(meta), path (mt_cc_genes)
    val nfeature_lower
    val nfeature_upper
    val ncount_lower
    val ncount_upper
    val max_mito_pct
    val vars_2_regress


    output:
    tuple val(meta), path ("*_NormQCSO.rds"),        emit: rds
    tuple val(meta), path("*Validation.log"), emit: log
    path("*.pdf")
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    if (mt_cc_genes.getClass().name == "nextflow.util.BlankSeparatedList") {
        mt_cc_genes = 'NULL'
    }
    if (vars_2_regress.toUpperCase().contains("S.SCORE") && vars_2_regress.toUpperCase().contains("G2M.SCORE")) {
        run_cc_score = "true"
    } else if (vars_2_regress.toUpperCase().contains("S.SCORE") && !(vars_2_regress.toUpperCase().contains("G2M.SCORE"))){
        throw new IllegalArgumentException("Only 1 cell cycle variable in regressed features, please include both G2M.Score, S.Score, or neither.")
    } else if (!(vars_2_regress.toUpperCase().contains("S.SCORE")) && vars_2_regress.toUpperCase().contains("G2M.SCORE")){
        throw new IllegalArgumentException("Only 1 cell cycle variable in regressed features, please include both G2M.Score, S.Score, or neither.")
    } else {
        run_cc_score = "false"
    }
    """
    Normalize_QC.R \\
        $mt_cc_genes \\
        $nfeature_lower \\
        $nfeature_upper \\
        $ncount_lower \\
        $ncount_upper \\
        $max_mito_pct \\
        $rds \\
        ${meta.id} \\
        $run_cc_score \\
        ${args}

    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(echo \$(Seurat --version) | sed "s/Seurat, version //g" )
    END_VERSIONS
    """

}
