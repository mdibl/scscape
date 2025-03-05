process NORMALIZE_QC {

    tag "${meta.id}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds), path (mt_cc_genes)
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
    path ("*FinalVersions.log"),                     emit: r_versions
    path("versions.yml"), emit: versions
    path("*Execution.log"), emit: exec

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
        ${args} 2>&1 | tee 01_${meta.id}_Execution.log

    perl -i -pe 's/"//g;s/\\[\\d\\d?\\d?\\] //g' 01_${meta.id}_NormQCValidation.log

    grep -i -E "R version " 01_${meta.id}_NormQCVersions.log | perl -pe 's/ version /: \\"/g;s/ \\(.*/\\"/g' >> 01_${meta.id}_FinalVersions.log
    perl -ne 'print if /other attached packages:/ .. /^\$/' 01_${meta.id}_NormQCVersions.log | grep -v "other" | perl -pe 's/\\\\[.*]\\s+//g;s/\\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: \\"/g;s/\$/\\"/' >> 01_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 01_${meta.id}_NormQCVersions.log
    touch 01_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
