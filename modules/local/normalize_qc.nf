process NORMALIZE_QC {

    tag "${meta.id}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    tuple val(meta), path (mito_genes)
    val nfeature_lower
    val nfeature_upper
    val ncount_lower
    val ncount_upper
    val max_mito_pct


    output:
    tuple val(meta), path ("*_QC.rds"),        emit: rds
    tuple val(meta), path("*.validation.log"), emit: log
    path("*.pdf")
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    if (mito_genes.getClass().name == "nextflow.util.BlankSeparatedList") {
        mito_genes = 'NULL'
    }
    """
    Normalize_QC.R \\
        $mito_genes \\
        $nfeature_lower \\
        $nfeature_upper \\
        $ncount_lower \\
        $ncount_upper \\
        $max_mito_pct \\
        $rds \\
        ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(echo \$( version) | sed "s/, version //g" )
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(echo \$(Seurat --version) | sed "s/Seurat, version //g" )
    END_VERSIONS
    """

}
