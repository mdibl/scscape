process MERGE_SO {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:v2.1.0' }"

    input:
    tuple val(meta), path (rds_files)
    val vars_2_regress

    output:
    tuple val(meta), path ("*Merged_SO.rds"), emit: rds
    path("*validation.log"),           emit: log
    path("*.pdf")
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def n_features = task.ext.args2 ?: 'ALL'
    def scale_method = task.ext.args3 ?: 'SCT'
    def args = task.ext.args  ?: ''
    """
    Merge.R \\
        ${meta} \\
        $vars_2_regress \\
        $n_features \\
        $scale_method \\
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
