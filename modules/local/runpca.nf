process RUN_PCA {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:v2.1.0' }"

    input:
    tuple val(meta), path (rds)

    output:
    tuple val(meta), path ("*_PCA.rds"), emit: rds
    tuple val(meta), path("*.validation.log"),           emit: log
    path("*.pdf")
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pcMax = task.ext.args  ?: 'null'
    def args  = task.ext.args2  ?: ''
    """
    RunPCA.R \\
        $rds \\
        ${pcMax} \\
        ${meta} \\
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
