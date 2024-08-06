process INTEGRATION {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    val integration_method

    output:
    tuple val(meta), path ("*_IntegrateSO.rds"), emit: rds
    tuple val(meta), path("*Validation.log"),           emit: log
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def scale_method = task.ext.args2 ?: 'SCT'
    """
    Integration.R \\
        $rds \\
        $integration_method \\
        ${meta} \\
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
