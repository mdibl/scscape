process FIND_NN_CLUSTER {

    tag "${meta.group}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:v2.1.0' }"

    input:
    tuple val(meta), path (rds)
    tuple val(meta), path (validation_log)
    val resolutions
    val integration_tool

    output:
    tuple val(meta), path ("*_Clustered.rds"), emit: rds
    tuple val(meta), path("*validation.log"),                   emit: log
    path("markers")
    path("*.pdf")
    //path ("versions.yml"),                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta.integrated) {
        integration_method = integration_tool
    } else {
        integration_method = "NULL"
    }
    def args = task.ext.args ?: ''
    """
    pcMax=\$(paste -s  <(grep PC $validation_log| grep -E -o "[0-9]") | sed 's|\\t||')

    FindNeighborsClustersMarkers.R \\
        $rds\\
        $resolutions \\
        \$pcMax \\
        ${integration_method} \\
        ${meta.group} \\
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
