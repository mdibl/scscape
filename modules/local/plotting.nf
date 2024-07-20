process DISPLAY_REDUCTION {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds_file)
    tuple val(meta), path(validation_log)
    val integrated
    val resolutions
    val makeLoupe
    val integration_tool



    output:
    tuple val(meta), path ("*Final.rds"), emit: rds
    path("*validation.log"),           emit: log
    path("*.cloupe")
    path("*.pdf")
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (integrated) {
        integration_method = integration_tool
    } else {
        integration_method = "NULL"
    }
    def args = task.ext.args1  ?: 'NULL'
    """
    pcMax=\$(paste -s  <(grep PC $validation_log| grep -E -o "[0-9]") | sed 's|\\t||')

    Plotting.R \
        $rds_file \
        $resolutions \
        \$pcMax \
        $integration_method \
        $meta \
        $makeLoupe \
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
