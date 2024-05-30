process INTEGRATION {

    tag "${meta.group}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    val integration_method
    
    output:
    tuple val(meta), path ("*_Integrated.rds"), emit: rds
    path("*validation.log"),           emit: log
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    Integration.R \\
        $rds \\
        $integration_method \\
        ${meta.group} \\
        ${args}      

    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
        //Seurat: \$(echo \$( version) | sed "s/, version //g" )
    //END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(echo \$(Seurat --version) | sed "s/Seurat, version //g" )
    END_VERSIONS
    """

} 
