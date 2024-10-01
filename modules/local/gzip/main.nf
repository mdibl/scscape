process GZIP {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(sample_files)

    output:
    tuple val(meta), path (sample_files), emit: zip
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for file in ${sample_files}/*; do
        if [[ \$file != *".gz"* ]]; then
            gzip \$file
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GZIP: \$(echo \$(gzip --version| head -n 1| sed 's/gzip //') )
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GZIP: \$(echo \$(gzip --version) )
    END_VERSIONS
    """

}

