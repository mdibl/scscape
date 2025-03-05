process FEATURE_NAMING {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/rgrindle/convert_gn_nms:latest' }"

    input:
    tuple val(meta), path(sample_files), path(gene_list)

    output:
    tuple val(meta), path (sample_files), path ("AuxGeneList.csv"), emit: data
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    while true; do
        sleep 1
    done
    convert.py ${sample_files}/features.tsv.gz $gene_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BASH: \$(echo \$(bash --version) )
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BASH: \$(echo \$(bash --version) )
    END_VERSIONS
    """

}

