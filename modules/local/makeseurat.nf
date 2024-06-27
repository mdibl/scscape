process MAKE_SEURAT {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (data_directory)
    tuple val(meta), path (genes_2_rm)
    val min_cells
    val min_features
    val gene_identifier

    output:
    tuple val(meta), path ("*_SO.rds"), emit: rds
    path("*.validation.log"),           emit: log
    path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    if (genes_2_rm.getClass().name == "nextflow.util.BlankSeparatedList") {
        genes_2_rm = 'null'
    }
    """
    MakeSeurat.R \\
        $data_directory \\
        $gene_identifier \\
        $genes_2_rm \\
        ${meta.id} \\
        ${meta.group} \\
        $min_cells \\
        $min_features \\
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
        Seurat: \$(Rscript -e "packageVersion('Seurat')")
    END_VERSIONS
    """

}
