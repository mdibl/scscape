process NORMALIZE_QC {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    tuple val(meta), path (mito_genes)
    val nfeature_lower
    val nfeature_upper
    val ncount_lower
    val ncount_upper
    val max_mito_pct
    

    output:
    tuple val(meta), path ("*_QC.rds"), emit: rds
    path("*.validation.log"),           emit: log
    //path ("versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
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
