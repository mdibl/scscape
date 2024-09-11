process FIND_DOUBLETS {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:

    tuple val(meta), path (rds)
    tuple val(meta), path (data_directory)
    val vars_2_regress



    output:
    tuple val(meta), path ("*_DoubletsRmSO.rds"),     emit: rds
    tuple val(meta), path("*Validation.log"),           emit: log
    path("*.pdf")
    path ("*FinalVersions.log"),                     emit: r_versions
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def n_features = task.ext.args2 ?: 'NULL'
    def scale_method = task.ext.args3 ?: 'NULL'
    def args = task.ext.args  ?: ''
    """
    DoubletFinder.R \\
        $rds \\
        $vars_2_regress \\
        $data_directory \\
        ${meta.id} \\
        $n_features \\
        $scale_method \\
        ${args}

    grep -i -E "R version " 02_${meta.id}_InitialVersions.log | perl -pe 's/ version /: "/g;s/ \(.*/"/g' >> 02_${meta.id}_FinalVersions.log
    perl -ne 'print if /other attached packages:/ .. /^\$/' 02_${meta.id}_InitialVersions.log | grep -v "other" | perl -pe 's/\[.*]\s+//g;s/\s+/\n/g' | grep -v "^\$" | perl -pe 's/_/: "/g;s/\$/"/' >> 02_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 02_${meta.id}_InitialVersions.log
    touch 02_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
