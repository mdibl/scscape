process MERGE_SO {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds_files)
    val vars_2_regress

    output:
    tuple val(meta), path ("*_MergedSO.rds"), emit: rds
    path("*Validation.log"),           emit: log
    path("*.pdf")
    path ("*FinalVersions.log"),                     emit: r_versions
    path("versions.yml"), emit: versions
    path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    def n_features = task.ext.args2 ?: 'NULL'
    def scale_method = task.ext.args3 ?: 'NULL'
    def args = task.ext.args  ?: ''
    """
    Merge.R \\
        ${meta} \\
        $vars_2_regress \\
        $n_features \\
        $scale_method \\
        ${args} 2>&1 | tee 03_${meta}_Execution.log

    perl -i -pe 's/"//g;s/\\[\\d\\d?\\d?\\] //g' 03_${meta}_MergeValidation.log

    grep -i -E "R version " 03_${meta}_MergeVersions.log | perl -pe 's/ version /: \\"/g;s/ \\(.*/\\"/g' >> 03_${meta}_FinalVersions.log
    perl -ne 'print if /other attached packages:/ .. /^\$/' 03_${meta}_MergeVersions.log | grep -v "other" | perl -pe 's/\\\\[.*]\\s+//g;s/\\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: \\"/g;s/\$/\\"/' >> 03_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 04_${meta}_MergeVersions.log
    touch 04_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
