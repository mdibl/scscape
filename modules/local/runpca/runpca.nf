process RUN_PCA {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)

    output:
    tuple val(meta), path ("*_PCASO.rds"), emit: rds
    tuple val(meta), path("*Validation.log"),           emit: log
    path("*.pdf")
    path ("*FinalVersions.log"),                     emit: r_versions
    path("versions.yml"), emit: versions
    path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    def pcMax = task.ext.args  ?: 'null'
    def args  = task.ext.args2  ?: ''
    """
    RunPCA.R \\
        $rds \\
        ${pcMax} \\
        ${meta} \\
        ${args} 2>&1 | tee > 05_${meta.id}_Execution.log

    grep -i -E "R version " 05_${meta}_InitialVersions.log | perl -pe 's/ version /: "/g;s/ \(.*/"/g' >> 05_${meta}_FinalVersions.log
    perl -ne 'print if /other attached packages:/ .. /^\$/' 05_${meta}_InitialVersions.log | grep -v "other" | perl -pe 's/\\[.*]\s+//g;s/\s+/\n/g' | grep -v "^\$" | perl -pe 's/_/: "/g;s/\$/"/' >> 05_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 05_${meta}_InitialVersions.log
    touch 05_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}