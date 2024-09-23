process INTEGRATION {

    tag "${meta}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    val integration_method

    output:
    tuple val(meta), path ("*_IntegrateSO.rds"), emit: rds
    //tuple val(meta), path("*Validation.log"),           emit: log
    //path ("*FinalVersions.log"),                     emit: r_versions
    //path("versions.yml"), emit: versions
    //path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def scale_method = task.ext.args2 ?: 'SCT'
    """
    Integration.R \\
        $rds \\
        $integration_method \\
        ${meta} \\
        $scale_method \\
        ${args} 2>&1 | tee > 04_${meta}_Execution.log

    ##grep -i -E "R version " 06_${meta}_InitialVersions.log | perl -pe 's/ version /: "/g;s/ \(.*/"/g' >> 06_${meta}_FinalVersions.log
    ##perl -ne 'print if /other attached packages:/ .. /^\$/' 06_${meta}_InitialVersions.log | grep -v "other" | perl -pe 's/\\[.*]\s+//g;s/\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: "/g;s/\$/"/' >> 06_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 06_${meta}_InitialVersions.log
    touch 06_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
