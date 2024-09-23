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
    tuple val(meta), path ("*_FinalSO.rds"), emit: rds
    //path("*Validation.log"),           emit: log
    path("*.cloupe")
    path("*.pdf")
    //path ("*FinalVersions.log"),                     emit: r_versions
    //path('versions.yml'), emit: versions
    //path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    if (integrated) {
        integration_method = integration_tool
    } else {
        integration_method = "NULL"
    }
    def args = task.ext.args  ?: ''
    def eula_agreement = task.ext.args2 ?: 'NULL'
    """
    pcMax=\$(grep -i -E "^PCs used" $validation_log | perl -pe "s/^PCs.*- //g")

    Plotting.R \
        $rds_file \
        $resolutions \
        \$pcMax \
        $integration_method \
        $meta \
        $makeLoupe \
        $eula_agreement \
        ${args} 2>&1 | tee > 04_${meta}_Execution.log

    ##grep -i -E "R version " 08_${meta}_InitialVersions.log | perl -pe 's/ version /: "/g;s/ \(.*/"/g' >> 08_${meta}_FinalVersions.log
    ##perl -ne 'print if /other attached packages:/ .. /^\$/' 08_${meta}_InitialVersions.log | grep -v "other" | perl -pe 's/\\[.*]\s+//g;s/\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: "/g;s/\$/"/' >> 08_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 08_${meta}_InitialVersions.log
    touch 08_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
