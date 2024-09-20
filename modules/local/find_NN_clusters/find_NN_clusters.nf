process FIND_NN_CLUSTER {

    tag "${meta.group}"
    label 'process_small'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (rds)
    tuple val(meta), path (validation_log)
    val resolutions
    val integration_tool

    output:
    tuple val(meta), path ("*_ClusterSO.rds"), emit: rds
    //tuple val(meta), path("*Validation.log"),  emit: log
    path("markers")
    path("*.pdf")
    //path ("*FinalVersions.log"),               emit: r_versions
    //path('versions.yml'), emit: versions
    //path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta.integrated) {
        integration_method = integration_tool
    } else {
        integration_method = "NULL"
    }
    def args = task.ext.args ?: ''
    def scale_method = task.ext.args2 ?: 'SCT'
    """
    pcMax=\$(grep -i -E "^PCs used" $validation_log | perl -pe "s/^PCs.*- //g")

    echo \$pcMax
    FindNeighborsClustersMarkers.R \\
        $rds \\
        $resolutions \\
        \$pcMax \\
        ${integration_method} \\
        ${meta.group} \\
        $scale_method \\
        ${args} 2>&1 | tee > 04_${meta.id}_Execution.log

    ##grep -i -E "R version " 07_${meta}_InitialVersions.log | perl -pe 's/ version /: "/g;s/ \(.*/"/g' >> 07_${meta}_FinalVersions.log
    ##perl -ne 'print if /other attached packages:/ .. /^\$/' 07_${meta}_InitialVersions.log | grep -v "other" | perl -pe 's/\\[.*]\s+//g;s/\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: "/g;s/\$/"/' >> 07_${meta}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 07_${meta}_InitialVersions.log
    touch 07_${meta}_FinalVersions.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
