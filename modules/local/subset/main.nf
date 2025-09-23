process SUBSET {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'docker.io/mdiblbiocore/seurat:latest' }"

    input:
    tuple val(meta), path (data_directory), path (genes_2_rm)
    val min_cells
    val min_features
    val gene_identifier

    output:
    tuple val(meta), path ("*SO.rds"),        emit: rds
    tuple val(meta), path("*Validation.log"), emit: log
    path ("*FinalVersions.log"),                     emit: r_versions
    path ('versions.yml'), emit: versions
    path("*Execution.log"), emit: exec

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    if (genes_2_rm.getClass().name == "nextflow.util.BlankSeparatedList") {
        genes_2_rm = 'NULL'
    }
    """
    MakeSeurat.R \\
        $data_directory \\
        $gene_identifier \\
        $genes_2_rm \\
        ${meta.id} \\
        $min_cells \\
        $min_features \\
        "$meta" \\
        ${args} 2>&1 | tee 00_${meta.id}_Execution.log

    perl -i -pe 's/"//g;s/\\[\\d\\d?\\d?\\] //g' 00_${meta.id}_InitialValidation.log

    grep -i -E "R version " 00_${meta.id}_InitialVersions.log | perl -pe 's/ version /: \\"/g;s/ \\(.*/\\"/g' >> 00_${meta.id}_FinalVersions.log
    perl -ne 'print if /other attached packages:/ .. /^\$/' 00_${meta.id}_InitialVersions.log | grep -v "other" | perl -pe 's/\\\\[.*]\\s+//g;s/\\s+/\\n/g' | grep -v "^\$" | perl -pe 's/_/: \\"/g;s/\$/\\"/' >> 00_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

    stub:
    """
    touch 00_${meta.id}_InitialVersions.log
    touch 00_${meta.id}_FinalVersions.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version| head -n 1| grep -Eo "[0-9]+[^ ]*"| head -n 1) )
    END_VERSIONS
    """

}
