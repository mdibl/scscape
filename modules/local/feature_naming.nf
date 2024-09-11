process FEATURE_NAMING {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(sample_files)
    tuple val(meta), path(gene_list)

    output:
    tuple val(meta), path (sample_files), path ("AuxGeneList.csv")
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash

    gzcat ${sample_files}/features.tsv.gz | perl -ane 'if (\$F[0] ne \$F[1]) { print "\$F[0]\t\$F[1]::\$F[0]\tExpression\n"; } else { print "\$F[0]\t\$F[1]\tExpression\n"; }' | gzip > features_new.tsv.gz

    rm ${sample_files}/features.tsv.gz
    mv features_new.tsv.gz ${sample_files}/features.tsv.gz

    echo "MTgenes" > MT.csv
    cut -f1 -d ","  $gene_list | grep -v "^\$" | tail -n +2 | perl -pe "s/^/\\t/;s/\$/::/" > origMT.csv
    gzcat ${sample_files}/features.tsv.gz | grep -f origMT.csv | cut -f2 >> MT.csv
    rm origMT.csv

    echo "G2Mgenes" > G2M.csv
    cut -f2 -d ","  $gene_list | grep -v "^\$" | tail -n +2 | perl -pe "s/^/\\t/;s/\$/::/"> origG2M.csv
    gzcat ${sample_files}/features.tsv.gz | grep -f origG2M.csv | cut -f2 >> G2M.csv
    rm origG2M.csv

    echo "Sgenes" > S.csv
    cut -f3 -d ","  $gene_list | grep -v "^\$" | tail -n +2 | perl -pe "s/^/\\t/;s/\$/::/"> origS.csv
    gzcat ${sample_files}/features.tsv.gz | grep -f origS.csv | cut -f2 >> S.csv
    rm origS.csv

    echo "RMgenes" > RM.csv
    cut -f4 -d ","  $gene_list | grep -v "^\$" | tail -n +2 | perl -pe "s/^/\\t/;s/\$/::/"> origRM.csv
    gzcat ${sample_files}/features.tsv.gz | grep -f origRM.csv | cut -f2 >> RM.csv
    rm UpdatedFiles/origRM.csv

    paste -d ',' MT.csv G2M.csv S.csv RM.csv > AuxGeneList.csv

    rm -r UpdatedFiles


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

