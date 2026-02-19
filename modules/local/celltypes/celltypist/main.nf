process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/celltype_assignment/celltypist", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/celltypist_scanpy:89a98f51262cfff4'
        : 'community.wave.seqera.io/library/celltypist_scanpy:44b604b24dd4cf33'}"

    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    val models
    val save_probabilities

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*_celltypist.pkl"), emit: obs
    path "*_celltypist_probabilities.parquet.gz", emit: probabilities
    path "*_celltypist_probabilities_metadata.csv", emit: probabilities_metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('celltypist.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch versions.yml
    """
}
