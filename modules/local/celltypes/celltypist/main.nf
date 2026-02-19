process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'
    // Conditional publishing for probability files
    publishDir path: "${params.outdir}/celltype_assignment/celltypist", 
        mode: 'copy',
        enabled: params.celltypist_save_probabilities == true,
        saveAs: { filename -> 
            filename.contains('_probabilities.') ? filename : null 
        }
    
    // Always publish main results
    publishDir path: "${params.outdir}/celltype_assignment/celltypist", 
        mode: 'copy',
        saveAs: { filename -> 
            (filename.endsWith('.h5ad') || filename.endsWith('_celltypist.pkl') || filename == 'versions.yml') ? filename : null
        }

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
    path("*_probabilities.parquet.gz"), emit: probabilities, optional: true
    path("*_probabilities_metadata.csv"), emit: probabilities_metadata, optional: true
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
