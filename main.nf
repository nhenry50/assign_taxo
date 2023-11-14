#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    taxonomic assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : ??????
----------------------------------------------------------------------------------------
*/

process GET_REF {

    label 'process_low'

    conda "r::r-tidyverse=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/r-tidyverse:1.2.1' }"

    output:
    path("refdb.rds")

    """
    #!/usr/bin/env Rscript

    options(timeout = 1200)

    if ("${params.refdb}" == "silva"){

        download.file(
            "http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData",
            "SILVA_SSU_r138_2019.RData",
            quiet = TRUE
        )

        load("SILVA_SSU_r138_2019.RData", ex <- new.env())

        tmp <- get(names(ex)[1],envir=ex)

        saveRDS(tmp,file="refdb.rds")

    } else if ("${params.refdb}" == "pr2") {

        download.file(
            "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU.decipher.trained.rds",
            "refdb.rds",
            quiet = TRUE
        )

    }
    
    """
}

process ASSIGN_IDTAXA {

    label 'process_low'

    conda "bioconda::bioconductor-decipher=2.28.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-decipher:2.28.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-decipher:2.28.0--r43ha9d7317_0' }"

    input:
    path(fasta_seq)
    path(refdb)

    output:
    path("taxonomy_${params.refdb}.tsv")

    """
    #!/usr/bin/env Rscript

    ##########################################
    # import sequences
    ##########################################

    sequences <- Biostrings::readDNAStringSet("${fasta_seq}", format = "fasta")
    names(sequences) <- gsub("^centroid=|;.+\$", "", names(sequences))

    ##########################################
    # assign with idtaxa
    ##########################################

    trainingSet <- readRDS("${refdb}")

    idtaxa_res <- DECIPHER::IdTaxa(
        sequences,
        trainingSet,
        strand = "top",
        threshold = ${params.idtaxa_thresh},
        minDescend = 0.9,
        processors = 1
    )

    taxonomy <- vapply(
        idtaxa_res,
        function(x) paste(x[["taxon"]], collapse = ";"),
        character(1)
    )

    confidence <- vapply(
        idtaxa_res,
        function(x) paste(round(x[["confidence"]], digits = 1), collapse = ";"),
        character(1)
    )


    ##########################################
    # assemble into one table
    ##########################################

    res <- data.frame(
        sequence = names(taxonomy),
        taxonomy = sub("^Root;", "", taxonomy),
        confidence = sub("^[^;]+;", "",confidence)
        )

    write.table(
        res,
        file = "taxonomy_${params.refdb}.tsv",
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

    """
}

workflow {

    GET_REF().set{ refbd }

    Channel.fromPath(params.inputfasta)
        .splitFasta( by: params.fastachunks, file: true )
        .set{ splitted_fasta }

    ASSIGN_IDTAXA(splitted_fasta, refbd)
        .collectFile(keepHeader: true, skip: 1, storeDir: params.outdir, name: "taxores")

}