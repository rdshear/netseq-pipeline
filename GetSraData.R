library(tidyverse)
library(xml2)
library(rentrez)
library(glue)

# filetype in "fastq", "run"
create_sample_grid <- function(bioproject,
            semantic_name = "fastq", org = "GCP", retmax = NULL) {

        bpid <- entrez_search(db = "bioproject", term = bioproject)
    if (bpid$count == 0) {
        stop(glue::glue("Bioproject <{bioproject}> not found"))
    }

    sra_ids <- entrez_link(dbfrom = "bioproject", id = bpid$ids, db = "sra")$links$bioproject_sra
    x <- read_xml(entrez_fetch(db = "sra", id = sra_ids, rettype = "xml", retmax = retmax))
    
    runs <- xml_find_all(x, '/EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN')
    run_id <- xml_attr(runs, "accession")
    expref <- xml_find_all(runs, "EXPERIMENT_REF")
    experiment_id <- xml_attr(expref, "accession")
    biosample_id <- xml_attr(expref, "refname")
    # Only reporting first member in pool
    member <- xml_find_first(runs, "Pool/Member")
    sample_title <- xml_attr(member, "sample_title")
    sra_sample_id <- xml_attr(member, "accession")
    sra_File <- xml_find_first(runs, glue("./SRAFiles/SRAFile[@semantic_name=\"{semantic_name}\"]"))
    sra_File_alt <- xml_find_first(sra_File, glue("./Alternatives[@org=\"{org}\"]"))
    filename <- xml_attr(sra_File, "filename")
    sra_url <- xml_attr(sra_File_alt, "url")
    tibble(bioproject, experiment_id, biosample_id, sra_sample_id, run_id, sample_title,  filename, sra_url)
}

# PRJNA669852 -- Churchman Lab - dozens of regulatory factors
# PRJNA668299 -- Mellor Lab - Spt4
result <- create_sample_grid("PRJNA669852")

# TODO Clean up sample_id's if we are going to allow multiple assays
# Infer strain and assay type from sample_title
result %>% 
    separate(sample_title, into = c("sample_id", "assay"), sep = "_", remove = FALSE) %>%
    separate(sample_id, into = "strain", sep = "-", extra = "drop", remove = FALSE) %>%
    filter(assay == "NETseq") %>%
    relocate(sample_id) %>%
    arrange(sample_title) %>%
    rename(inputFastQ = sra_url, "entity:sample_id" = sample_id) -> out


write.table(out, file = "~/Downloads/sample.tsv", sep = "\t", row.names = FALSE)




