setwd('ITS_Fungi_Metabarcoding/RESULTS/')

# Required packages
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

source("../../UK_SURVEY/ANALYSIS_SCRIPTS/get_cutoffs.R")
source("../../UK_SURVEY/ANALYSIS_SCRIPTS/taxa_match.R")

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
taxa_file_path <- 'dnabarcoder/dnabarcoder/all_nochimeras_formatted.Fungi.ITS2.unite2024ITS2_BLAST.classification.dynamic.format.nosize'
cutoff_file_path <- 'dnabarcoder/taxonomic_cutoffs.csv'
asv_sequences_path <- 'ITSx/all_nochimeras_formatted.Fungi.ITS2.nosize.fasta'
asv_table_path <- 'ALLimportedsequences/ASV_frequency_table.dynamic.format.csv'
otu_output_path <- 'OTUs'

threads <- 8
minlen <- 50

# Files to clear after processing each rank
files <- list.files("tmp", full.names = TRUE) # Files to clear after each rank

# Load data from files
taxa_file <- fread(taxa_file_path)
cutoff_file <- fread(cutoff_file_path) %>%
  select(rank = TaxonRank, taxa = TaxonomicName, cutoff = "Cut-off")
asv_sequences <- readDNAStringSet(asv_sequences_path)

# Create temporary directories for intermediate files
dir.create("tmp")
dir.create("tmp_clusters")

###############################################################################*
# (1) Match cut-offs to classified ASVs #######################################
###############################################################################*

# Match unique taxa to cut-offs
unique_taxa_cutoffs <- taxa_file %>%
  select(kingdom, phylum, class, order, family, genus, species, rank) %>%
  pivot_longer(
    -"rank", 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  filter(
    !taxa %in% c("unidentified", "") & !is.na(taxa)
  ) %>%
  select(-level) %>%
  unique(.) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  print(.)

# Define ranks and superranks for looping
ranks <- c("species", "genus", "family", "order", "class", "phylum")
superranks_species <- c("genus", "family", "order", "class", "phylum", "kingdom")
superranks_genus <- c("family", "order", "class", "phylum", "kingdom")
superranks_family <- c("order", "class", "phylum", "kingdom")
superranks_order <- c("class", "phylum", "kingdom")
superranks_class <- c("phylum", "kingdom")
superranks_phylum <- c("kingdom")

taxa_cutoffs <- get_species_cutoff(
  taxa_file, unique_taxa_cutoffs, cutoff_file, superranks_species
) %>%
  get_genus_cutoff(., unique_taxa_cutoffs, cutoff_file, superranks_genus) %>%
  get_family_cutoff(., unique_taxa_cutoffs, cutoff_file, superranks_family) %>%
  get_order_cutoff(., unique_taxa_cutoffs, cutoff_file, superranks_order) %>%
  get_class_cutoff(., unique_taxa_cutoffs, cutoff_file, superranks_class) %>%
  get_phylum_cutoff(., unique_taxa_cutoffs, cutoff_file, superranks_phylum) %>%
  group_by(kingdom) %>%
  mutate(kingdom_cutoff = phylum_cutoff) %>%
  ungroup(.)

# ###############################################################################*
# # (2) Kingdom and phylum clusters #############################################
# ###############################################################################*
#
# ### Kingdom clusters ####
#
this_rank <- "kingdom"
this_subrank <- "phylum"
this_supertaxon <- "Fungi"
#
# # Write the taxonomy file as the kingdom file:
# # Here I add a kingdom cuttoff column based on the phylum values as part of the
# # temporary solution
taxa_cutoffs %>%
  group_by(!!sym(this_rank)) %>%
  mutate(
    !!paste0(this_rank, "_cutoff") := first(!!sym(paste0(this_subrank, "_cutoff")))
  ) %>%
  ungroup(.) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# # Kingdom cut-offs
this_cutoff <- taxa_cutoffs %>%
  group_by(!!sym(this_rank)) %>%
  mutate(
    !!paste0(this_rank, "_cutoff") := first(!!sym(paste0(this_subrank, "_cutoff")))
  ) %>%
  ungroup(.) %>%
  pull(!!paste0(this_rank, "_cutoff")) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

repeat {
  #### Create cluster cores ####
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

  # Attempt to cluster kingdom fungi to subkingdom groups
  identified_asvs <- taxa_cutoffs %>%
    filter(get(this_rank) != "Fungi" & get(this_rank) != "unidentified" & get(this_rank) != "")  %>%
    select(OTU_ID) %>%
    pull(.)

  # Skip the loop for a given taxon if there are no identified ASVs
  initial_identified_count <- length(identified_asvs)
  if (initial_identified_count == 0) {
    break
  }

  unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_rank) == "Fungi" | get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # Skip the loop for a given taxon if there are no unidentified ASVs
  initial_unidentified_count <- length(unidentified_asvs)
  if (initial_unidentified_count == 0) {
    break
  }

  # Filter identified and unidentified sequences
  identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
  unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

  # Write the fasta files
  writeXStringSet(identified_sequences, "./tmp/identified_fasta")
  writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

  ## Cluster unidentified ASVs ####

  message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

  # Build the BLAST database using identified ASVs
  system2("makeblastdb",
          args = c(
            "-in", "./tmp/identified_fasta",
            "-dbtype", "nucl"
          )
  )

  # BLAST unidentified ASVs against the cluster cores
  system2("blastn",
          args = c(
            "-query", "./tmp/unidentified_fasta",
            "-db", "./tmp/identified_fasta",
            "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
            "-task", "blastn-short",
            "-num_threads", as.character(threads),
            "-out", "./tmp/unidentified_blast.out"
          )
  )

  # Read and filter BLAST results, and update taxonomy
  new_clusters <- fread("./tmp/unidentified_blast.out") %>%
    select(
      OTU_ID = V1,
      reference_ID = V2,
      sim = V3,
      alignment_length = V4,
      qlen = V5,  # Full query sequence length
      slen = V6   # Full reference sequence length
    ) %>%
    group_by(OTU_ID) %>%
    dplyr::slice(1) %>%
    ungroup(.) %>%
    # Compute score and coverage for both reference and query
    mutate(
      # Query sequence coverage
      query_coverage = case_when(
        qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
        TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
      ),
      # Reference sequence coverage
      reference_coverage = case_when(
        slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
        TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
      ),
      score = sim / 100,
      # Adjust the score for short sequence overlap following the dnabarcoder approach
      score = case_when(
        alignment_length < minlen ~ (score * alignment_length) / minlen,
        TRUE ~ score
      )
    ) %>%
    # Ensure either query or reference sequence coverage >= 90%
    filter(score >= this_cutoff & (query_coverage >= 90 | reference_coverage >= 90)) %>%
    # Initialise taxonomy columns with placeholders (including dynamic rank)
    mutate(
      phylum = "unidentified",
      class = "unidentified",
      order = "unidentified",
      family = "unidentified",
      genus = "unidentified",
      species = "unidentified",
      rank = this_rank,
      cutoff = this_cutoff,
      score = score
    ) %>%
    # Join based on reference_ID in new assignments
    left_join(
      taxa_cutoffs %>% select(
        reference_ID = OTU_ID, !!sym(this_rank),
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      ),
      by = "reference_ID"
    ) %>%
    # Select only the relevant columns for the output
    select(
      OTU_ID, reference_ID,
      kingdom, phylum, class, order, family, genus, species,
      rank, cutoff, score,
      paste0(this_rank, "_cutoff"),
      paste0(this_subrank, "_cutoff")
    )

  # Update taxonomy
  taxa_cutoffs <- bind_rows(
    new_clusters,
    taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
  )

  # Check if the number of unidentified ASVs has changed
  remaining_unidentified_count <- taxa_cutoffs %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "" | get(this_subrank) == "Fungi_phy_Incertae_sedis") %>%
    select(OTU_ID) %>%
    pull(.) %>%
    length(.)

  # Write the updated taxonomy
  fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t") %>%
  unique(.)
  ## Break or repeat the process #####
  # Break the loop if there are no more unidentified ASVs
  if (remaining_unidentified_count == 0) {
    message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
    break
    # Break the loop if no new assignments on a repeated clustering round
  } else if (remaining_unidentified_count == initial_unidentified_count) {
    message(paste0("No new assignments were made for ", this_supertaxon, "..."))
    break
    # Otherwise, continue clustering until one of the above conditions is met
  } else {
    message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
  }

}

# If for some reason there any unidentified seqences at rank kingdom, write them
# to a file for external processing

# Write the phylum file with undentified ASVs
fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt")) %>%
  filter(get(this_rank) == "unidentified" | get(this_rank) == "" | get(this_subrank) == "Fungi_phy_Incertae_sedis") %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_unidentified.txt"), sep = "\t")
# Write the phylum file with identified clusters
fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt")) %>%
  filter(get(this_rank) != "unidentified" & get(this_rank) != "" & get(this_subrank) != "Fungi_phy_Incertae_sedis") %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")


### Phylum clusters ####

# !!! Temporary solution for clustering to sub-kingdom groups finishes here !!! #

# !!! Start phylum clustering here in the future !!! #

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "kingdom"
this_rank <- "phylum"
this_subrank <- "class"

# Write the kingdom file as the phylum file: Update taxa and cutoffs
taxa_match_phylum(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_phylum_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Kingdom cut-offs for phylum clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))


### (2a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        alignment_length = V4,
        qlen = V5,  # Full query sequence length
        slen = V6   # Full reference sequence length
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score and coverage for both reference and query
      mutate(
        # Query sequence coverage
        query_coverage = case_when(
          qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
          TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
        ),
        # Reference sequence coverage
        reference_coverage = case_when(
          slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
          TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
        ),
        score = sim / 100,
        # Adjust the score for short sequence overlap following the dnabarcoder approach
        score = case_when(
          alignment_length < minlen ~ (score * alignment_length) / minlen,
          TRUE ~ score
        )
      ) %>%
      # Ensure either query or reference sequence coverage >= 90%
      filter(score >= this_cutoff & (query_coverage >= 90 | reference_coverage >= 90)) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        kingdom = this_supertaxon,
        class = "unidentified",
        order = "unidentified",
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        kingdom, phylum, class, order, family, genus, species,
        rank, cutoff, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      ) %>%
      unique(.)

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (2b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t") %>%
    unique(.)
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (3) Class clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "phylum"
this_rank <- "class"
this_subrank <- "order"

# Write the phylum file as the class file: Update taxa and cutoffs
taxa_match_class(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_class_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Stop if there are any NA values for a given taxon
if (any(is.na(taxa_cutoffs[["class_cutoff"]]))) {
  taxa_cutoffs %>% filter(is.na(class_cutoff)) %>% print(.)
  stop("Error: class_cutoff contains NA values. Exiting...")
}

# Phylum cut-offs for class clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (3a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        alignment_length = V4,
        qlen = V5,  # Full query sequence length
        slen = V6   # Full reference sequence length
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score and coverage for both reference and query
      mutate(
        # Query sequence coverage
        query_coverage = case_when(
          qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
          TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
        ),
        # Reference sequence coverage
        reference_coverage = case_when(
          slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
          TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
        ),
        score = sim / 100,
        # Adjust the score for short sequence overlap following the dnabarcoder approach
        score = case_when(
          alignment_length < minlen ~ (score * alignment_length) / minlen,
          TRUE ~ score
        )
      ) %>%
      # Ensure either query or reference sequence coverage >= 90%
      filter(score >= this_cutoff & (query_coverage >= 90 | reference_coverage >= 90)) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        phylum = this_supertaxon,
        order = "unidentified",
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        phylum, class, order, family, genus, species,
        rank, cutoff, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      ) %>%
      unique(.)

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (3b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (4) Order clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "class"
this_rank <- "order"
this_subrank <- "family"

# Write the class file as the order file: Update taxa and cutoffs
taxa_match_order(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_order_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Stop if there are any NA values for a given taxon
if (any(is.na(taxa_cutoffs[["order_cutoff"]]))) {
  taxa_cutoffs %>% filter(is.na(order_cutoff)) %>% print(.)
  stop("Error: order_cutoff contains NA values. Exiting...")
}

# Class cut-offs for order clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (4a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        alignment_length = V4,
        qlen = V5,  # Full query sequence length
        slen = V6   # Full reference sequence length
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score and coverage for both reference and query
      mutate(
        # Query sequence coverage
        query_coverage = case_when(
          qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
          TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
        ),
        # Reference sequence coverage
        reference_coverage = case_when(
          slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
          TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
        ),
        score = sim / 100,
        # Adjust the score for short sequence overlap following the dnabarcoder approach
        score = case_when(
          alignment_length < minlen ~ (score * alignment_length) / minlen,
          TRUE ~ score
        )
      ) %>%
      # Ensure either query or reference sequence coverage >= 90%
      filter(score >= this_cutoff & (query_coverage >= 90 | reference_coverage >= 90)) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        class = this_supertaxon,
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        class, order, family, genus, species,
        rank, cutoff, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (4b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (5) Family clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "order"
this_rank <- "family"
this_subrank <- "genus"

# Write the order file as the family file: Update taxa and cutoffs
taxa_match_family(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_family_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Stop if there are any NA values for a given taxon
if (any(is.na(taxa_cutoffs[["family_cutoff"]]))) {
  taxa_cutoffs %>% filter(is.na(family_cutoff)) %>% print(.)
  stop("Error: family_cutoff contains NA values. Exiting...")
}

# Order cut-offs for family clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (5a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        alignment_length = V4,
        qlen = V5,  # Full query sequence length
        slen = V6   # Full reference sequence length
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score and coverage for both reference and query
      mutate(
        # Query sequence coverage
        query_coverage = case_when(
          qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
          TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
        ),
        # Reference sequence coverage
        reference_coverage = case_when(
          slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
          TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
        ),
        score = sim / 100,
        # Adjust the score for short sequence overlap following the dnabarcoder approach
        score = case_when(
          alignment_length < minlen ~ (score * alignment_length) / minlen,
          TRUE ~ score
        )
      ) %>%
      # Ensure either query or reference sequence coverage >= 90%
      filter(score >= this_cutoff & (query_coverage >= 90 | reference_coverage >= 90)) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        order = this_supertaxon,
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        order, family, genus, species,
        rank, cutoff, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (5b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (6) Genus clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "family"
this_rank <- "genus"
this_subrank <- "species"

# Write the family file as the genus file: Update taxa and cutoffs
taxa_match_genus(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_genus_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Stop if there are any NA values for a given taxon
if (any(is.na(taxa_cutoffs[["genus_cutoff"]]))) {
  taxa_cutoffs %>% filter(is.na(genus_cutoff)) %>% print(.)
  stop("Error: genus_cutoff contains NA values. Exiting...")
}

# Family cut-offs for genus clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (6a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        alignment_length = V4,
        qlen = V5,  # Full query sequence length
        slen = V6   # Full reference sequence length
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score and coverage for both reference and query
      mutate(
        # Query sequence coverage
        query_coverage = case_when(
          qlen > 0 ~ (alignment_length / qlen) * 100,   # When qlen is valid
          TRUE ~ 0                                      # Default to 0 when qlen is invalid or zero
        ),
        # Reference sequence coverage
        reference_coverage = case_when(
          slen > 0 ~ (alignment_length / slen) * 100,   # When slen is valid
          TRUE ~ 0                                      # Default to 0 when slen is invalid or zero
        ),
        score = sim / 100,
        # Adjust the score for short sequence overlap following the dnabarcoder approach
        score = case_when(
          alignment_length < minlen ~ (score * alignment_length) / minlen,
          TRUE ~ score
        )
      ) %>%
      # Ensure both query or reference sequence coverage >= 90%
      filter(score >= this_cutoff & (query_coverage >= 90 & reference_coverage >= 90)) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        family = this_supertaxon,
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        family, genus, species,
        rank, cutoff, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (6b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-b", "T",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (7) Species clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "genus"
this_rank <- "species"

# Write the genus file as the species file: Update taxa and cutoffs
taxa_match_species(this_superrank) %>%
  # Update cut-offs
  mutate(
    !!sym(paste0(this_rank, "_cutoff")) := case_when(
      is.na(!!sym(paste0(this_rank, "_cutoff"))) ~ get_species_cutoff(cur_data(), unique_taxa_cutoffs, cutoff_file, !!sym(paste0("superranks_", this_rank))) %>% pull(!!sym(paste0(this_rank, "_cutoff"))),
      TRUE ~ !!sym(paste0(this_rank, "_cutoff"))
    )
  ) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt")) %>%
  unique(.)

# Stop if there are any NA values for a given taxon
if (any(is.na(taxa_cutoffs[["species_cutoff"]]))) {
  taxa_cutoffs %>% filter(is.na(species_cutoff)) %>% print(.)
  stop("Error: species_cutoff contains NA values. Exiting...")
}

# Genus cut-offs for species clustering
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (7a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt")) %>%
    unique(.)
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]

  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))

  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }

    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)

    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }

    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]

    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")

    ## Cluster unidentified ASVs ####

    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))

    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )

    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length qlen slen'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )

    if (file.size("./tmp/unidentified_blast.out") == 0) {
      message("No BLAST results found; skipping reference-based clustering for ", this_supertaxon)
      break
    } else {
      # Continue with reading and processing the BLAST results
      new_clusters <- fread("./tmp/unidentified_blast.out") %>%
        select(
          OTU_ID = V1,
          reference_ID = V2,
          sim = V3,
          alignment_length = V4,
          qlen = V5,  # Full query sequence length
          slen = V6   # Full reference sequence length
        ) %>%
        group_by(OTU_ID) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        # Compute score and coverage for both reference and query
        mutate(
          query_coverage = case_when(
            qlen > 0 ~ (alignment_length / qlen) * 100,
            TRUE ~ 0
          ),
          reference_coverage = case_when(
            slen > 0 ~ (alignment_length / slen) * 100,
            TRUE ~ 0
          ),
          score = sim / 100,
          # Adjust the score for short sequence overlap
          score = case_when(
            alignment_length < minlen ~ (score * alignment_length) / minlen,
            TRUE ~ score
          )
        ) %>%
        # Filter for alignment score and coverage meeting the cutoff criteria
        filter(score >= this_cutoff & (query_coverage >= 95 & reference_coverage >= 95)) %>%
        # Initialize taxonomy columns for the new clusters
        mutate(
          genus = this_supertaxon,
          rank = this_rank,
          cutoff = this_cutoff,
          score = score
        ) %>%
        # Join on reference_ID to include additional taxonomy info
        left_join(
          taxa_cutoffs %>% select(
            reference_ID = OTU_ID, !!sym(this_rank),
            paste0(this_rank, "_cutoff")
          ),
          by = "reference_ID"
        ) %>%
        # Select only the necessary columns for the final output
        select(
          OTU_ID, reference_ID,
          genus, species,
          rank, cutoff, score,
          paste0(this_rank, "_cutoff")
        ) %>%
        unique(.)

    }

    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )

    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)

    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")

    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }

  ### (7b) De novo clustering ###########################################

  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)

  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))

    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")

    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-L", "0.95",
      "-b", "T",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))

    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")

    # Initialise pseudo_id
    pseudo_id <- 1

    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name) %>%
      unique(.)

    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (8) Merge clusters ###########################################################
###############################################################################*
source("../../UK_SURVEY/ANALYSIS_SCRIPTS/merge_clusters.R")

#### (8a) Merge the clusters ####
taxonomy <- merge_clusters() %>%
  # Replace sub-kingdom groups to "Fungi": Because I used sub-kingdom groups for
  # determining global cutoffs I will rename these groups back to Fungi: 
  mutate(
    across(
      everything(),
      ~ str_replace_all(., "dikarya_fungi|terrestrial_basal_fungi|zoosporic_basal_fungi", "Fungi")
    )
  ) %>%
  print(.)

#### (8b) Generate OTU table ####
OTU_table <- fread(asv_table_path) %>% 
  select(ASV_ID = OTU_ID, everything()) %>%
  pivot_longer(
    cols = -ASV_ID,
    names_to = "sample_ID",
    values_to = "ASV_abundance"
  ) %>%
  left_join(
    taxonomy %>%
      select(ASV_ID, OTU_ID),
    by = "ASV_ID",
    relationship = "many-to-many"
  ) %>%
  group_by(OTU_ID, sample_ID) %>%
  summarise(
    OTU_abundance = sum(ASV_abundance)
  ) %>%
  ungroup(.) %>%
  pivot_wider(
    names_from = sample_ID,
    values_from = OTU_abundance,
    values_fill = list(OTU_abundance = 0)
  ) %>%
  arrange(desc(rowSums(select(., -OTU_ID)))) %>%
  print(.)

# !!! DO BETTER HERE !!!
# Seven OTUs are in the taxonomy table but not in the OTU table, and two OTUs 
# are replicated. I need to figure out how this can happen in the pipeline,
# but will continue for now by filtering the seven OTUs before filtering  and summing the 
# for the replicated OTUs
# Check OTU_ID lengths match:
taxonomy %>%
  select(-ASV_ID) %>%
  unique(.) %>%
  nrow(.)
taxonomy %>%
  select(OTU_ID) %>%
  unique(.) %>%
  nrow(.)
OTU_table %>%
  select(OTU_ID) %>%
  unique(.) %>%
  nrow(.)
# Check which OTUs are in the taxonomy file but not the OTU file
taxonomy_updated <- taxonomy %>%
  filter(OTU_ID %in% OTU_table$OTU_ID)
# Still two replicated OTUs that I deal with when saving to not complicate the 
# fasta filtering
taxonomy_updated %>%
  select(-ASV_ID) %>%
  unique(.) %>%
  nrow(.)


#### (8c) Filter fasta file ####

# Select the representative sequence for each OTU based on ASV abundance
fasta_headers <- taxonomy_updated %>%
  select(
    OTU_ID, ASV_ID, 
    kingdom, phylum, class, order, family, genus, species
  ) %>%
  unique(.) %>%
  # I can have taxonomy in sequence headers too using "header" column. This is 
  # here for reference but for now I leave as OTU_ID (see step 3 below to change)
  mutate(
    taxonomy = str_c(kingdom, phylum, class, order, family, genus, species, sep = ";"),
    header = paste0(OTU_ID, " ", taxonomy)
  ) %>%
  select(OTU_ID, header) %>%
  unique(.) %>%
  print(.)

# Generate the final fasta file:
# Step 1: Read the FASTA file
fasta <- readDNAStringSet(asv_sequences_path)

# Step 2: Filter the FASTA sequences by to representative sequences
fasta_tibble <- tibble(
  OTU_ID = names(fasta),
  sequence = as.character(fasta)
) %>%
  unique(.) %>%
  inner_join(
    fasta_headers,
    by = "OTU_ID"
  ) %>%
  print(.)

# Step 3: Write the fatsa file
filtered_fasta <- DNAStringSet(fasta_tibble$sequence)
# Keep OTU_ID for now: names(filtered_fasta) <- fasta_tibble$header
names(filtered_fasta) <- fasta_tibble$OTU_ID
#### (8d) Save OTU information ####

# Taxonomy table
taxonomy_updated %>%
  # Remove the ASV_ID column
  select(-ASV_ID) %>%
  unique(.) %>%
  group_by(
    OTU_ID, ReferenceID, 
    kingdom, phylum, class, order, family, genus,
    rank, score, cutoff
    ) %>%
  summarise(
    species = first(species),
    abundance = sum(as.numeric(abundance))
  ) %>%
  ungroup(.) %>%
  unique(.) %>%
  select(
    OTU_ID, ReferenceID,
    kingdom, phylum, class, order, family, genus, species,
    rank, score, cutoff, abundance
  ) %>%
  fwrite(., paste0(otu_output_path, "taxonomy.txt"), sep = "\t")
# OTU table
OTU_table %>%
  fwrite(., paste0(otu_output_path, "otu_table.txt"), sep = "\t")
# Fasta file
writeXStringSet(filtered_fasta, paste0(otu_output_path, "sequences.fasta"))

