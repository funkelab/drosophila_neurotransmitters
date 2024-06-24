# Load libraries
library(fafbseg)
library(tidyverse)
library(readr)
library(catmaid)

# Get the data we have already built
gt.nt.orig <- readr::read_csv(file = "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters/gt_data.csv")

# Transmitters we care about
fast.nts <- c("acetylcholine", "gaba", "glutamate",
              "dopamine", "serotonin", "octopamine",
              "nitric oxide", "histamine", "tyramine", "glycine")
neg.fast.nts <- c("acetylcholine-negative", "gaba-negative", "glutamate-negative",
              "dopamine-negative", "serotonin-negative", "octopamine-negative",
              "nitric oxide-negative", "histamine-negative", "tyramine-negative", "glycine-negative")
all.fast.nts <- c(fast.nts, neg.fast.nts)

# Function to process known_nt column
filter_words <- function(input_string, words_to_keep, invert = FALSE){
  words <- unlist(strsplit(input_string, ",|, |;|; "))
  words <- gsub("^ | $","",words)
  if (invert){
    filtered_words <- words[! words %in% words_to_keep]
  }else{
    filtered_words <- words[words %in% words_to_keep]
  }
  paste(filtered_words, collapse = ", ")
}

# Query and organse flytable data, from midbrain and optic lobe tables
ft <- fafbseg::flytable_query("select _id, root_id, root_630, root_783, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, soma_x, soma_y, soma_z, side, ito_lee_hemilineage, hartenstein_hemilineage, top_nt, flow, super_class, cell_class, cell_type, hemibrain_match, hemibrain_type, malecns_type, cb_type, root_duplicated, morphology_group, known_nt, known_nt_source, notes from info")
ft$region <- 'midbrain'
ft.optic <- fafbseg::flytable_query("select * from optic")
ft.optic <- ft.optic[, intersect(colnames(ft.optic),
                                 colnames(ft))]
ft.optic <- ft.optic[!ft.optic$root_id %in% ft$root_id,]
ft.optic$region = 'optic_lobes'
ft.all <- plyr::rbind.fill(ft, ft.optic)
ft <- ft.all %>% dplyr::filter(!duplicated(root_id))
ft$species = 'adult_drosophila_melanogaster'

# Make cross typing sheet
ft.cross <- ft %>%
  tidyr::separate_longer_delim(hemibrain_type, delim = ", ") %>%
  tidyr::separate_longer_delim(hemibrain_type, delim = ",") %>%
  dplyr::mutate(cell_type = case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(hemibrain_type) ~ hemibrain_type,
    !is.na(cb_type) ~ cb_type,
    !is.na(morphology_group) ~ morphology_group,
    !is.na(malecns_type) ~ malecns_type,
    TRUE ~ cell_type
  )) %>%
  dplyr::distinct(cell_type, hemibrain_type, malecns_type, morphology_group, ito_lee_hemilineage, hartenstein_hemilineage, region) %>%
  dplyr::mutate(in_fafb = TRUE,
                in_hemibrain = !is.na(hemibrain_type),
                in_banc = 'to_be_found',
                in_mcns = !is.na(malecns_type),
                in_l1 = FALSE) %>%
  dplyr::arrange(ito_lee_hemilineage)
readr::write_csv(x = ft.cross, file = "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters/exdata/cell_type_cross_matching.csv")

# Take only the entries with a known_nt column
ft.nt <- ft %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(hemibrain_type) ~ hemibrain_type,
    !is.na(cb_type) ~ cb_type,
    !is.na(morphology_group) ~ morphology_group,
    TRUE ~ cell_type
  )) %>%
  dplyr::select(cell_type, ito_lee_hemilineage, hemibrain_type, notes, known_nt, known_nt_source, species, region) %>%
  dplyr::filter(!is.na(known_nt), !known_nt%in%c(""," ","NA","unknown")) %>%
  tidyr::separate_longer_delim(c(known_nt, known_nt_source), delim = ";") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt = filter_words(known_nt, all.fast.nts)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(known_nt!="") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt_evidence = gsub("\\(|\\)", "",
                                         regmatches(known_nt_source, gregexpr("\\((.*?)\\)", known_nt_source))[[1]])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(known_nt_source = gsub("\\(.*?\\)", "", known_nt_source),
                known_nt_source = gsub(" $", "", known_nt_source),
                known_nt_source = gsub("et al |et al,", "et al., ", known_nt_source),
                known_nt_source = gsub("^ ", "", known_nt_source),
                known_nt_source = gsub("\\)\\)", ")", known_nt_source)) %>%
  tidyr::separate_longer_delim(cell_type, delim = ", ") %>%
  # tidyr::separate_longer_delim(cell_type, delim = ",") %>%
  dplyr::arrange(cell_type,known_nt, known_nt_source, known_nt_evidence) %>%
  dplyr::distinct(species, region, cell_type, ito_lee_hemilineage, known_nt, known_nt_source, known_nt_evidence) %>%
  dplyr::rename(hemilineage=ito_lee_hemilineage) %>%
  dplyr::mutate(known_nt_confidence = dplyr::case_when(
    known_nt_evidence%in%c("immuno","immuno, MARCM","immuno, intersection","immuno, RNAi") ~ 5,
    known_nt_evidence%in%c("TAPIN","intersection","transgenic","EASI-FISH") ~ 4,
    known_nt_evidence%in%c("RT-PCR","scRNA-seq","immuno, tract based", "FISH", "MCFO") ~ 3,
    known_nt_evidence%in%c("RNAi","MARCM", "MCFO, unsure","FACS RNA-seq", "immuno, unsure") ~ 1,
    known_nt_evidence%in%c("educated guess", "scRNA-seq, unsure") ~ 0,
    TRUE ~ 2
  ))

# Add other missing data from maleCNS, not matched up yet
malecns.extra <- readr::read_csv("gt_sources/male_cns/malecns_extra.csv")
ft.nt <- plyr::rbind.fill(ft.nt,malecns.extra)

# Turn into a matrix
ft.nt.m <- ft.nt %>%
  tidyr::separate_longer_delim(known_nt, delim = ", ") %>%
  dplyr::distinct(cell_type, known_nt, known_nt_source, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_type), !is.na(known_nt)) %>%
  dplyr::mutate(value = dplyr::case_when(
    grepl("negative",known_nt) ~ -1,
    TRUE ~ 1
  )) %>%
  dplyr::mutate(known_nt = gsub("-.*","",known_nt)) %>%
  tidyr::spread(key = known_nt, value = value, fill = 0) %>%
  as.data.frame()

# Get MANC summary
malevnc:::choose_malevnc_dataset('MANC')
mc.find <- neuprintr::neuprint_search("Traced",field="status",dataset="manc:v1.0")
mc.ids <- unique(mc.find$bodyid)
mc.meta <- malevnc::manc_neuprint_meta(mc.ids)
mc.meta <- subset(mc.meta, !is.na(hemilineage))
vnc.hls <- table(mc.meta$hemilineage, mc.meta$predictedNt)
vnc.hls <- as.matrix(vnc.hls)
vnc.hls <- t(apply(vnc.hls, 1, function(row) row/sum(row)))
good.vnc.hls <- t(apply(vnc.hls, 1, function(row) any(row>0.9)))
good.hls <- rownames(vnc.hls)[good.vnc.hls]
good.hls <- setdiff(good.hls,"NA")
mvnc.nt <- data.frame()
for(ghl in good.hls){
  nt <- colnames(vnc.hls)[which.max(vnc.hls[ghl,])]
  dat <- mc.meta %>%
    dplyr::filter(predictedNt==nt, hemilineage == ghl) %>%
    dplyr::rename(known_nt = predictedNt, cell_type =  type) %>%
    dplyr::mutate(known_nt_source = 'Lacin et al. 2019',
                  known_nt_evidence = "FISH",
                  known_nt_confidence = 3,
                  region = "VNC",
                  species = "adult_drosophila_melanogaster") %>%
    dplyr::distinct(cell_type, hemilineage, known_nt, known_nt_source, known_nt_evidence, known_nt_confidence, region, species)
  mvnc.nt <- rbind(mvnc.nt, dat)
}
mvnc.nt.m <- mvnc.nt %>%
  dplyr::filter(!is.na(cell_type), !is.na(hemilineage), !is.na(known_nt)) %>%
  tidyr::separate_longer_delim(known_nt, delim = ", ") %>%
  dplyr::distinct(cell_type, known_nt, known_nt_source, .keep_all = TRUE) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::spread(key = known_nt, value = value, fill = 0)

# Get information from the larval connectome
conn <- catmaid_connection(server = "https://l1em.catmaid.virtualflybrain.org/")
# al <- catmaid::catmaid_get_annotationlist(conn=conn)
# ans <- unique(al$annotations$name)
# nt.annotations <- ans[grepl("ach|chol|gaba|ergic|npf|NPF|ser|dop|oct",ans)]
good.ans <- c(
  "cholinergic", "acetylcholine","mw cholinergic","Cholinergic","ChaT", "mw cholinergic", "mALT","Rh5", "Rh6", "KC", "KCs", "Kenyon_cell", "KC", "Kenyon cell",
  "GABAergic", "gaba", "GABA", "GAD1", "mw GABAergic", "APL", "mlALT",
  "Glutamatergic", "glutamate", "Vglut", "mw glutamatergic", "glutamatergic",
  "dopamine", "mw dopaminergic","dopaminergic","Dopaminergic",
  "Serotonergic", "serotonin", "serT", "mw serotonergic", "serotonergic", "CSD",
  "mw octopaminergic", "octopaminergic", "Octopaminergic", "octopamine", "VUM", "IAL-1",
  "mw tyraminergic", "Tyramine", "tyramine","Tyraminergic", "tyraminergic",
  "mw glycinergic", "Glycine", "glycinergic","Glycinergic", "glycine", "glyT",
  "mw nitric oxide", "NOS", "nos", "nitric oxide",
  "mw histaminergic", "histamine", 'Histamine', "histaminergic", "Histaminergic")
good.ans <- paste0("^",good.ans,"$")
l1.nts <- catmaid_query_by_annotation(good.ans, conn = conn)
papers <- catmaid_query_by_annotation("papers", conn = conn)
l1.nts.dat <- data.frame()
for(i in 1:length(l1.nts)){
  ldf <- l1.nts[[i]]
  deepannos <- subset(ldf, type!="neuron")
  if(nrow(deepannos)){
    extras <- catmaid_query_by_annotation(paste0("^", deepannos$name,"$"), conn = conn)
    extras <- do.call(rbind, extras)
    ldf <- rbind(ldf, extras)
  }
  ldf <- subset(ldf, type=="neuron")
  nt <- unique(subset(attr(l1.nts,"annotations"),id==names(l1.nts)[i])$name)
  known_nt_confidence <- 2
  if(grepl("^mw",nt)){
    known_nt_confidence <- 4
  }
  ldf$known_nt <- nt
  ldf$known_nt_source <- NA
  ldf$known_nt_evidence <- NA
  ldf$known_nt_confidence <- known_nt_confidence
  for(ii in 1:nrow(ldf)){
    annos <- catmaid_get_annotations_for_skeletons(ldf$skid[ii])
    paps <- subset(annos, annos$id %in% papers$id)
    if(!nrow(paps)){
      paps <- "Winding et al., 2024"
    }
    ldf$known_nt_source[ii] <- paste(paps$annotation, collapse =", ")
  }
  ldf <- dplyr::distinct(ldf, name, known_nt, known_nt_source, known_nt_confidence)
  l1.nts.dat <- rbind(l1.nts.dat, ldf)
}
l1.nts.df <- l1.nts.dat %>%
  dplyr::rename(cell_type = name) %>%
  dplyr::arrange(dplyr::desc(known_nt_confidence)) %>%
  dplyr::mutate(known_nt = dplyr::case_when(
    known_nt %in% c("cholinergic", "acetylcholine","mw cholinergic","Cholinergic","ChaT", "mw cholinergic", "mALT", "Rh5", "Rh6", "KC", "KCs", "Kenyon_cell", "KC", "Kenyon cell") ~ "acetylcholine",
    known_nt %in% c("Glutamatergic", "glutamate", "Vglut", "mw glutamatergic", "glutamatergic") ~ "glutamate",
    known_nt %in% c("GABAergic", "gaba", "GABA", "GAD1", "mw GABAergic", "APL", "mlALT") ~ "gaba",
    known_nt %in% c("dopamine", "mw dopaminergic","dopaminergic","Dopaminergic") ~ "dopamine",
    known_nt %in% c("Serotonergic", "serotonin", "serT", "mw serotonergic", "serotonergic", "CSD") ~ "serotonin",
    known_nt %in% c("mw octopaminergic", "octopaminergic", "Octopaminergic","octopamine", "VUM", "IAL-1") ~ "octopamine",
    known_nt %in% c("mw tyraminergic", "Tyramine", "tyramine","Tyraminergic", "tyraminergic") ~ "tyramine",
    known_nt %in% c("mw glycinergic", "Glycine", "glycinergic","Glycinergic", "glycine", "glyT") ~ "glycine",
    known_nt %in% c("mw nitric oxide", "NOS", "nos", "nitric oxide") ~ "nitric oxide",
    known_nt %in% c("mw histaminergic", "histamine", 'Histamine', "histaminergic", "Histaminergic") ~ "histamine",
    TRUE ~ NA
  )) %>% # some manual fixes
  dplyr::rowwise() %>%
  dplyr::mutate(known_nt = dplyr::case_when(
    grepl(cell_type,"picky") ~ "glutamate",
    TRUE ~ NA
  )) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cell_type, known_nt, .keep_all = TRUE) %>%
  dplyr::arrange(cell_type, known_nt) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(unitary = length(unique(known_nt))==1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(unitary|grepl("picky",cell_type)) %>%
  dplyr::select(cell_type, known_nt, known_nt_source, known_nt_confidence) %>%
  dplyr::mutate(
    species = "larval_drosophila_melanogaster",
    region = "cns",
    hemilineage = "unpublished",
  )
l1.nts.m <- l1.nts.df %>%
  dplyr::filter(!is.na(cell_type), !is.na(hemilineage), !is.na(known_nt)) %>%
  tidyr::separate_longer_delim(known_nt, delim = ", ") %>%
  dplyr::distinct(cell_type, known_nt, known_nt_source, .keep_all = TRUE) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::spread(key = known_nt, value = value, fill = 0)

# Order the data appropriately
gt.nt.df <- plyr::rbind.fill(ft.nt, mvnc.nt, l1.nts.df)
gt.nt <- plyr::rbind.fill(ft.nt.m, mvnc.nt.m, l1.nts.m)
gt.nt[is.na(gt.nt)] <- 0
column.order <-c("species", "region", "hemilineage", "cell_type",
                 "known_nt_source", "known_nt_evidence", "known_nt_confidence",
                 "acetylcholine", "glutamate", "gaba", "glycine",
                 "dopamine","serotonin", "octopamine", "tyramine", "histamine",
                 "nitric oxide")
gt.nt <- gt.nt[,column.order]
colnames(gt.nt) <- snakecase::to_snake_case(colnames(gt.nt))

# We know for certain papers, that Ach, vGlut and ChaT were all tested
all.fast.nt.papers <- c("Dolan et al., 2019","Aso et al., 2014",
                        "Cheong et al. 2023","Davis et al., 2020",
                        "Eschbach, Fushiki et al. 2020",
                        "Janelia","Lacin et al. 2019",
                        "Nern et al., 2024",
                        "Turner-Evans et al. 2020",
                        "Aso et al., 2019")
all.fast.nt.papers <- paste(all.fast.nt.papers,collapse="|")
for(i in 1:nrow(gt.nt)){
  pap <- gt.nt$known_nt_source[i]
  if(grepl(all.fast.nt.papers,pap)){
    if (any(gt.nt[i,c("acetylcholine","glutamate","gaba")]>0)){
      gt.nt[i,c("acetylcholine","glutamate","gaba")][gt.nt[i,c("acetylcholine","glutamate","gaba")]==0]=-1
    }
  }
}

# In these papers, if one monoamine was tested, they all likely were
all.moa.nt.papers <- c("Janelia",
                      "Nern et al., 2024",
                      "Turner-Evans et al. 2020",
                      "Aso et al., 2019")
all.moa.nt.papers <- paste(all.moa.nt.papers,collapse="|")
for(i in 1:nrow(gt.nt)){
  pap <- gt.nt$known_nt_source[i]
  if(grepl(all.moa.nt.papers,pap)){
    if (any(gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")]>0)){
      gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")][gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")]==0]=-1
    }
  }
}

# Combine old and new
gt.nt.new <- rbind(gt.nt.orig, gt.nt) %>%
  dplyr::distinct()
# dupes <- duplicated(paste0(gt.nt.new$cell_type,gt.nt.new$known_nt_source))

# Leave unknown hemilinegaes blank
gt.nt.new$hemilineage[is.na(gt.nt.new$hemilineage)] <- ""
gt.nt.new$known_nt_source[is.na(gt.nt.new$known_nt_source)] <- ""
gt.nt.new$known_nt_evidence[is.na(gt.nt.new$known_nt_evidence)] <- ""
gt.nt.new$known_nt_evidence[is.na(gt.nt.new$known_nt_confidence)] <- 0

# Save data
readr::write_csv(x = gt.nt.df,
                 file = "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters/gt_sources/bates_2024/202405-starting_gt_data.csv")
readr::write_csv(x = gt.nt.new,
                 file = "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters/gt_data.csv")
readr::write_csv(x = l1.nts.df,
                 file = "/Users/GD/LMBD/Papers/synister/drosophila_neurotransmitters/gt_sources/bates_2024/202405-starting_larval_gt_data.csv")








