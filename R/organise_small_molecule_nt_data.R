# Load libraries
library(fafbseg)
library(tidyverse)
library(readr)
library(catmaid)

# Get the data we have already built
gt.nt.orig <- data.frame() #readr::read_csv(file = "gt_data.csv")

# Transmitters we care about
fast.nts <- c("acetylcholine", "glutamate",  "gaba", "glycine",
          "dopamine", "serotonin", "octopamine", "tyramine",
              "histamine", "nitric oxide")
neg.fast.nts <- c("acetylcholine-negative", "gaba-negative", "glutamate-negative",
              "dopamine-negative", "serotonin-negative", "octopamine-negative",
              "nitric oxide-negative", "histamine-negative", "tyramine-negative", "glycine-negative")
all.fast.nts <- c(fast.nts, neg.fast.nts)

# Function to process neurotransmitter_verified column
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
#ft <- fafbseg::flytable_query("select _id, root_id, root_630, root_783, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, soma_x, soma_y, soma_z, side, hemilineage, hartenstein_hemilineage, top_nt, flow, super_class, cell_class, cell_type, hemibrain_match, hemibrain_type, malecns_type, cb_type, root_duplicated, morphology_group, neurotransmitter_verified, neurotransmitter_verified_source, notes from info")
ft <- bancr::franken_meta()
# ft$region <- 'midbrain'
# ft.optic <- fafbseg::flytable_query("select * from optic")
# ft.optic <- ft.optic[, intersect(colnames(ft.optic),
#                                  colnames(ft))]
# ft.optic <- ft.optic[!ft.optic$root_id %in% ft$root_id,]
# ft.optic$region = 'optic_lobes'
# ft.all <- plyr::rbind.fill(ft, ft.optic)
# ft <- ft.all %>% dplyr::filter(!duplicated(root_id))
ft$species <- 'adult_drosophila_melanogaster'
ft <- bind_rows(
  ft,
  ft %>%
    filter(cell_type != hemibrain_type, !is.na(hemibrain_type)) %>%
    mutate(cell_type = hemibrain_type, dataset = "hemibrain")
)

# Check all cell types are consistent
problem_types <- ft %>%
  dplyr::filter(!is.na(neurotransmitter_verified)) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::filter(n_distinct(neurotransmitter_verified) > 1 |
           n_distinct(neurotransmitter_verified_source) > 1) %>%
  dplyr::distinct(cell_type) %>%
  dplyr::pull(cell_type)
message("Types with inconsistent labels: ", length(problem_types))

# Make cross typing sheet
ft.cross <- ft %>%
  tidyr::separate_longer_delim(hemibrain_type, delim = ", ") %>%
  tidyr::separate_longer_delim(hemibrain_type, delim = ",") %>%
  dplyr::mutate(cell_type = case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(hemibrain_type) ~ hemibrain_type,
    !is.na(morphology_group) ~ morphology_group,
    TRUE ~ cell_type
  )) %>%
  dplyr::distinct(cell_type, hemibrain_type, morphology_group, hemilineage, region) %>%
  dplyr::mutate(in_fafb = TRUE,
                in_hemibrain = !is.na(hemibrain_type),
                in_banc = TRUE,
                in_l1 = FALSE) %>%
  dplyr::arrange(hemilineage)
readr::write_csv(x = ft.cross, file = "inst/extdata/cell_type_cross_matching.csv")

# Take only the entries with a neurotransmitter_verified column
ft.nt <- ft %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(hemibrain_type) ~ hemibrain_type,
    !is.na(morphology_group) ~ morphology_group,
    TRUE ~ cell_type
  )) %>%
  dplyr::select(cell_type, hemilineage, hemibrain_type, neurotransmitter_verified, neurotransmitter_verified_source, species, region) %>%
  dplyr::filter(!is.na(neurotransmitter_verified),
                !cell_type%in%c("AN_4_None","unknown","columnar"),
                !grepl("AN_|SA_",cell_type),
                !is.na(cell_type),
                !neurotransmitter_verified%in%c(""," ","NA","unknown")) %>%
  tidyr::separate_longer_delim(c(neurotransmitter_verified, neurotransmitter_verified_source), delim = ";") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(neurotransmitter_verified = filter_words(neurotransmitter_verified, all.fast.nts)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(neurotransmitter_verified!="") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(neurotransmitter_verified_evidence = gsub("\\(|\\)", "",
                                         regmatches(neurotransmitter_verified_source, gregexpr("\\((.*?)\\)", neurotransmitter_verified_source))[[1]])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(neurotransmitter_verified_source = gsub("\\(.*?\\)", "", neurotransmitter_verified_source),
                neurotransmitter_verified_source = gsub(" $", "", neurotransmitter_verified_source),
                neurotransmitter_verified_source = gsub("et al |et al,", "et al., ", neurotransmitter_verified_source),
                neurotransmitter_verified_source = gsub("^ ", "", neurotransmitter_verified_source),
                neurotransmitter_verified_source = gsub("\\)\\)", ")", neurotransmitter_verified_source)) %>%
  tidyr::separate_longer_delim(cell_type, delim = ", ") %>%
  # tidyr::separate_longer_delim(cell_type, delim = ",") %>%
  dplyr::arrange(cell_type,neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_evidence) %>%
  dplyr::distinct(species, region, cell_type, hemilineage, neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_evidence) %>%
  dplyr::rename(hemilineage=hemilineage) %>%
  dplyr::mutate(neurotransmitter_verified_confidence = dplyr::case_when(
    neurotransmitter_verified_evidence%in%c("immuno","immuno, MARCM","immuno, intersection","immuno, RNAi") ~ 5,
    neurotransmitter_verified_evidence%in%c("TAPIN","intersection","transgenic","EASI-FISH") ~ 4,
    neurotransmitter_verified_evidence%in%c("RT-PCR","scRNA-seq","immuno, tract based", "FISH", "MCFO") ~ 3,
    neurotransmitter_verified_evidence%in%c("RNAi","MARCM", "MCFO, unsure","FACS RNA-seq", "immuno, unsure") ~ 1,
    neurotransmitter_verified_evidence%in%c("educated guess", "scRNA-seq, unsure","unknown") ~ 0,
    TRUE ~ 2
  ))

# All flywire with neurotransmitter_verified
ft.nt.all <- ft %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(hemibrain_type) ~ hemibrain_type,
    !is.na(morphology_group) ~ morphology_group,
    TRUE ~ cell_type
  )) %>%
  dplyr::filter(!is.na(neurotransmitter_verified), !neurotransmitter_verified%in%c(""," ","NA","unknown")) %>%
  dplyr::distinct(neuron_id, top_nt, cell_type, hemilineage, neurotransmitter_verified, neurotransmitter_verified_source)
readr::write_csv(x = ft.nt.all,
                 file = "gt_sources/bates_2024/202508-franken_gt_data.csv")

# Add other missing data from maleCNS, not matched up yet
malecns.extra <- readr::read_csv("gt_sources/extra.csv")
ft.nt <- plyr::rbind.fill(ft.nt,malecns.extra)

# Turn into a matrix
ft.nt.m <- ft.nt %>%
  tidyr::separate_longer_delim(neurotransmitter_verified, delim = ", ") %>%
  dplyr::distinct(cell_type, neurotransmitter_verified, neurotransmitter_verified_source, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(cell_type), !is.na(neurotransmitter_verified)) %>%
  dplyr::mutate(value = dplyr::case_when(
    grepl("negative",neurotransmitter_verified) ~ -1,
    TRUE ~ 1
  )) %>%
  dplyr::mutate(neurotransmitter_verified = gsub("-.*","",neurotransmitter_verified)) %>%
  dplyr::filter(!is.na(neurotransmitter_verified)) %>%
  tidyr::spread(key = neurotransmitter_verified, value = value, fill = 0) %>%
  as.data.frame()

# # Get MANC summary
# malevnc:::choose_malevnc_dataset('MANC')
# mc.find <- neuprintr::neuprint_search("Traced",field="status",dataset="manc:v1.0")
# mc.ids <- unique(mc.find$bodyid)
# mc.meta <- malevnc::manc_neuprint_meta(mc.ids)
# mc.meta <- subset(mc.meta, !is.na(hemilineage))
# vnc.hls <- table(mc.meta$hemilineage, mc.meta$predictedNt)
# vnc.hls <- as.matrix(vnc.hls)
# vnc.hls <- t(apply(vnc.hls, 1, function(row) row/sum(row)))
# good.vnc.hls <- t(apply(vnc.hls, 1, function(row) any(row>0.9)))
# good.hls <- rownames(vnc.hls)[good.vnc.hls]
# good.hls <- setdiff(good.hls,"NA")
# good.hls <- good.hls[!is.na(good.hls)]
# mvnc.nt <- data.frame()
# for(ghl in good.hls){
#   nt <- colnames(vnc.hls)[which.max(vnc.hls[ghl,])]
#   dat <- mc.meta %>%
#     dplyr::filter(predictedNt==nt, hemilineage == ghl) %>%
#     dplyr::rename(neurotransmitter_verified = predictedNt, cell_type =  type) %>%
#     dplyr::mutate(neurotransmitter_verified_source = 'Lacin et al. 2019',
#                   neurotransmitter_verified_evidence = "FISH",
#                   neurotransmitter_verified_confidence = 3,
#                   region = "VNC",
#                   species = "adult_drosophila_melanogaster") %>%
#     dplyr::distinct(cell_type, hemilineage, neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_evidence, neurotransmitter_verified_confidence, region, species)
#   mvnc.nt <- rbind(mvnc.nt, dat)
# }
# mvnc.nt.m <- mvnc.nt %>%
#   dplyr::filter(!is.na(cell_type), !is.na(hemilineage), !is.na(neurotransmitter_verified)) %>%
#   tidyr::separate_longer_delim(neurotransmitter_verified, delim = ", ") %>%
#   dplyr::distinct(cell_type, neurotransmitter_verified, neurotransmitter_verified_source, .keep_all = TRUE) %>%
#   dplyr::mutate(value = 1) %>%
#   tidyr::spread(key = neurotransmitter_verified, value = value, fill = 0)

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
  neurotransmitter_verified_confidence <- 2
  if(grepl("^mw",nt)){
    neurotransmitter_verified_confidence <- 4
  }
  ldf$neurotransmitter_verified <- nt
  ldf$neurotransmitter_verified_source <- NA
  ldf$neurotransmitter_verified_evidence <- NA
  ldf$neurotransmitter_verified_confidence <- neurotransmitter_verified_confidence
  for(ii in 1:nrow(ldf)){
    annos <- catmaid_get_annotations_for_skeletons(ldf$skid[ii])
    paps <- subset(annos, annos$id %in% papers$id)
    if(!nrow(paps)){
      paps <- "Winding et al., 2024"
    }
    ldf$neurotransmitter_verified_source[ii] <- paste(paps$annotation, collapse =", ")
  }
  ldf <- dplyr::distinct(ldf, name, neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_confidence)
  l1.nts.dat <- rbind(l1.nts.dat, ldf)
}
l1.nts.df <- l1.nts.dat %>%
  dplyr::rename(cell_type = name) %>%
  dplyr::arrange(dplyr::desc(neurotransmitter_verified_confidence)) %>%
  dplyr::mutate(neurotransmitter_verified = dplyr::case_when(
    neurotransmitter_verified %in% c("cholinergic", "acetylcholine","mw cholinergic","Cholinergic","ChaT", "mw cholinergic", "mALT", "Rh5", "Rh6", "KC", "KCs", "Kenyon_cell", "KC", "Kenyon cell") ~ "acetylcholine",
    neurotransmitter_verified %in% c("Glutamatergic", "glutamate", "Vglut", "mw glutamatergic", "glutamatergic") ~ "glutamate",
    neurotransmitter_verified %in% c("GABAergic", "gaba", "GABA", "GAD1", "mw GABAergic", "APL", "mlALT") ~ "gaba",
    neurotransmitter_verified %in% c("dopamine", "mw dopaminergic","dopaminergic","Dopaminergic") ~ "dopamine",
    neurotransmitter_verified %in% c("Serotonergic", "serotonin", "serT", "mw serotonergic", "serotonergic", "CSD") ~ "serotonin",
    neurotransmitter_verified %in% c("mw octopaminergic", "octopaminergic", "Octopaminergic","octopamine", "VUM", "IAL-1") ~ "octopamine",
    neurotransmitter_verified %in% c("mw tyraminergic", "Tyramine", "tyramine","Tyraminergic", "tyraminergic") ~ "tyramine",
    neurotransmitter_verified %in% c("mw glycinergic", "Glycine", "glycinergic","Glycinergic", "glycine", "glyT") ~ "glycine",
    neurotransmitter_verified %in% c("mw nitric oxide", "NOS", "nos", "nitric oxide") ~ "nitric oxide",
    neurotransmitter_verified %in% c("mw histaminergic", "histamine", 'Histamine', "histaminergic", "Histaminergic") ~ "histamine",
    TRUE ~ NA
  )) %>% # some manual fixes
  dplyr::rowwise() %>%
  dplyr::mutate(neurotransmitter_verified = dplyr::case_when(
    grepl(cell_type,"picky") ~ "glutamate",
    TRUE ~ neurotransmitter_verified
  )) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cell_type, neurotransmitter_verified, .keep_all = TRUE) %>%
  dplyr::arrange(cell_type, neurotransmitter_verified) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(unitary = length(unique(neurotransmitter_verified))==1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(unitary|grepl("picky",cell_type)) %>%
  dplyr::select(cell_type, neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_confidence) %>%
  dplyr::mutate(
    species = "larval_drosophila_melanogaster",
    region = "cns",
    hemilineage = "unpublished",
  ) %>%
  dplyr::filter(!is.na(neurotransmitter_verified))
l1.nts.m <- l1.nts.df %>%
  dplyr::filter(!is.na(cell_type), !is.na(hemilineage), !is.na(neurotransmitter_verified)) %>%
  tidyr::separate_longer_delim(neurotransmitter_verified, delim = ", ") %>%
  dplyr::distinct(cell_type, neurotransmitter_verified, neurotransmitter_verified_source, .keep_all = TRUE) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::spread(key = neurotransmitter_verified, value = value, fill = 0)

# Order the data appropriately
gt.nt.df <- plyr::rbind.fill(ft.nt, l1.nts.df) # mvnc.nt
gt.nt <- plyr::rbind.fill(ft.nt.m, l1.nts.m) # mvnc.nt.m
gt.nt[is.na(gt.nt)] <- 0
column.order <-c("species", "region", "hemilineage", "cell_type",
                 "neurotransmitter_verified_source", "neurotransmitter_verified_evidence", "neurotransmitter_verified_confidence",
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
                        "Nern et al., 2024", "Wolff et al., 2024",
                        "Turner-Evans et al. 2020",
                        "Aso et al., 2019")
all.fast.nt.papers <- paste(all.fast.nt.papers,collapse="|")
for(i in 1:nrow(gt.nt)){
  pap <- gt.nt$neurotransmitter_verified_source[i]
  if(grepl(all.fast.nt.papers,pap)){
    if (any(gt.nt[i,c("acetylcholine","glutamate","gaba")]>0)){
      gt.nt[i,c("acetylcholine","glutamate","gaba")][gt.nt[i,c("acetylcholine","glutamate","gaba")]==0]=-1
    }
  }
}

# In these papers, if one monoamine was tested, they all likely were
all.moa.nt.papers <- c("Janelia", "Wolff et al., 2024",
                      "Nern et al., 2024",
                      "Turner-Evans et al. 2020",
                      "Aso et al., 2019")
all.moa.nt.papers <- paste(all.moa.nt.papers,collapse="|")
for(i in 1:nrow(gt.nt)){
  pap <- gt.nt$neurotransmitter_verified_source[i]
  if(grepl(all.moa.nt.papers,pap)){
    if (any(gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")]>0)){
      gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")][gt.nt[i,c("dopamine","serotonin","octopamine", "tyramine")]==0]=-1
    }
  }
}

# Combine old and new
gt.nt.new <- rbind(gt.nt.orig, gt.nt) %>%
  dplyr::distinct()
# dupes <- duplicated(paste0(gt.nt.new$cell_type,gt.nt.new$neurotransmitter_verified_source))

# Leave unknown hemilinegaes blank
gt.nt.new$hemilineage[is.na(gt.nt.new$hemilineage)] <- ""
gt.nt.new$neurotransmitter_verified_source[is.na(gt.nt.new$neurotransmitter_verified_source)] <- ""
gt.nt.new$neurotransmitter_verified_evidence[is.na(gt.nt.new$neurotransmitter_verified_evidence)] <- ""
gt.nt.new$neurotransmitter_verified_evidence[is.na(gt.nt.new$neurotransmitter_verified_confidence)] <- 0

# Save data
readr::write_csv(x = gt.nt.df,
                 file = "gt_sources/bates_2024/202508-gt_data.csv")
readr::write_csv(x = gt.nt.new,
                 file = "gt_data.csv")
readr::write_csv(x = l1.nts.df,
                 file = "gt_sources/bates_2024/202508-starting_larval_gt_data.csv")

#######################################
### Make explicit hemibrain mapping ###
#######################################

hb.find <- neuprintr::neuprint_search("Traced",field="status",dataset="hemibrain:v1.2.1")
hb.ids <- unique(hb.find$bodyid)
hb.meta.orig <- neuprintr::neuprint_get_meta(hb.ids)
hb.nt <- hb.meta.orig %>%
  dplyr::filter(!is.na(type)) %>%
  dplyr::left_join(ft.nt,
                   by = c("type"="cell_type"),
                   relationship = "many-to-many")
hb.nt <- hb.nt %>%
  dplyr::filter(!is.na(neurotransmitter_verified)) %>%
  dplyr::arrange(type,bodyid) %>%
  dplyr::select(bodyid, pre, type, species, region, hemilineage, neurotransmitter_verified, neurotransmitter_verified_source, neurotransmitter_verified_evidence, neurotransmitter_verified_confidence)
readr::write_csv(x = hb.nt,
                 file = "gt_sources/bates_2024/202508-hemibrain_gt_data.csv")

#############################
### Make plots for README ###
#############################
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)

# Organise data
ft.plot <- ft %>%
  dplyr::mutate(known_fast_nts = gsub("; |:|, |,",", ",neurotransmitter_verified)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(known_fast_nts = paste(sort(filter_words(unique(unlist(strsplit(known_fast_nts, split=", "))),fast.nts)),
                                       collapse = " ")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(known_fast_nts = ifelse(grepl(",",known_fast_nts), "co-expression",known_fast_nts)) %>%
  dplyr::mutate(known_fast_nts = ifelse(is.na(known_fast_nts)|known_fast_nts=="", "unknown",known_fast_nts)) %>%
  dplyr::mutate(super_class = ifelse(is.na(super_class),flow,super_class)) %>%
  dplyr::mutate(super_class = ifelse(is.na(super_class),"other",super_class)) %>%
  dplyr::select(region, super_class, known_fast_nts)

# Read the color data from CSV
color_data <- read_csv("settings/paper_colours.csv")

# Create a named vector of colors
color_vector <- setNames(color_data$hex, color_data$label)

# Prepare the data: count, calculate percentages, and order factors
nt.order <- c("acetylcholine", "glutamate",  "gaba", "glycine",
              "dopamine", "serotonin", "octopamine", "tyramine",
              "histamine", "nitric oxide", "co-expression", "unknown")
plot_data <- ft.plot %>%
  count(super_class, known_fast_nts) %>%
  group_by(super_class) %>%
  mutate(percentage = n / sum(n),
         total_count = sum(n)) %>%
  ungroup() %>%
  mutate(known_fast_nts = factor(known_fast_nts,
                                 levels = rev(nt.order),
                                 ordered = TRUE))

# Create the plot
g <- ggplot(plot_data, aes(x = super_class, y = percentage, fill = known_fast_nts)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(y = 1.05, label = total_count, group = super_class),
            color = "black", size = 3, fontface = "bold") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values = color_vector, drop = FALSE) +
  labs(title = "FAFB-FlyWire 783 / MANC v1.2.1 neuron distribution by super class\n and known fast-acting neurotransmitter",
       x = "super class",
       y = "percentage",
       fill = "fast-acting neurotransmitter") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  coord_flip()  # Flip coordinates for horizontal bars

# Save
ggsave(g, filename = "inst/images/fafb_783_known_nts.png", width = 8, height = 8)

#############################
### Save GT data for BANC ###
#############################
library(bancr)

# get GT
poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin")
gt.data <- read_csv("gt_sources/bates_2024/202508-gt_data.csv")

# Get BANC neurons
banc.meta <- banctable_query("SELECT root_id, supervoxel_id, position, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, manc_cell_type, hemibrain_cell_type from banc_meta")

# Join
banc.gt <- gt.data %>%
  dplyr::right_join(banc.meta %>%
                      dplyr::filter(cell_type %in% !!gt.data$cell_type,
                                    !is.na(cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type),
                    relationship = "many-to-many",
                    by = "cell_type") %>%
  dplyr::right_join(banc.meta %>%
                      dplyr::filter(fafb_cell_type %in% !!gt.data$cell_type,
                                    !is.na(fafb_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = fafb_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type") %>%
  dplyr::right_join(banc.meta %>%
                      dplyr::filter(hemibrain_cell_type %in% !!gt.data$cell_type,
                                    !is.na(hemibrain_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = hemibrain_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type") %>%
  dplyr::right_join(banc.meta %>%
                      dplyr::filter(manc_cell_type %in% !!gt.data$cell_type,
                                    !is.na(manc_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = manc_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type") %>%
  dplyr::mutate(root_id = coalesce(root_id.x,root_id.x.x,root_id.y,root_id.y.y),
                supervoxel_id = coalesce(supervoxel_id.x,supervoxel_id.x.x,supervoxel_id.y,supervoxel_id.y.y),
                position = coalesce(position.x,position.x.x,position.y,position.y.y)) %>%
  dplyr::select(-root_id.x,-root_id.y,-root_id.x.x,-root_id.y.y,
                -supervoxel_id.x,-supervoxel_id.x.x,-supervoxel_id.y,-supervoxel_id.y.y,
                -position.x,-position.x.x,-position.y,-position.y.y) %>%
  dplyr::filter(!is.na(root_id), !is.na(cell_type)) %>%
  dplyr::distinct()

# Save
write_csv(x = banc.gt, file = "gt_sources/banc/202508-banc_gt_data.csv")

# Double check optic lobe types




















