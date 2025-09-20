# Load libraries
library(bancr)
library(tidyverse)

# get GT
poss.nts <- c("acetylcholine","gaba","glutamate","dopamine","histamine","octopamine","serotonin")
gt.data <- readr::read_csv("gt_sources/bates_2024/202508-gt_data.csv") %>%
  dplyr::select(cell_type, neurotransmitter_verified, neurotransmitter_verified_source,
                neurotransmitter_verified_evidence, neurotransmitter_verified_confidence)

# Get BANC neurons
banc.meta <- banctable_query("SELECT root_id, supervoxel_id, position, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, manc_cell_type, hemibrain_cell_type from banc_meta")

# Join
banc.gt <- dplyr::bind_rows(
  dplyr::right_join(gt.data, banc.meta %>%
                      dplyr::filter(cell_type %in% !!gt.data$cell_type,
                                    !is.na(cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type),
                    relationship = "many-to-many",
                    by = "cell_type") %>%
  dplyr::right_join(gt.data, banc.meta %>%
                      dplyr::filter(fafb_cell_type %in% !!gt.data$cell_type,
                                    !is.na(fafb_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = fafb_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type"),
  dplyr::right_join(gt.data, banc.meta %>%
                      dplyr::filter(hemibrain_cell_type %in% !!gt.data$cell_type,
                                    !is.na(hemibrain_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = hemibrain_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type"),
  dplyr::right_join(gt.data, banc.meta %>%
                      dplyr::filter(manc_cell_type %in% !!gt.data$cell_type,
                                    !is.na(manc_cell_type)) %>%
                      dplyr::select(root_id, supervoxel_id, position, cell_type = manc_cell_type),
                    relationship = "many-to-many",
                    by = "cell_type")) %>%
  dplyr::mutate(neurotransmitter_verified = dplyr::coalesce(neurotransmitter_verified,
                                                            neurotransmitter_verified.x,
                                                            neurotransmitter_verified.y),
                neurotransmitter_verified_source = dplyr::coalesce(neurotransmitter_verified_source,
                                                                   neurotransmitter_verified_source.x,
                                                                   neurotransmitter_verified_source.y),
                neurotransmitter_verified_evidence = dplyr::coalesce(neurotransmitter_verified_evidence,
                                                                     neurotransmitter_verified_evidence.x,
                                                                     neurotransmitter_verified_evidence.y),
                neurotransmitter_verified_evidence = dplyr::coalesce(neurotransmitter_verified_evidence,
                                                                     neurotransmitter_verified_evidence.x,
                                                                     neurotransmitter_verified_evidence.y)) %>%
dplyr::select(  "root_id",
                "supervoxel_id",
                "position",
                "cell_type",
                "neurotransmitter_verified",
                "neurotransmitter_verified_source",
                "neurotransmitter_verified_evidence",
                "neurotransmitter_verified_confidence") %>%
  dplyr::filter(!is.na(root_id), !is.na(cell_type)) %>%
  dplyr::distinct()

# Save
readr::write_csv(x = banc.gt, file = "gt_sources/banc/202508-banc_gt_data.csv")

