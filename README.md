# Drosophila neurotransmitters

This repository contains data on neurotransmitters used by different cell type in *Drosophila melanogaster*. 
It is stored in the form of a .csv file with version control through GitHub PRs.

Note, this rpeository is only for fast-acting small molecule transmitters, such as acetylcholine, GABA, glutamate, 
serotonin, dopamine, octopamine, nitric oxide, glycine, tyramine, etc. We have another repo for neuropeptide annotation:

## How to change data (but I do not use Git)
1. Download the `gt_data.csv` file and open with your favorite spreadsheet reader.
2. Add your data to the bottom of the `gt_data.csv` file, and save **as a CSV**.
3. Return to this page and click on the plus sign next to `Code`. You should see a dropdown menu, choose `Upload Files`.
4. Upload your version of the `gt_data.csv` file you modified. Make sure that it has the same name!
5. Fill in the form below.
    - Change `Add files via upload` to a more precise description of your changes. Keep it short and simple! 
    - Add more detail into the `Add an optional extended description...` form field below it.
6. Click on `Create a new branch for this commit and start a pull request`. You can keep the branch name that comes up by default.
7. You will be brought to a `Open a pull request` page. Review your changes and click on `Create pull request` once you're done.
8. Wait for review from one of the maintainers! You may be asked follow-up questions to make sure your added data matches with our expectations, so please be patient.

## How to change data (but I know how to use Git)
Clone this repository and create a new branch.
Modify the data file, commit with a meaningful name and push.
Create a PR for your branch once you're happy with it.

## About the data
We aim to collate as much data from the literature as possible, that can be linked to a neuronal cell type from a  connectoem dataset. 
At present, those datasets are FAFB-FlyWire (whole brain), HemiBrain (partial midbrain), FANC (ventral nerve cord), MANC (ventral nerve cord), optic-lobe (optic lobe), maleCNS (whole nervous system), BANC (whole nervous system) and the L1 (whole larval nervous system).

The file `gt_data.csv` contains one row per cell type + study, where the given study has identified a small molecule transmitter (or multiple transmitters) used by the given cell type. 

It uses a single cell type name, that can be linked between connectomes using the file `exdata/ell_type_cross_matching.csv`.

The data columns are:

*species* - the species name for the observation, for now this is only d. melanogaster

*region* - the gross subregion for the cell type, i.e. midbrain, optic lobes, ventral nerve cord

*hemilineage* - the hemilineage bundle to which the cell type belongs, in the nomenclature of Ito et al., 2013 and (midbrain), or (ventral nerve cord)

*cell_type* - a cell type name relevant to one of the connectomic datasets. In general, we prefer a FAFB-FlyWire (brain) or MANC (ventrla nerve cord) name.

*known_nt_source* - the name of the study from which the observation this row records has originated. Note, rows are unique combinations of cell_type and known_nt_source, so each can repeat over multiple rows is many studies look at the same cell type, or information on many cell types has been reported by the same study.
   
*known_nt_evidence* - the method used by the given study to determine transmission.

*known_nt_confidence* - an expression of how confident you are that the study has correctly identified the right transmitter for the given cell type, and how well that cell type has been matched to connectome data. Scores ~indicate:
  - 5: evidence for protein expression in the given cell type, cell type specific labelling.
  - 4: evidence for protein expression with coarser anatomical detail / reliable transcript expression usign in situ hybridisation, and ideally for which some negative data is available (different transmitter options tried per cell type) 
  - 3: identification of RNA transcipts related ot transmitter expression,
  - 2: Unreliable moirphological match to the EM / more bulk RNA sequecning / gross neuroanatomy based on immunohistochemistry 
  - 1: Genetic knockdown, e.g. RNAi of transmitter pathways / speculative morphological matches to EM
  - 0: Educated guesses at transmission based on any of the above, but lacking anatomical precision in matching to the EM. 

 *acetylcholine, glutamate, gaba, glycine,dopamine,serotonin, octopamine, tyramine, histamine, nitric oxide* - Each small molecule transmitter column contains a -1, 0 or 1. 1 = positive evidence for transmitter usage, 0 = no evidence either way for transmitter usage, -1 = negative evidence for transmitter usage. Due to the way wetlab reports are gathered and conveyed, we have very little negative data from the literaure but it is useful - and so we really encourage you to add it, if you have it!
