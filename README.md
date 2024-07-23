# Drosophila Neurotransmitters

This repository contains curated data on neurotransmitters used by different cell types (Bates et al. 2019) in *Drosophila melanogaster*. The data is stored in a version-controlled .csv file named `gt_data.csv`, with changes managed through GitHub Pull Requests.

![fafb_783_known_nts](https://github.com/funkelab/drosophila_neurotransmitters/blob/main/inst/images/fafb_783_known_nts.png?raw=true)

## Scope

This repository focuses on fast-acting small molecule transmitters, including:

- **Acetylcholine:** The primary excitatory neurotransmitter in the insect central nervous system, crucial for learning and memory.
- **Glutamate:** Can be either excitatory of inhibitory in the central nervous system, depending on the downstream receptor.
- **GABA (γ-Aminobutyric acid):** The main inhibitory neurotransmitter in the adult fly brain.
- **Glycine:** While less common in Drosophila, it can act as an inhibitory neurotransmitter in some neurons.
- **Dopamine**: Involved in reward, motivation, and motor control; important for learning and memory.
- **Serotonin:** Regulates various behaviors including sleep, feeding, and aggression.
- **Octopamine:** Often described as the insect analogue of norepinephrine, involved in arousal and aggression.
- **Tyramine:** Precursor to octopamine, also functions as a neurotransmitter affecting various behaviors.
- **Histamine:** Primary neurotransmitter in photoreceptors and some mechanosensory neurons.
- **Nitric oxide:** A gaseous signaling molecule involved in various physiological processes including memory formation.

For neuropeptide annotation, please refer to our separate repository: [funkelab/drosophila_neuropeptides](https://github.com/funkelab/drosophila_neuropeptides).

## How to Contribute Data

### For Git Novices

1. Download the `gt_data.csv` file and open it with your preferred spreadsheet application.
2. Add your data to the bottom of the file and save it as a CSV.
3. On the repository page, click the "+" sign next to "Code" and select "Upload Files".
4. Upload your modified `gt_data.csv` file, ensuring the filename remains unchanged.
5. Fill in the commit form:
   - Provide a concise description of your changes in the first field.
   - Add more detailed information in the "Add an optional extended description..." field.
6. Select "Create a new branch for this commit and start a pull request".
7. Review your changes on the "Open a pull request" page and click "Create pull request".
8. Wait for review from a maintainer. Be prepared to answer follow-up questions about your data.

### For Git Users

1. Clone this repository and create a new branch.
2. Modify the `gt_data.csv` file.
3. Commit your changes with a meaningful commit message and push to your branch.
4. Create a Pull Request for your branch when you're satisfied with your changes.

## Our Goal

Our goal is to collate as much data from the literature as possible, linking neurotransmitter information to neuronal cell types from connectomic datasets. Current datasets include:

- FAFB-FlyWire (whole brain)
- HemiBrain (partial midbrain)
- FANC (ventral nerve cord)
- MANC (ventral nerve cord)
- optic-lobe (optic lobe)
- maleCNS (whole nervous system)
- BANC (whole nervous system)
- L1 (whole larval nervous system)

The `gt_data.csv` file contains one row per cell type and study, where a given study has identified one or more small molecule transmitters used by a specific cell type.

Cross data set cell type mapping is given in the file: `/inst/extdata/cell_type_cross_matching.csv`

We have predicted a smaller set of transmitter in MANC, FAFB-FlyWire and Hemibrain (Eckstein and Bates et al., 2024).

## About the data

We aim to collate as much data from the literature as possible, that can be linked to a neuronal cell type from a  connectoem dataset. 
At present, those datasets are FAFB-FlyWire (whole brain), HemiBrain (partial midbrain), FANC (ventral nerve cord), MANC (ventral nerve cord), optic-lobe (optic lobe), maleCNS (whole nervous system), BANC (whole nervous system) and the L1 (whole larval nervous system).

The file `gt_data.csv` contains one row per cell type + study, where the given study has identified a small molecule transmitter (or multiple transmitters) used by the given cell type. 

It uses a single cell type name, that can be linked between connectomes using the file `exdata/cell_type_cross_matching.csv`.

The data columns are:

*species* - the species name for the observation, for now this is only d. melanogaster

*region* - the gross subregion for the cell type, i.e. midbrain, optic lobes, ventral nerve cord

*hemilineage* - the hemilineage bundle to which the cell type belongs, in the nomenclature of Ito et al., 2013 and (midbrain), or (ventral nerve cord)

*cell_type* - a cell type name relevant to one of the connectomic datasets. In general, we prefer a FAFB-FlyWire (brain) or MANC (ventral nerve cord) name.

*known_nt_source* - the name of the study from which the observation this row records has originated. Note, rows are unique combinations of `cell_type` and `known_nt_source`, so each can repeat over multiple rows if many studies look at the same cell type, or information on many cell types has been reported by the same study.
   
*known_nt_evidence* - the method used by the given study to determine transmission.

*known_nt_confidence* - an expression of how confident you are that the study has correctly identified the right transmitter for the given cell type, and how well that cell type has been matched to connectome data. Scores ~indicate:
  - 5: evidence for protein expression in the given cell type, cell type specific labelling.
  - 4: evidence for protein expression with coarser anatomical detail / reliable transcript expression using in-situ hybridisation, and ideally for which some negative data is available (different transmitter options tried per cell type) 
  - 3: identification of RNA transcripts related to transmitter expression,
  - 2: Unreliable morphological match to the EM / more bulk RNA sequencing / gross neuroanatomy based on immunohistochemistry 
  - 1: Genetic knockdown, e.g. RNAi of transmitter pathways / speculative morphological matches to EM
  - 0: Educated guesses at transmission based on any of the above, but lacking anatomical precision in matching to the EM. 

 *acetylcholine, glutamate, gaba, glycine,dopamine,serotonin, octopamine, tyramine, histamine, nitric oxide* - Each small molecule transmitter column contains a -1, 0 or 1. 1 = positive evidence for transmitter usage, 0 = no evidence either way for transmitter usage, -1 = negative evidence for transmitter usage. Due to the way wetlab reports are gathered and conveyed, there is relatively little negative data from the literature but it is useful - and so we really encourage you to add it, if you have it!

This project started in June 2024 with >900 cell types identified across 71 studies in adult D. melanogaster, perhaps ~1-5% of the central nervous system. The folder `gt_sources` contain some tables from a few of these studies, which we used to help build this data. We encourage you to submit such tables to this folder as well.

## Data overview

As of 20th July 2024, the data we have include, in the FAFB-783 dataset:

```
Synapses
-------------------------------
acetylcholine  27888474
dopamine         695904
gaba            7872711
glutamate       7095890
serotonin        320187
tyramine         283593
histamine       1103457
glycine            6430
nitric_oxide    2461618
octopamine       295673

Neurons
-------------------------------
acetylcholine  51559
dopamine        1411
gaba            8962
glutamate      11408
serotonin        102
tyramine         104
histamine      12120
glycine            7
nitric_oxide    4499
octopamine        62

Cell types
-------------------------------
acetylcholine  442
dopamine        58
gaba           199
glutamate      166
serotonin       30
tyramine        16
histamine        5
glycine          1
nitric_oxide    10
octopamine      21
```

Here we can see the synapse counts (note log scale) by confidence bin for out different ground-truth classes. The bars are in order of increasing confidence threshold (so left-most blue is Ach with confidence >= 1 and right-most blue is Ach with confidence >=5): 

![synapses_known_nt](https://github.com/funkelab/drosophila_neurotransmitters/blob/main/inst/images/synapses_known_nt.png?raw=true)

The same, for neuron counts:

![neurons_known_nt](https://github.com/funkelab/drosophila_neurotransmitters/blob/main/inst/images/neurons_known_nt.png?raw=true)

The same, for cell type counts:

![cell_types_known_nt](https://github.com/funkelab/drosophila_neurotransmitters/blob/main/inst/images/cell_types_known_nt.png?raw=true)

## Acknowledgements

*As of July 20th 2024 -- repo private* 

This data was collated by [Alexander Bates](https://as-bates.netlify.app/portfolio/) at Harvard Medical School while in the group of Prof. Rachel Wilson
It is manageed and curated together with [Diane Adjavon](https://adjavon.github.io/) in the laboratory of [Jan Funke](https://www.hhmi.org/scientists/jan-funke) at Janelia Research Campus. 
We thank [Mark Eddison](https://www.janelia.org/people/mark-eddison) at Janelia for sharing early assignment of EASI-FISH data for central complex neurons. 

If you use this collected data in your research please liaise with Alex, Diane and Jan on the appropriate ways to acknowledge this resource.

## Citations

1. Bates, A. S., Janssens, J., Jefferis, G. S., & Aerts, S. (2019). Neuronal cell types in the fly: single-cell anatomy meets single-cell genomics. Current Opinion in Neurobiology, 56, 125-134.

2. Eckstein, N., Bates, A. S., Champion, A., Du, M., Yin, Y., Schlegel, P., ... & Funke, J. (2024). Neurotransmitter classification from electron microscopy images at synaptic sites in Drosophila melanogaster. Cell, 187(10), 2574-2594.

For a complete list of references, please see our [citations.md](citations.md) file.

## License

[Include your chosen license information here]

## Contact

For questions or concerns, please open an issue in this repository or contact alexander_bates[at]hms.harvard.edu.
