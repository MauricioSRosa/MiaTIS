# MiaTIS
A simple and fast pipeline for metabolic interaction analysis based on Tn-Seq data

MaIS is a simple pipeline for analyzing potential key metabolites during the interaction between a bacterium and its host. The pipeline uses EggNog-Mapper annotation data to predict the metabolic networks of both the bacterium and the host, and cross-references these networks with gene essentiality profiles—derived from Tn-Seq data analyzed via TRANSIT—observed during the interaction.

The pipeline requires two inputs:
1 - An EggNog-Mapper genome annotation file (in ".tabular" format).

2 - A gene essentiality file (in ".genes.txt" format) for the host interaction, analyzed via TRANSIT (HMM method).


3 - [Optional] A gene essentiality file (in ".genes.txt" format) for growth in rich medium, analyzed via TRANSIT (HMM method), to serve as a control.
