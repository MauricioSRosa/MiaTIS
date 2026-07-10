# MiaTIS


MiaTIS is a simple Python pipeline for analyzing potential key metabolites during bacteria-host interactions using transposon insertion sequencing (Tn-Seq) data.

>>>> HOW IT WORKS

The pipeline uses genome annotation files for both the bacteria and the host to predict the metabolic networks of both organisms and cross-references this information with bacterial gene essentiality data obtained via Tn-Seq during the host interaction. Based on this, the script identifies critical bottlenecks in metabolic pathways, specifically looking for non-essential genes whose products serve as the basis for genes essential to the interaction—making these products potential key metabolites for future in-depth studies of the specific bacteria-host interaction. 

- The pipeline operates in three stages:

  * Standardized annotation of bacterial and host genomes using EggNog Mapper
    - We recommend using the galaxy.eu platform for annotation
  * Analysis of bacterial gene essentiality during host interaction using TRANSIT (HMM method)
    - [Optional] Analysis of bacterial gene essentiality in rich media as an essentiality control
  * Analysis and integration of annotation and essentiality data to identify metabolic bottlenecks
    - KEGG search for annotated metabolic pathways (enzymes and metabolites involved)
    - Construction of the predicted metabolic network (Metabolite -> Enzyme -> Metabolite)
    - Mapping essentiality classifications onto the predicted network
    - If data is available, evaluation of genes essential in rich media to classify those that are strictly essential regardless of the condition. Distinguishing from essential genes based solely on interaction conditions
    - Removal of disproportionately abundant network nodes to eliminate biases involving energy currencies and other compounds used in multiple metabolic reactions
    - Identification of bottlenecks
    - Annotation and classification of bottlenecks (Predicted compounds and their classification)

- Required Files

  * Bacterial genome annotation file (.tabular format; obtained via EggNog Mapper)
  * Host genome annotation file (.tabular format; obtained via EggNog Mapper)
  * Bacterial gene essentiality file (.genes.txt format; obtained via TRANSIT) for interaction conditions
    - [Optional] Bacterial gene essentiality file (.genes.txt format; obtained via TRANSIT) for rich media conditions

- Required Python packages
  * pandas
  * requests
  * tkinter

>>>> LIMITATIONS AND FUTURE IMPROVEMENTS

  - The script considers only genes with an assigned KEGG Reaction in the annotation; this may exclude certain hydrolases and other generalist enzymes from the analysis
    * To be improved in future versions
  - During predicted metabolic network assembly, the script excludes overrepresented nodes based on a fixed cutoff value (default: 50) to avoid biases associated with compounds used in multiple metabolic reactions, such as energy currencies and water
    * Cutoff values ​​and methods for handling these overrepresented metabolites are to be improved in future versions
  - Initially, the script focuses primarily on metabolic bottlenecks
    * General analyses of essentiality within the interaction context—extending beyond just enzymes—are to be implemented
  - By default, the script utilizes all enzymatic reactions assigned to a gene during network assembly
    * To be improved to weigh which reactions actually proceed within the organism's network versus those that are merely possible but not necessarily part of the organism's metabolism
