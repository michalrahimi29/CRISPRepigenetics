# CRISPRepi
We share here CRISPRepi which is a tool used to predict CRISPR/Cas9 on-target efficiency by adding epigenetic information. CRISPRepi is a multi-layered preceptron, i.e. a cascade of fully connected layers of 64,
8 neurons, terminating with a final neuron, which outputs the on-target efficiency. The input is a fixed RNA sequence of 32nt length and four possible epigenetic marks: Chromatin accessibility, CTCF binding, H3K4me3 and DNA methylation. 
We utilized the Leenay et al. dataset to re-evaluate the contribution of flanking sequences, and epigenetic marks to on-target efficiency prediction and found that only downstream nucleotides improve prediction performance and that the epi-
genetic marks are highly informative of on-target efficiency with open-chromatin regions and H3K4me3 modification being the most informative. Moreover, we successfully show the generalizability of our trained model by 
successfully predicting the on-target efficiency to other cell types. 

# Requirements
The model is implemented with Keras 2.8.4 using Tensorflow backend and numpy 1.24.3.

# Run
To run CRISPRepi follow the instruction bellow:

### Model Training and evaluation:
