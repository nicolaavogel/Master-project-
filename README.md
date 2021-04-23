# Master-project
Master project: "Hierarchical multiclass supervised classification of metabarcoded environmental DNA in a forensic context"

ABSTRACT:

Forensic ecology is a constantly developing field and provides new possibilities to find evidence
in criminal investigations. To obtain the evidence, a combination of different approaches to
analysis are necessary:
Supervised classification aims to predict categorical response variables into predefined classes.
The underlying data patterns are used to build a prediction model to identify the classes of
unknown test samples. Hierarchical multiclass supervised classification uses multiple multiclass
classifications to classify samples on different levels of the hierarchy. This method reduces
false-negative predictions from data imbalance. One of the data types used for supervised
classification models is environmental DNA (eDNA). The different detected species create a
unique pattern that can be used for predictions of environmental gradients, ecological habitats,
or the geographical origin. These mechanisms are used for forensic analysis. Metabarcoding
eDNA of ecological elements (e.g. soil, pollen, settled dust) that can be found at a crime scene
can reveal information that is usable as evidence. In forensic geoscience, the provenancing of
samples aims to find the geographical origin of a sample, matching analysis compares
discovered soil samples to soil samples from a crime scene to establish a connection, and
habitat prediction aims to find potential crime scenes.
In this thesis, the aim was to combine the most commonly used analysis tools in forensics to
obtain more accurate predictions and information that can be used for investigations. The
combination of sample matching for provenancing and habitat classification was computed over
three hierarchical levels to account for imbalanced data and improve multiclass classification
results. Simultaneously, a distance-based assessment of the classification computed a likelihood
ratio based on a Bayesian framework for a statistical assessment of the data. The training data
stems from the “Biowide” project and describes different habitats throughout Denmark. The test
data consisted of single soil core samples taken from the “Biowide” sites and urban samples
taken in the Copenhagen area. The results were mixed: The hierarchical multiclass supervised
classification with a k-nearest neighbour classifier and relative abundance read count data as
predictors showed high potential to accurately predict even low samples size classes and raise
the sensitivity of habitat classifications. All unknown test samples that were used, were classified
correctly. The distance-based assessment of the classification was difficult to compute due to the
high intra-class difference and low inter-class difference. Those difficulties, combined with the
wide-spread samples, common in ecological datasets, did not allow for a convincing likelihood
ratio approach. To overcome those difficulties, methods have to improve the inter-class and
intra-class difference, as well as a weighted multidimensional distance to account for variance
differences after ordination. A raise of the sample size to improve statistical significance is
vitally important. The matching analysis showed accurate results, but comparison pictures need
distinguishable ecological features to be helpful for investigations.
