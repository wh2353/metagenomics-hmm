Mothur Pipeline Analysis:


The Mothur pipeline needs to be first installed:


Mothur can be found through the following link:

https://github.com/mothur/mothur/releases/tag/v1.37.3



Following the protocol in the link: http://www.mothur.org/wiki/MiSeq_SOP


Input and output files used for the analysis can be found at:

https://github.com/wh2353/metagenomics-hmm/tree/master/data


Datasets: 19 sample datasets plus 1 mock dataset as control, each dataset contains two files: which labeled as R1_001 or R2_002, thus, there are in total 19*2 + 2 = 40 fasta data files


Stability.files: txt document contains names of fasta data files, used for distinguishing among different samples


reference database: too large to be uploaded, please download through the following link: http://www.mothur.org/w/images/9/98/Silva.bacteria.zip


trainset9_032012.pds.fasta: traning datasets for naive bayesian RDP classifier for taxonomic classification


output_OTU_assingment.xlsx: output results from Mothur pipeline, including abundance and taxonomic information of all identified OTUs.


KronoaExccelTemplate-2.2.zip: visualize tool


Krona_Visualization.pdf: visualized results for abundance and taxonomic information of OTUs based on the data in output_OTU_assingment.xlsx.

