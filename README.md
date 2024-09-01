# deepCNNvalid
Code accompanying the publication "Validation of genetic variants from NGS data using Deep Convolutional Neural Networks".

To find out about the utilisied architecture and training hyperparameters, you will be mainly interested in the Python scripts 
contained in the "training" subdirectory, which contains the actual neural network learning.

If you want to reproduce the results from the publication, the recommended way is to download and decompress the version archived
at (TODO: Zenodo), as it already includes binary dumps in NumPy format of all processed datasets. In this case, run
```
python ./training/5fold_cv.py
```

from shell in the root directory to replicate the results of the cross-validation,

```
python ./training/longterm_training.py
```

for the 100 split evaluation, and

```
python ./training/compare_datasets.py
```

for the evaluation on the independent held-out validation dataset from Kotani et al. (https://www.nature.com/articles/s41375-018-0253-3).

The code in this work was executed using Python 3.8.8 and relies on the following packages:
- NumPy 1.19.2
- Pandas 1.2.3
- TensorFlow 2.4.1
- Keras 2.4.3
- scikit-learn 0.24.1


If, on the other hand, you wish to rebuild the data from scratch, run 

```
sh download_and_process_data.sh
```

This will download the sequences from the Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra) and European Nucleotide Archive (https://www.ebi.ac.uk/ena/browser/)
before processing them. Be advised, that this will download approximately 1000GB of sequencing data before processing it. The required tools for data processing, whose binaries should be findable in PATH, are:
- sratoolkit v2.9.2 
- bwa-mem v0.7.15 
- picardtools v2.2.2 
- GATK v3.7 
- samtools v1.3.1 
- varscan v2.4.2 
- annovar v2015Dec14 