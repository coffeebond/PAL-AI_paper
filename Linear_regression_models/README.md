# Multiple linear regression models for predicting poly(A) tail-length change

1.  Generate kmer metrics with `Kmer_feature_extraction.R` 

The global variable `fa` define the input fasta file, which contains the 3' UTR sequences of all mRNAs.
Use global variables `kmer_len_min` and `kmer_len_max` for setting lengths of kmers included.
Use global variable `len_max` to set the 3'-UTR length.

This will genereate two `gz` files, with one containing the kmer count metric and the other containing the kmer positional metric.
Two files generated with `kmer_len_min = 3`, `kmer_len_max = 7`, and `len_max = 300` were provided. 

Configure the feature file `../Data/Linear_regression/feature_file_list.txt` as the input for running the linear regression.

2.  Tune hyperparameters with Optuna (please use the same virtual environment as PAL-AI). 

```         
python -u Mutivariate_linear_regression.py \
           -i ../Data/Frog_oocytes/XL_oocyte_pp_7h_tail_length_change_tag_cutoff_10.txt \
           -u ../Data/Linear_regression/XENLA_10.1_GCF_XBmodels_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa \
           -f ../Data/Linear_regression/feature_file_list.txt \
           -o SGD_L_300_3_7mer_opti \
           -m opti \
           -op 1000
```

3.  Perform linear regression

```         
python -u Mutivariate_linear_regression.py \
           -i ../Data/Frog_oocytes/XL_oocyte_pp_7h_tail_length_change_tag_cutoff_10.txt \
           -u ../Data/Linear_regression/XENLA_10.1_GCF_XBmodels_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa \
           -f ../Data/Linear_regression/feature_file_list.txt \
           -o SGD_L_300_3_7mer_opti \
           -m cv
```
