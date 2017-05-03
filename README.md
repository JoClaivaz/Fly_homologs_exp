## Fly_homologs_exp
To study effects of domain modifications in fly.

This folder is synchronized with our github repo : https://github.com/JoClaivaz/Fly_homologs_exp.git
This repo is public.

scripts folder contains all scripts we created for data analysis.

### Scripts folder

`data_extraction.R`:
- Inputs: species and database, allowing to recover gene expression specific to tissues from different available experiments present in the considered database.
- Outputs: `.RData`immage, containing for each tissue from one experiment: Ensembl_Gene_ID_long/short_domain/long_domain/domain_number_long/domain_difference/Ensembl_Gene_ID_short/RPKM_long/RPKM_short/RPKM_short_norm.

`DataframeCreation_TSpec_GroupShift.R`
- Inputs: `.RData` output from `data_extraction.R`, output directory where the output file will be stored, threshold to define ubiquity and specific genes.
- Outputs: data frame in `.txt` format containing: Ensemble_gene_ID_short/Ensembl_gene_ID_long/num_domain_long/num_domain_modif/position_modif/tspec_short/tspec_long/maximal_tissue_expression_short/maximal_tissue_expression_long/specificity_short/specificity_long/specificity_shift.
The script computes the different tissue specificity index (tau), splits domain_difference (from `data_extraction.R`, be carefful only works for maximal two domain changes, need modifications in function of the fly data), discardes the domain modifications undefined, defines which tissues present the maximal expression, defines if the gene expression is ubiquitous or specific and which kind of shifts occured between considered homologs.


#### Mirjam_old_scripts subfolder
This folder contains all the scripts developped by Mirjam through her work.
