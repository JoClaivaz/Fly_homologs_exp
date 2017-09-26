#### Fly_homologs_exp
To study effects of domain modifications in fly.

This folder is synchronized with our github repo : https://github.com/JoClaivaz/Fly_homologs_exp.git
This repo is public.

Scripts folder contains all scripts we created for data analysis.

6 pairwise species comparisons are considered (DROME vs Specie2):
DROME vs DROAN
DROME vs DROMO
DROME vs DROPS
DROME vs DROSI
DROME vs DROVI
DROME vs DROYA


###Pipeline
##Domain modification inference
*'extract_ortholog_fasta_sequence.py'
extraction of the amino acid fasta sequence in function of the considered ortholog pairs.

*'pfamscan'
inference of the domain architecture from the amino acid sequence presents in 'ortholog_DROME_Specie2'. The generated output is `ortholog_DROME_Specie2_domain`

*'DomainDiff'
inference of the domain modification amongst ortholog pairs (Carsten's script, obtained from Münster University). Take as input `ortholog_DROME_Specie2_domain`, the output is `ortholog_DROME_Specie2_domain_modifications`.

*'extract_ortholog_modification.py'
extraction of pair without domain modification (control group in analysis), complex domain modiffication with the same number of domain in both gene of a pair or more than 1 domain difference, domain modification with 1 domain difference. All considered pair have at least one domain inferred by pfamscan.

*'extract_all_ortholog_from_pfamscan.py'
parsing the output of pfamscan in function of the considered ortholog pairs, each considered ortholog has a unique pair number identifier.

*'analysis_pair_OMAvsPFAMSCAN' (shell script DROME vs DROAN as example)
inference whether a gene has at least one domain inferred by pfamscan (yes or no). 
NB: the last pair should be checked manually, doesn't work for the last pair

*'status_inference.R'
inference whether an ortholog pair have at least one domain inferred by pfamscan for each gene of the pair.

*'extract_both_status_pfamscan.py'
extract domain inference from pfamscan, for pair which are both inferred by pfamscan at least for one domain.

*'data_analysis_repetition_domain.R'
inference whether loss domain is repeted in the protein architecture (no informative) or not (informative, the pair is kept). Parsing pair with domain loss not repeted in the protein architecture.

*'extract_ortholog_domain_loss_domainDiff.py'
intersect of ortholog pair with one domain modification in N- or C- termini of the protein inferred by domainDIFF and the loss domain is not part of a repetition inferred by data_analysis_repetition_domain.R.

##Analysis
*'extract_fbgn_OMAname.py'
extract FBgn identifier corresponding to the OMA identifier, allowing the bound between analysis and domain modification inference.

*'data_expression_extraction_final_pair.py'
extraction of expression data obtained from BGee specific to ortholog pairs group (1 domain modification and control)

*'analysis_states_pair'
for each species, the script the different states available in the data expression recovered from BGee

*extraction of protein length
using EMBOSS software, as input 'ortholog_DROME_Specie2' (recovered from OMA) and as output 'DROME_Specie2_protein_length.txt'

*'edgeR_topGO_analysis.R'
for each pairwise species comparison, the applied edgeR pipeline allows the logFC differential gene expression inference by quasi-likelihood method, after library size and gene length normalization. The script store logFC value and p-value independently in different files.
analysis of the modification effect on logFC by ANOVA
topGO analysis allows the determination of gene ontology enrichment ('BP' and 'MF') in the domain modification group.
different plots for results visualization and report

*'DROME_id_extraction.py'
extraction of the different FBgn DROME identifier present in the different pairwise species comparisons datasets.

*obtention of CG identifier (not mandatory, if it isnt done 'DROME_ID_conversion.py' should be adapted)
in http://flybase.org/static_pages/downloads/IDConv.html convert FBgnXXX in CGXXXX using 'DROME_only_FBgn_names' file as input, and 'FlyBase_IDs.txt' as output

*'DROME_ID_conversion.py'
generate a table allowing the conversion between DROMEXXXXX (OMA), FBgnXXXXX (BGee) and CGXXXXX (GO term), usefull for the analysis between the different kind of data

*'pair_species_trimming.py'
to allow a multispecies analysis, the script keeps only ortholog family if the gene is present in the 1 domain modification or no modfication group, for each species

*'expression_analysis_all_species.R'
logFC analysis for the whole pairwise species comparison
with or without z-score normalization due to the utilization of different species and logFC fitted on different trend model
analysis of the modification effect on logFC by Kruskal-Wallis rank sum test
different plots for results visualization and report

*'expression_analysis_all_species_WhichSpeciesModified_topGO.R'
generate a table allowing to group the different gene in function of which species are in the modified group
topGO analysis allows the determination of gene ontology enrichment ('BP' and 'MF') in the domain modification group.

*'flybase_formatage.py'
formatage of the flybase file allowing to be read in R (removing warning/missing data). 
Previously, the file 'gene_rpkm_report_fb_2017_04.tsv' was download from ftp://ftp.flybase.net/releases/FB2017_04/precomputed_files/genes/.

*'GainLossDomain_Tspec_DROME.R'
'GainLoss step'
as gain or loss domain can't be inferred, the analysis of the effect on logFC in function of the domain length in DROME in comparison with the domain length of other species is done with this script.
analysis of the modification effect on logFC by ANOVA and different plots for results visualization and report
'Tspec step'
to see whether a DROME tissue has more domain modification in the species than other (tissue specific genes contribute more to the evolution of the specie)
data visualization
kruskal-wallis rank sum test and chi-squared test


###'Scripts' folder
##'Domain_modification_inference' subfolder
#'extract_ortholog_fasta_sequence.py'
Input: `pairwise_ortholog_DROME_Specie2`, `DROME.fa` and `Specie2.fa` obtained from http://omabrowser.org/oma/genomePW/ and http://omabrowser.org/oma/export/.
Output: `ortholog_DROME_Specie2` conatining the gene identifier and the fasta sequences of each considered pair.

#'extract_ortholog_modification.py'
Input: 'ortholog_DROME_Specie2_domain' and 'pairwise_ortholog_DROME_Specie2'
Output: 'ortholog_DROME_DROPS_domain_nomodif', 'ortholog_DROME_DROPS_domain_complex_modif' and 'ortholog_DROME_DROPS_1_domain_modif'.

#'extract_all_ortholog_from_pfamscan.py'
Input: 'pairwise_ortholog_DROME_DROYA' (from OMA) and 'ortholog_DROME_DROYA_domain' (output pfamscan)
Output: 'all_ortholog_DROME_DROYA_domain_parsed'

#'analysis_pair_OMAvsPFAMSCAN'
Input: 'pairwise_ortholog_DROME_Specie2' and 'ortholog_DROME_DROAN_domain'
Output: 'status_output_DROME_Specie2'

#'status_inference.R'
Input: 'status_output_DROME_Specie2'
Output: 'status_output_DROME_Specie2'

#'extract_both_status_pfamscan.py'
Input: 'DROME_Specie2_status.csv' and 'all_ortholog_DROME_Specie2_domain_parsed'
Output: 'all_both_orthologPair_DROME_Specie2'

#'data_analysis_repetition_domain.R'
Input: 'all_ortholog_DROME_Specie2_domain_modif_parsed'
Output: 'DROME_Specie2_domain_loss'

#'extract_ortholog_domain_loss_domainDiff.py'
Input: 'DROME_Specie2_domain_loss' and 'ortholog_DROME_Specie2_domain_modifications'
Output: 'final_pair_DROME_Specie2_domain_loss'

##'Analysis' subfolder
#'extract_fbgn_OMAname.py'
Input: 'final_pair_DROME_Specie2_domain_loss', 'final_pair_DROME_Specie2_domain_nomodif' and 'ortholog_DROME_Specie2'
Output: 'final_pair_DROME_Specie2_domain_loss_fbgn' and 'ortholog_DROME_Specie2_domain_nomodif_fbgn'

#'data_expression_extraction_final_pair.py'
Input: 'path_folder' containing expression data, recovered from http://bgee.org/, 'Species_names' for a given pairwise species comparisons, 'final_pair_DROME_Specie2_domain_loss_fbgn' and 'ortholog_DROME_Specie2_domain_nomodif_fbgn'.
Output: 'data_expression_DROME_Specie2' and 'data_expression_DROME_Specie2_nomodif'

#'analysis_states_pair'
Input: 'data_expression_DROME_specie2'

#'edgeR_topGO_analysis.R'
Input: 'data_expression_DROME_Specie2', 'data_expression_DROME_Specie2_nomodif', 'DROME_Specie2_protein_length.txt' and 'DROME_id_converter' (in order to have DROME gene ID converter, logFC table should previously be inferred and stored).
Output: 'logFC_drome_specie2_QL', 'pvalue_drome_specie2_QL', 'specie2_topGO_BP', 'specie2_topGO_MF' and 'GO_graph'

#'DROME_id_extraction.py'
Input: 'logFC_drome_specie2_QL' tables for each pairwise species, 'ortholog_DROME_DROAN_domain_nomodif_fbgn' and 'final_pair_DROME_DROAN_domain_loss_fbgn'
Output: 'DROME_only_FBgn_names'

#'DROME_ID_conversion.py'
Input: 'FlyBase_IDs.txt' and 'DROME_FBgn_names'
Output: 'DROME_id_converter'

#'pair_species_trimming.py'
Input: 'final_pair_DROME_Specie2_domain_loss_fbgn' and 'ortholog_DROME_Specie2_domain_nomodif_fbgn' (for the 6 Specie2)
Output: 'all_species_genes'

#'expression_analysis_all_species.R'
Input: 'logFC_drome_specie2_QL' (for the 6 Specie2), 'all_species_genes' and 'DROME_id_converter'

#'expression_analysis_all_species_WhichSpeciesModified_topGO.R'
Input: 'all_species_genes'

#'flybase_formatage.py'
Input: 'gene_rpkm_report_fb_2017_04.tsv'
Output: 'gene_rpkm_report_fb_2017_04_formated'

#'GainLossDomain_Tspec_DROME.R'
'GainLoss step'
Input: 'logFC_drome_specie2_QL', 'final_pair_DROME_Specie2_domain_loss' and 'ortholog_DROME_Specie2_domain' (for the 6 Specie2), 'all_species_genes' and 'DROME_id_converter'
'Tspec step'
Input: 'gene_rpkm_report_fb_2017_04_formated' 