
############################################################################
#### based on rule annotate_best_centrimo_experiment from Snakefile

# TF_Experiment_map : /home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/ChIP-atlas/data/sacCer3/experiment_table/experimentList_TFs_and_others.tab

# ## NOTE 4: the awk variables are relative to the ChIP-atlas experiment table.
# rule annotate_best_centrimo_experiment:
#     """
#     Assing the TF name to the selected motif.

#     This rule is executed for the best experiment of each dataset (experiment).
#     """
#     input:
#         tf_jaspar_map = config["TF_Experiment_map"],
#         best_exp = os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best")
#     output:
#         os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated")
#     message:
#         "; Assigning TF name to the experiment: {wildcards.TF} "
#     priority:
#         86
#     shell:
#         """
#         motif_name="$(more {input.best_exp} | cut -f1 | xargs basename | cut -d'.' -f1 | perl -lne '$_ =~ s/_m(\\d+)$/_peak-motifs_m$1/g; print $_ ')";

#         motif_folder="$(more {input.best_exp} | cut -f1 | xargs dirname | perl -lne '$_ =~ s/central_enrichment/motifs\/jaspar\/logos/gi; print $_')";

#         ##  -v: Pass variable as arguments
#         ##  
#         awk -v motif_n=$motif_name -v motif_f=$motif_folder 'NR==FNR{{mapping[$11] = $1}} NR!=FNR{{split($1, arr, "/"); split(arr[3], enc, "_"); print enc[1]"\\t"enc[3]"\\t"enc[3]"\\t"$0"\\t"motif_f"/"motif_n"_logo.png"}}' {input.tf_jaspar_map} {input.best_exp} > {output}
#         """
############################################################################


# tf_jaspar_map_fly=/home/rvdlee/JASPAR/ModERN/data/ModERN_fly_experimentList_TFs_and_others.tab


########### sacCer
TF_CURRENT=SRX977486_BY4741_TAF4

tf_jaspar_map=/home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/ChIP-atlas/data/sacCer3/experiment_table/experimentList_TFs_and_others.tab
# SRX977483       sacCer3 TAF4    Yeast strain    BY4741
# SRX977484       sacCer3 TAF4    Yeast strain    BY4741
# SRX977485       sacCer3 TAF4    Yeast strain    BY4741
# SRX977486       sacCer3 TAF4    Yeast strain    BY4741
# SRX977487       sacCer3 TAF5    Yeast strain    BY4741
# SRX977488       sacCer3 TAF5    Yeast strain    BY4741
# SRX977489       sacCer3 TAF5    Yeast strain    BY4741

best_exp=/home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/ChIP-atlas/ChIP-atlas_results__test_run_screen/sacCer3/${TF_CURRENT}/central_enrichment/selected_motif/${TF_CURRENT}.501bp.fa.sites.centrimo.best










# ########### FLY/WORM
# TF_CURRENT=ce_OP18_LIN-39_L1_WA

# ##INPUT
# tf_jaspar_map_worm=/home/rvdlee/JASPAR/ModERN/data/ModERN_worm_experimentList_TFs_and_others.tab
tf_jaspar_map=/home/rvdlee/JASPAR/ModERN/data/ModERN_fly_peaks_experiment_map.txt

# best_exp=/home/rvdlee/JASPAR/ModERN/results/ModERN_worm/output/${TF_CURRENT}/central_enrichment/selected_motif/${TF_CURRENT}.501bp.fa.sites.centrimo.best
# #         best_exp = os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best")
best_exp=/home/rvdlee/JASPAR/ModERN/results/ModERN_fly/output/dm_Stat92E-GFP_Stat92E_E12-24/central_enrichment/selected_motif/dm_Stat92E-GFP_Stat92E_E12-24.501bp.fa.sites.centrimo.best



##OUTPUT
output=~/TMP_TEST
# output=/home/rvdlee/JASPAR/ModERN/results/ModERN_worm/output/${TF_CURRENT}/central_enrichment/selected_motif/${TF_CURRENT}.501bp.fa.sites.centrimo.best.TF_associated
#         os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated")










####################CODE
motif_name="$(more $best_exp | cut -f1 | xargs basename | cut -d'.' -f1 | perl -lne '$_ =~ s/_m(\\d+)$/_peak-motifs_m$1/g; print $_ ')";

motif_folder="$(more $best_exp | cut -f1 | xargs dirname | perl -lne '$_ =~ s/central_enrichment/motifs\/jaspar\/logos/gi; print $_')";

echo $best_exp
cat $best_exp
echo $motif_name
# cat $motif_name
echo $motif_folder
# ls $motif_folder

##  -v: Pass variable as arguments
##  
awk -v motif_n=$motif_name -v motif_f=$motif_folder 'NR==FNR{{mapping[$11] = $1}} NR!=FNR{{split($1, arr, "/"); split(arr[3], enc, "_"); print enc[1]"\\t"enc[3]"\\t"enc[3]"\\t"$0"\\t"motif_f"/"motif_n"_logo.png"}}' $tf_jaspar_map $best_exp

## OUTPUT:
# SRX977486	TAF4	TAF4	ChIP-atlas_results/sacCer3/SRX977486_BY4741_TAF4/central_enrichment/SRX977486_BY4741_TAF4_m1.501bp.fa.sites.centrimo	-1.37858	ChIP-atlas_results/sacCer3/SRX977486_BY4741_TAF4/motifs/jaspar/logos/SRX977486_BY4741_TAF4_m1_logo.png

