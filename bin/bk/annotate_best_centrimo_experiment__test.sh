
# TEST FOR MODERN
perl annotate_best_centrimo_experiment.pl --best /home/rvdlee/JASPAR/ModERN/results/ModERN_fly/output/dm_Stat92E-GFP_Stat92E_E12-24/central_enrichment/selected_motif/dm_Stat92E-GFP_Stat92E_E12-24.501bp.fa.sites.centrimo.best --map /home/rvdlee/JASPAR/ModERN/data/ModERN_fly_peaks_experiment_map.txt --output TMP.TF_associated_modern

# TEST FOR CHIP-ATLAS
perl annotate_best_centrimo_experiment.pl --best /home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/ChIP-atlas/ChIP-atlas_results__test_run_screen/sacCer3/SRX977484_BY4741_TAF4/central_enrichment/selected_motif/SRX977484_BY4741_TAF4.501bp.fa.sites.centrimo.best --map /home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/ChIP-atlas/data/sacCer3/experiment_table/experimentList_TFs_and_others.tab --output TMP.TF_associated_chipatlas
