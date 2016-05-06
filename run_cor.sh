main(){
#    go_base_correlation
#    CDS_base_correlation
#    merge_fold_change
    hierachy
#    go_table
}

go_base_correlation(){
    python bin/convert_go.py -i input/all_strains_uniprot.csv -g input/go.obo -s input/goslim_generic.obo > input/goslim2CDS.csv
    python3 bin/go_sRNA_cor_1.py -g input/goslim2CDS.csv -i output/gene_wise_quantifications_combined_cpm_notrrna_srna_separate.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv > test2
}

CDS_base_correlation(){
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -t sRNA > gene_raw_srna
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -t CDS > gene_raw_cds
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -t both > gene_raw_read
####### After filtering or converting, run norm_deseq.R to get log2 fold change
}

merge_fold_change(){
    python bin/combine_gff_deseq.py -i output/fold_change/TSB_OD_0.2_vs_TSB_OD_0.5.csv -g input/gene_wise_quantifications_combined.csv -t fold_first
    mv tmp_fold_change.csv fold_change.csv
    for FILE in TSB_OD_1 TSB_t0 TSB_t1 TSB_t2 TSB_ON pMEM_OD_0.2 pMEM_OD_0.5 pMEM_OD_1 pMEM_t0 pMEM_t1 pMEM_t2
    do
        python bin/combine_gff_deseq.py -i output/fold_change/TSB_OD_0.2_vs_${FILE}.csv -g fold_change.csv -t fold_middle
        mv tmp_fold_change.csv fold_change.csv
    done
    python bin/combine_gff_deseq.py -i output/fold_change/TSB_OD_0.2_vs_pMEM_ON.csv -g fold_change.csv -t fold_last
    mv tmp_fold_change.csv fold_change.csv
}

hierachy(){
#    cat output/gene_wise_quantifications_combined_deseq2_fold_CDS.csv > output/gene_wise_quantifications_combined_deseq2_fold_srna_cds_separate.csv
#    cat output/gene_wise_quantifications_combined_deseq2_fold_srna.csv >> output/gene_wise_quantifications_combined_deseq2_fold_srna_cds_separate.csv
    python3 bin/hierarchy.py -d 7.5 -i output/gene_wise_quantifications_combined_deseq2_fold_srna_cds_together.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -g input/Staphylococcus_aureus_HG003.gff > test2    
}

#go_table(){
#    python bin/gen_go_table.py -i test2 -g input/all_strains_uniprot.csv -o input/go.obo -s input/goslim_generic.obo -f input/Staphylococcus_aureus_HG003.gff -t goslim -c input/goslim2cds.csv -pg 0.5 -pc 0.33 -r input/Staphylococcus_aureus_HG003_sRNA.csv > test3
#}
main

