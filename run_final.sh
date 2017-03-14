main(){
#    filter_CDS
#    merge_fold_change
#    hierachy
#    gen_test_group
#    gen_group_pop
#    compute_enrich_go
#    merge_go_srna_cds
#    gen_output_table_figure
    compute_percent
#    print_all
############ ranking sRNA to get the constant high expressed ones
#    pre_rank_srna
#    rank_srna
}

filter_CDS(){
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined_new.csv -s input/Staphylococcus_aureus_HG003_sRNA_new.csv -t sRNA > gene_raw_srna
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined_new.csv -s input/Staphylococcus_aureus_HG003_sRNA_new.csv -t CDS > gene_raw_cds
    python bin/filter_cds.py -i input/gene_wise_quantifications_combined_new.csv -s input/Staphylococcus_aureus_HG003_sRNA_new.csv -t both > gene_raw_read
###### After filtering or converting, run norm_deseq.R to get log2 fold change
}

merge_fold_change(){
    python bin/combine_gff_deseq.py -i output/fold_change/TSB_OD_0.2_vs_TSB_OD_0.5.csv -g input/gene_wise_quantifications_combined_new.csv -t fold_first
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
    mv fold_change.csv output/gene_wise_quantifications_combined_deseq2_new.csv
    ####### run it in cuba. hawaii scipy path is fail
    python3 bin/hierarchy.py -d 5.5 -i output/gene_wise_quantifications_combined_deseq2_new.csv -s input/Staphylococcus_aureus_HG003_sRNA_new.csv -g input/Staphylococcus_aureus_HG003_new.gff > test2
    ####### open test2 and delete first line
}

gen_test_group(){
    mkdir groups
    mkdir group_go
    mv test2 output/class_5.5_new
    python bin/get_pro_id.py -i output/class_5.5_new -g input/all_strains_uniprot_new.csv -d go-basic.obo
    mv test_* groups
    mv go_* group_go
}

gen_group_pop(){
    python bin/get_protein_id.py -t output/class_7.5_new -g input/all_strains_uniprot_new.csv
    mv all bin/goatools/data/group_pop
    mv all_go in/goatools/data/go_database
#    python bin/gen_group_pop.py -t input/all_strains_uniprot_new.csv > bin/goatools/data/group_pop
#    python bin/gen_go_database.py -t input/all_strains_uniprot_new.csv > bin/goatools/data/go_database
}

compute_enrich_go(){
	#### copy source of fisher(cfisher.so(filename may not the same), cfisher.pyx->cfisher.py and cfisher.c) to bin/goatools/scripts. And change some code from python2.7 to python3.5
    mkdir group_go_enrichment
    for FILE in $(ls groups)
    do
	python3 bin/goatools/scripts/find_enrichment.py --pval=0.05 --indent groups/$FILE bin/goatools/data/group_pop bin/goatools/data/go_database > group_go_enrichment/$FILE
    done
}

merge_go_srna_cds(){
    mkdir group_table
    for FILE in $(ls group_go_enrichment)
    do
        GO=$(echo $FILE | sed "s/test/go/")
        python bin/merge_srna_cds.py -c output/class_5.5_new -o group_go/$GO -g group_go_enrichment/$FILE -s input/Staphylococcus_aureus_HG003_sRNA_new.csv > group_table/$FILE
    done
    rm group_table/all.csv
    for FILE in $(ls group_table)
    do
	cat group_table/$FILE >> group_table/all.csv
    done
}

gen_output_table_figure(){
    mkdir group_final_output
    python3 bin/gen_final_output.py \
        -i group_table/all.csv \
        -g input/Staphylococcus_aureus_HG003_new.gff \
        -n output/gene_wise_quantifications_combined_deseq2_new.csv \
        -c output/class_5.5_new -l 3  \
        -p input/GeneSpecificInformation_NCTC8325_del.csv > group_final_output/all.csv
    mv group_*.png group_final_output/
}

compute_percent(){
    python bin/percent_go.py -i group_final_output/all.csv -s input/Staphylococcus_aureus_HG003_sRNA_new.csv > group_final_output/all_percent.csv
}

print_all(){
    for FILE in $(ls group_go)
    do
       echo "group_${FILE}" >> all.csv
       cat group_go/$FILE >> all.csv
    done
}

pre_rank_srna(){
    python bin/sort_srna_rpkm.py \
        -i input/gene_wise_quantifications_combined_rpkm_new.csv \
        -s input/Staphylococcus_aureus_HG003_sRNA_new.csv \
        > input/rpkm_name.csv
######### use excel to get the name of sort
}

rank_srna(){
    python bin/count_sort_rpkm.py -i input/sort_rpkm_name.csv > output/rank_srna_rpkm
    python3 bin/plot_srna.py -r input/rpkm_name.csv -s output/rank_srna_rpkm -l -o output/high_express_srna_log.png
    python3 bin/plot_srna.py -r input/rpkm_name.csv -s output/rank_srna_rpkm -o output/high_express_srna.png
}

main
