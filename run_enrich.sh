main(){
#    gen_test_group
#    compute_enrich_go
#    compute_enrich_goslim
#    merge_go_srna_cds
    gen_output_table_figure
}

gen_test_group(){
    python bin/get_pro_id.py -i input/class_7.5_new -g input/all_strains_uniprot.csv -d go-basic.obo
    mv test_* groups
    mv go_* group_go
}

compute_enrich_go(){
    for FILE in $(ls groups)
    do
        python bin/goatools/scripts/find_enrichment.py --pval=0.05 --indent groups/$FILE bin/goatools/data/group_pop bin/goatools/data/go_database > group_go_enrichment/$FILE
    done
}

compute_enrich_goslim(){
    for FILE in $(ls groups)
    do
        python bin/goatools/scripts/find_enrichment.py --pval=0.05 --indent groups/$FILE bin/goatools/data/group_pop bin/goatools/data/goslim_database > group_goslim_enrichment/$FILE
    done
}

merge_go_srna_cds(){
    for FILE in $(ls group_go_enrichment)
    do
        GO=$(echo $FILE | sed "s/test/go/")
        python bin/merge_srna_cds.py -c input/class_7.5_new -o group_go/$GO -g group_go_enrichment/$FILE -s input/Staphylococcus_aureus_HG003_sRNA.csv > group_table/$FILE
    done
    rm group_table/all.csv
    for FILE in $(ls group_table)
    do
        cat group_table/$FILE >> group_table/all.csv
    done
}

gen_output_table_figure(){
    python bin/gen_final_output.py \
        -i group_table/all.csv \
        -g input/Staphylococcus_aureus_HG003.gff \
        -n output/gene_wise_quantifications_combined_deseq2_fold_srna_cds_together.csv \
        -c input/class_7.5_new -l 3  \
        -p input/GeneSpecificInformation_NCTC8325.tsv > group_final_output/all.csv
    mv group_*.png group_final_output/
    python bin/percent_go.py -i group_final_output/all.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv > group_final_output/all_percent.csv
#    python bin/gen_all_protein.py -c input/class_7.5_new -i group_final_output/all_percent.csv > group_final_output/all_include_nonenrich.csv
}
main
