main(){

    correlation_one
#    correlation_all
#    network
}

correlation_one(){
#    mkdir query
#    python3 bin/corr_query_old.py -pc 0.77 -nc -0.77 -p 1832755:1832968:+ -g input/GeneSpecificInformation_NCTC8325_del.csv -i output/gene_wise_quantifications_combined_deseq2_new.csv
     python3 bin/corr_query.py -go bin/goatools/data/go_database \
                               -gp bin/goatools/scripts/find_enrichment.py \
                               -ga bin/goatools/data/go_database \
                               -po bin/goatools/data/group_pop \
                               -gb go-basic.obo \
                               -pc 0.77 -nc -0.77 -p 2244481:2244535:- \
                               -g input/GeneSpecificInformation_NCTC8325_del.csv \
                               -i output/gene_wise_quantifications_combined_deseq2_final.csv
#    mv *_positive* query
#    mv *_negative* query
}

correlation_all(){
    mkdir query
#    python3 bin/corr_srna_old.py -pc 0.8 -nc -0.8 -g input/GeneSpecificInformation_NCTC8325_del.csv -i output/gene_wise_quantifications_combined_deseq2_new.csv -t input/Staphylococcus_aureus_HG003_overlap.csv 
    python3 bin/corr_srna.py -go bin/goatools/data/go_database \
                             -gp bin/goatools/scripts/find_enrichment.py \
                             -ga bin/goatools/data/go_database \
                             -po bin/goatools/data/group_pop \
                             -gb go-basic.obo \
                             -pc 0.77 -nc -0.77 \
                             -g input/GeneSpecificInformation_NCTC8325_del.csv \
                             -i output/gene_wise_quantifications_combined_deseq2_final.csv \
                             -t input/Staphylococcus_aureus_HG003_overlap.csv \
		             -o query
}

network(){
    python3 bin/network_bokeh.py -i output/gene_wise_quantifications_combined_deseq2_final.csv
    python3 bin/network_bokeh2.py -i output/gene_wise_quantifications_combined_deseq2_final.csv
}

main
