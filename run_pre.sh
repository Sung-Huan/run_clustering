main(){

#    regulon_compute
#    kmeans
#    hierachy
#    plot_kmeans
#    gen_gff
#    get_protein_names
#    GO_term
#    pre_rank_srna
    rank_srna
#    anti_kmeans
#    anti_plot_kmeans
#    anti_get_protein_names
#    anti_GO_term
}

regulon_compute(){
    python bin/regulon.py -i input/gene_wise_quantifications_combined_rpkm.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -o output
    python bin/delete_nan.py -i output/pmem_p_all > pmem_p_all
    mv pmem_p_all output/pmem_p_all
    python bin/delete_nan.py -i output/pmem_s_all > pmem_s_all
    mv pmem_s_all output/pmem_s_all
    python bin/delete_nan.py -i output/tsb_p_all > tsb_p_all
    mv tsb_p_all output/tsb_p_all
    python bin/delete_nan.py -i output/tsb_s_all > tsb_s_all
    mv tsb_s_all output/tsb_s_all
}
## need to run R kmeans after all_rpkm.py
#data <- read.table("gene_rpkm", sep="\t")
#> kc <- kmeans(data, 40, iter=25)
#> lapply(kc[1], write, "kmeans.txt", append=TRUE)
kmeans(){
    python bin/all_rpkm.py -i output/gene_wise_quantifications_combined_cpm_notrrna_srna.csv -g input/Staphylococcus_aureus_HG003.gff -s input/Staphylococcus_aureus_HG003_sRNA.csv > gene_cpm
#    python bin/classify.py -r gene_rpkm -k kmeans.txt -n name_rpkm
#    mv kmeans_* output/kmeans_all14_normal
#    mv name_* output/kmeans_all14_normal
#    rm kmeans.txt
}

hierachy(){
    python3 bin/hierarchy.py -g input/Staphylococcus_aureus_HG003.gff \
        -s input/Staphylococcus_aureus_HG003_sRNA.csv \
        -i input/gene_wise_quantifications_combined_rpkm.csv
    mv hie* output/hierachy
    for FILE in $(ls output/hierachy/*.gff)
    do
        cp ${FILE}.gff ANNOgesic/output/target/annotation
        sh anno.sh
        INDEX=$(echo $FILE | cut -d'_' -f 5)
#        mkdir output/kmeans_all14_normal/fig_kmeans_$INDEX
        mv ANNOgesic/output/Go_term/all_CDS/statistics/Staphylococcus_aureus_HG003/stat* output/hierachy/hie_$INDEX
    done
}

plot_kmeans(){
    rm output/kmeans_all14_normal/*.png 
    for FILE in $(ls output/kmeans_all14_normal/kmeans_*)
    do
        NAME=$(echo $FILE | sed "s/al\/kmeans/al\/name_kmeans/")
        echo $FILE
        echo $NAME
        python bin/plot_all_kmeans.py -i $FILE -n $NAME -f all -s input/Staphylococcus_aureus_HG003_sRNA.csv
    done
    mv *.png output/kmeans_all14_normal/
}

gen_gff(){
    for FILE in $(ls output/kmeans_all14_normal/name_*)
    do
        python bin/get_gff.py -g input/Staphylococcus_aureus_HG003.gff -n $FILE > ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff
        cp ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff ${FILE}.gff
    done
}

get_protein_names(){
    for FILE in $(ls output/kmeans_all14_normal/name_*.gff)
    do
        NAME_FILE=$(echo $FILE | sed s"/.gff//")
        NAME=$(echo $FILE | sed "s/name/protein_name/")
        PRO=$(echo $NAME | sed "s/.gff/.txt/")
        python bin/get_protein_name.py -i $FILE -s input/Staphylococcus_aureus_HG003_sRNA.csv -n $NAME_FILE> $PRO
    done
}

GO_term(){
    for FILE in $(ls output/kmeans_all14_normal/name_*)
    do
        python bin/get_gff.py -g input/Staphylococcus_aureus_HG003.gff -n $FILE > ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff
        cp ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff ${FILE}.gff
        sh anno.sh
        INDEX=$(echo $FILE | cut -d'_' -f 5)
#        mkdir output/kmeans_all14_normal/fig_kmeans_$INDEX
        mv ANNOgesic/output/Go_term/all_CDS/statistics/Staphylococcus_aureus_HG003/figs/* output/kmeans_all14_normal/fig_kmeans_$INDEX
    done
}

pre_rank_srna(){
    python bin/sort_srna_rpkm.py \
	-i input/gene_wise_quantifications_combined_rpkm.csv \
	-s input/Staphylococcus_aureus_HG003_sRNA.csv \
	> input/rpkm_name.csv
}
# use excel to get the name of sort
rank_srna(){
    python bin/count_sort_rpkm.py -i input/sort_rpkm_name.csv > output/rank_srna_rpkm
    python3 bin/plot_srna.py -r input/rpkm_name.csv -s output/rank_srna_rpkm -l -o output/high_express_srna_log.png
    python3 bin/plot_srna.py -r input/rpkm_name.csv -s output/rank_srna_rpkm -o output/high_express_srna.png
}

anti_kmeans(){
#    python bin/all_rpkm.py -i input/gene_wise_quantifications_combined_rpkm.csv -s input/Staphylococcus_aureus_HG003_sRNA.csv -a > gene_rpkm
    python bin/classify.py -r gene_rpkm -k kmeans.txt -n name_rpkm
    mv kmeans_* output/kmeans_all14_normal_repress
    mv name_* output/kmeans_all14_normal_repress
    rm kmeans.txt
}

anti_plot_kmeans(){
    rm output/kmeans_all14_normal_repress/*.png
    for FILE in $(ls -d output/kmeans_all14_normal_repress/kmeans_*)
    do
        NAME=$(echo $FILE | sed "s/ss\/kmeans/ss\/name_kmeans/")
        python bin/plot_all_kmeans.py -i $FILE -n $NAME -f all -s input/Staphylococcus_aureus_HG003_sRNA.csv -a
    done
    mv *.png output/kmeans_all14_normal_repress/
}

anti_GO_term(){
    for FILE in $(ls output/kmeans_all14_normal_repress/name_*)
    do
        python bin/get_gff.py -g input/Staphylococcus_aureus_HG003.gff -n $FILE > ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff
        cp ANNOgesic/output/target/annotation/Staphylococcus_aureus_HG003.gff ${FILE}.gff
#        sh anno.sh
#        INDEX=$(echo $FILE | cut -d'_' -f 5)
#        mkdir output/kmeans_all14_normal/fig_kmeans_$INDEX
#        mv ANNOgesic/output/Go_term/all_CDS/statistics/Staphylococcus_aureus_HG003/figs/* output/kmeans_all14_normal/fig_kmeans_$INDEX
    done
}

anti_get_protein_names(){
    for FILE in $(ls output/kmeans_all14_normal_repress/name_*.gff)
    do
        NAME=$(echo $FILE | sed "s/name/protein_name/")
        PRO=$(echo $NAME | sed "s/.gff/.txt/")
        python bin/get_protein_name.py -i $FILE > $PRO
    done
}
main
