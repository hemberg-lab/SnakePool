configfile : "config.yaml"

############ Params #############################



##################################################




import glob, os
import random
import csv
from collections import defaultdict

#DATA = [os.path.basename(f).split(".")[0] for f in glob.glob('data/input/*')]
DATA = [os.path.basename(f).split(".")[0] for f in config["samples"]]

ME_clusters = ['NM1_7', 'NM1_5', 'NM1_4', 'NM1_14', 'NM1_13', 'NM1_6', 'NM1_10', 'NM1_9', 'NM1_3', 'NM1_8', 'NM1_11', 'NM1_2', 'N1_5', 'N1_14', 'N1_6', 'N1_10', 'N1_9', 'N1_3', 'N1_8', 'N1_11', 'N1_2', 'NM2_4', 'NM2_14', 'NM2_6', 'NM2_10', 'NM2_9', 'NM2_3', 'NM2_8', 'NM2_11', 'NM2_2', 'N2_4', 'N2_13', 'N2_6', 'N2_10', 'N2_9', 'N2_3', 'N2_8', 'N2_11', 'N2_2', 'NM3_4', 'NM3_13', 'NM3_6', 'NM3_10', 'NM3_9', 'NM3_3', 'NM3_8', 'NM3_11', 'NM3_2', 'N3_6', 'N3_10', 'N3_9', 'N3_3', 'N3_8', 'N3_11', 'N3_2', 'N4_6', 'N4_10', 'N4_9', 'N4_3', 'N4_8', 'N4_11', 'N4_2', 'N5_4', 'N5_13', 'N5_6', 'N5_10', 'N5_9', 'N5_3', 'N5_8', 'N5_11', 'N5_2', 'NN1_13', 'NN1_10', 'NN1_9', 'NN1_3', 'NN1_8', 'NN1_11', 'NN1_2']



####TEMPAll ####

present_quant_files = set([])
quant_files = [os.path.basename(f).split(".")[0] for f in glob.glob('Whippet/Quant/*.psi.gz')]

for f in quant_files:
    present_quant_files.add(f)


# rule all:
#     input:
#         expand("Whippet/BAM/{sample}.sorted.bam.bai", sample=quant_files)

#
#
# rule all:
#    input:
#     #diff_files
#     #diff_files + expand("BAM/{sample}.sorted.bam", sample=REPLICATES)
#     #expand("Whippet/BAM/{sample}.sorted.bam", sample= [x.split("/")[-1].split(".")[0] for x in config["samples"]] )
#     target_pools
#     #expand("Whippet/Delta/{ME_cluster}.ME.diff.gz", ME_cluster=ME_clusters)

#print("asda")

def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

present_quant_files = set([])
quant_files = [os.path.basename(f).split(".")[0] for f in glob.glob('Whippet/Quant/*.psi.gz')]

for f in quant_files:
    present_quant_files.add(f)

#cluster_compare = [ ('GABA-ergic Neuron', 'Glutamatergic Neuron') ]

repeats = 10
#np = 50  #number of pools





#cluster_compare = { "Neuronal-vs-Non_Neuronal" :( ["GABA-ergic Neuron",  "Glutamatergic Neuron"], ["Endothelial Cell", "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte Precursor Cell" ] ),
#                    "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( ["GABA-ergic Neuron"], ['Glutamatergic Neuron'] )     }


#cluster_compare_np = { "Neuronal-vs-Non_Neuronal" :( 100, 10 ),
#                    "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( 50, 50 )     }


cluster_compare = dict()
cluster_compare_np = dict()


#Compare_ID          A.cluster_names A.number_of_pools B.cluster_names B.number_of_pools

with open(config["run_metadata"]) as run:

    run_metadata = csv.DictReader(run, delimiter="\t")

    for row in run_metadata:

        A_cluster_names = []
        B_cluster_names = []

        for c in row["A.cluster_names"].split(","):

            if c in quant_files:

                A_cluster_names.append(c)

        for c in row["B.cluster_names"].split(","):

            B_cluster_names.append(c)

        cluster_compare[row["Compare_ID"]] = (A_cluster_names, B_cluster_names)
        cluster_compare_np[row["Compare_ID"]] = (int(row["A.number_of_pools"]), int(row["B.number_of_pools"]))


print(expand("Whippet/BAM/Merge/{compare_name}.A.sort.bam", compare_name=cluster_compare.keys()))

rule all:
    input:
        expand("Whippet/BAM/Merge/{compare_name}.{x}.sort.bam", compare_name=cluster_compare.keys(), x=["A", "B"])



# cluster_compare = { "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( ["GABA-ergic Neuron"], ['Glutamatergic Neuron'] )}             }


#cluster_compare = {"GABA-ergic_Neuron_vs_Glutamatergic_Neuron" : ( ('GABA-ergic Neuron', "Endothelial Cell"), ('Glutamatergic Neuron', "Endothelial Cell") )}

# cluster_compare = [ ( ('GABA-ergic Neuron', 'Glutamatergic Neuron') ]
#
# [( ('GABA-ergic Neuron'), ('Glutamatergic Neuron') ),
# (("GABA-ergic Neuron",  "Glutamatergic Neuron"), ("Endothelial Cell", "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte Precursor Cell" ) ) ]

cluster_files = defaultdict(list)

#with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/Tasic_clustering.txt") as Tasic:
with open(config["cluster_medatada"]) as Tasic:

    Tasic_clustering = csv.DictReader(Tasic, delimiter="\t")

    for row in Tasic_clustering:

        if row["Run_s"] in present_quant_files:

            cluster_files[row[config["cluster_name"]]].append(row[config["file_basename"]])
            #cluster_files[row["broad_type"]].append(row["Run_s"])


target_pool_psi = []
target_pool_delta = []



for compare_name, c in cluster_compare.items():

    g1, g2 = c

    c1_names = []
    for c1 in g1:

        c1_names += cluster_files[c1]

    c2_names = []
    for c2 in g2:
        c2_names += cluster_files[c2]


    #c1, c2 = c

    #c1_names = cluster_files[c1]
    #c2_names = cluster_files[c2]


    #compare_name = "_vs_".join([x.replace(" ", "_") for x in c ])



    for r in range(repeats):


        delta_name = "Whippet/Delta/Tasic/" + compare_name +  "_pool_" +  str(r+1)

        #print(delta_name)

        target_pool_delta.append( delta_name + ".diff.gz")

# print(target_pool_delta)

# rule all:
#    input:
#     target_pool_delta
#
#
#










rule Splice_Junction_Library:
    input:
        config["Genome_fasta"],
        config["Gene_anontation_fasta"],
        config["Gene_anontation_bed12"]
    params:
        ME_len = config["ME_len"]
    output:
        "Round1/ME_TAGs.fa"
    shell:
        "python2 src/SJ_tags_generator_for_micro_exons.py {input} {params.ME_len} > {output}"


# rule sra_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.sra"
#     output:
#         temp("data/fastq_paired/{sample}.fastq")
#     shell:
#         "fastq-dump {input} -O data/fastq_paired/"


# rule fastq_gz_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.fastq.gz"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "gzip -dc {input} > {output}"
#
# rule fastq_input:
#     input:
#         config["input_dir"] + "/{sample}.fastq"
#     output:
#         "data/fastq/{sample}.fastq"
#     shell:
#         "ln -s {input} {output}"

#rule download_to_fastq:
#    input:
#        "download/{sample}.download.sh"
#    output:
#        "data/fastq/{sample}.fastq"
#    shell:
#        "bash {input}"


# rule split_fastq:
#     input:
#         "data/fastq_paired/{sample}.fastq"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "python2 src/split_paired_end.py {input} > {output}"


rule bwa_index:
    input:
        "Round1/ME_TAGs.fa"
    output:
        "Round1/ME_TAGs.fa.amb"
    shell:
        "bwa index {input}"

rule Round1_bwa_mem_to_tags:
    input:
        "Round1/ME_TAGs.fa",
        "data/fastq/{sample}.fastq",
        "Round1/ME_TAGs.fa.amb"
    output:
        temp("Round1/{sample}.sam")
    threads: 5
    shell:
        "bwa mem -t {threads} -O 2,6 -L 25 {input[0]} {input[1]} > {output}"


rule Round1_alingment_pre_processing:
    input:
        "Round1/{sample}.sam"
    output:
        "Round1/{sample}.sam.pre_processed"
    shell:
        "python2 src/alingment_pre_processing.py {input} F > {output}"



#############################################################################

                     #Round1 - POST-Processing###

#############################################################################




rule row_Micro_Exon_reads:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.pre_processed"
    output:
        "Round1/{sample}.sam.row_ME",
        "Round1/{sample}.sam.row_ME.fastq"
    shell:
        "python2 src/row_ME.py {input} > {output[0]}"


rule hisat2_Genome_index:
    input:
        config["Genome_fasta"]
    output:
        "data/Genome.1.ht2"
    threads: 5
    shell:
        "hisat2-build {input} data/Genome"


rule hisat2_to_Genome:
    input:
        "Round1/{sample}.sam.row_ME.fastq",
        "data/Genome.1.ht2"
    output:
        "Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam"
    threads: 1
    shell:
        "hisat2 -x data/Genome -U {input[0]} > {output}"


rule Round1_filter:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.row_ME",
        "Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam",
        config["GT_AG_U2_5"],
        config["GT_AG_U2_3"],
        config["vertebrates_phylop"]
    params:
        ME_len = config["ME_len"]
    output:
        "Round1/{sample}.sam.row_ME.filter1"
    shell:
        "python2 src/ME_filter1.py {input} {params.ME_len} > {output}"


rule Micro_Exon_table:
    input:
        expand("Round1/{sample}.sam.row_ME.filter1", sample=DATA )
    output:
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric"
    shell:
        "cat {input} > Round1/TOTAL/TOTAL.sam.row_ME.filter1 && python2 src/ME_centric_table.py Round1/TOTAL/TOTAL.sam.row_ME.filter1 > {output}"


rule Micro_Exon_Tags:
    input:
        "Round1/ME_TAGs.fa",
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric"
    output:
        "Round2/ME_canonical_SJ_tags.de_novo.fa"
    shell:
        "python2 src/Micro_exons_tags.py  {input} > {output}"







#############################################################################

                     #Round2#

#############################################################################


rule Get_ME_from_annotation:
    input:
        config["Genome_fasta"],
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        config["Gene_anontation_bed12"],
        config["GT_AG_U2_5"],
        config["GT_AG_U2_3"],
        config["vertebrates_phylop"],
        config["ME_DB"]
    params:
        ME_len = config["ME_len"]
    output:
        "ME_canonical_SJ_tags.DB.fa",
        "DB.ME_centric"
    shell:
        "python2 src/Get_annotated_microexons.py  {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {params.ME_len} {input[6]} "


rule merge_tags:
    input:
        "Round2/ME_canonical_SJ_tags.de_novo.fa",
        "ME_canonical_SJ_tags.DB.fa"
    output:
        "Round2/ME_canonical_SJ_tags.fa"
    shell:
        "cat {input[0]} {input[1]} > {output}"


rule merge_ME_centric:
    input:
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        "DB.ME_centric"
    output:
        "Round2/TOTAL.ME_centric.txt"
    shell:
        "cat {input[0]} {input[1]} > {output}"


rule Round2_bowtie_tags_index:
    input:
        "Round2/ME_canonical_SJ_tags.fa"
    output:
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    shell:
        "bowtie-build {input} {input}"


rule Round2_bowtie_to_tags:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        "data/fastq/{sample}.fastq",
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    output:
        temp("Round2/{sample}.sam")
    threads: 5
    shell:
        "bowtie {input[0]} -p {threads} -q {input[1]} -S -v 2 > {output}"


rule Round2_alingment_pre_processing:
    input:
        "Round2/{sample}.sam"
    output:
        "Round2/{sample}.sam.pre_processed"
    shell:
        "python2 src/alingment_pre_processing_round2_bowtie.py {input} F > {output}"



#############################################################################

                     #Round2 - POST-Processing###

#############################################################################


rule ME_reads:
    input:
        "Round2/{sample}.sam.pre_processed"
    output:
        "Round2/{sample}.sam.pre_processed.fastq"
    shell:
        "python2 src/round2_ME_reads_fastq.py {input} > {output}"


rule bowtie_Genome_index:
    input:
        config["Genome_fasta"]
    output:
        config["Genome_fasta"] + ".1.ebwt"
    shell:
        "bowtie-build {input} {input}"

rule bowtie_to_genome:
    input:
        "Round2/{sample}.sam.pre_processed.fastq",
        config["Genome_fasta"],
        config["Genome_fasta"] + ".1.ebwt"
    output:
        "Round2/{sample}.sam.pre_processed.hg19.sam"
    shell:
        "bowtie {input[1]} -p 1 -q {input[0]} -S -v 2| awk '$2==0 || $2==16'> {output}"


rule Round2_filter:
    input:
        "Round2/{sample}.sam.pre_processed",
        "Round2/{sample}.sam.pre_processed.hg19.sam",
    output:
        "Round2/{sample}.sam.pre_processed.filter1"
    shell:
        "python2 src/Filter1_round2.py {input} > {output}"


rule ME_SJ_coverage:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        "Round2/TOTAL.ME_centric.txt",
        config["Gene_anontation_bed12"],
        "Round2/{sample}.sam.pre_processed.filter1"
    params:
        ME_len = config["ME_len"]
    output:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
    shell:
        "python2 src/ME_SJ_coverage.py {input} {params.ME_len} > {output}"


rule GetPSI:
    input:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
    output:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage.PSI"
    shell:
      "python3 src/GetPSI.py {input} 5 > {output}"


rule Total_sample_exon_counts:
    input:
        expand("Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage", sample=DATA )
    output:
        "Round2/TOTAL.filter1.ME_SJ_coverage"
    shell:
      "cat Round2/*.filter1.ME_SJ_coverage > {output}"



rule Output:
    input:
        "Round2/TOTAL.ME_centric.txt",
        "Round2/TOTAL.filter1.ME_SJ_coverage"
    params:
        wd = config["working_directory"]
    output:
        "Report/out_filtered_ME.txt",
        "Report/out_low_scored_ME.txt",
        "Report/out_shorter_than_3_ME.txt",
        "Report/report.html",
        "Report/out_filtered_ME.cov.txt"
    shell:
        '''R -e  'rmarkdown::render("src/final_filters2.Rmd",params = list(ME_table="{params.wd}{input[0]}", ME_coverage="{params.wd}{input[1]}", out_filtered_ME="{params.wd}{output[0]}", out_low_scored_ME="{params.wd}{output[1]}", out_shorter_than_3_ME="{params.wd}{output[2]}", min_number_files_detected=1, out_filtered_ME_cov="{params.wd}{output[4]}" ), output_file="{params.wd}{output[3]}")'  '''


############## Gene Count #######


#
# rule total_hisat2_to_Genome:
#      input:
#          "data/fastq/{sample}.fastq",
#          "data/Genome.1.ht2"
#      output:
#          "Genome_aligments/{sample}.sam"
#      threads: 5
#      shell:
#          "hisat2 -x data/Genome -U {input[0]} -p 5 > {output}"
#
#
# rule gene_count:
#     input:
#         "/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.gtf",
#         "Genome_aligments/{sample}.sam"
#     output:
#         "Genome_aligments/{sample}.gene_count.txt"
#     threads: 1
#     shell:
#         "featureCounts -a {input[0]} -o {output} {input[1]}"
#
#
#
# rule done_gene_count:
#     input:
#         expand("Genome_aligments/{sample}.gene_count.txt", sample=DATA )
#     output:
#         "Round2/done.txt"
#     shell:
#         "echo done > {output}"
######




###### Whippet ########

rule whippet_index:
    input:
        config["Genome_fasta"],
        "/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/out_filtered_ME.no_chrM.gtf"
    params:
        bin = config["whippet_bin_folder"]
    output:
        "Whippet/Index/whippet.jls",
        "Whippet/Index/whippet.jls.exons.tab.gz"
    shell:
        "julia {params.bin}/whippet-index.jl --fasta {input[0]} --gtf {input[1]} --index {output[0]}"


rule FASTQ_compress:
    input:
        "data/fastq/{sample}.fastq"
    output:
        "data/fastq/{sample}.fastq.gz"
    shell:
        "gzip {input}"


rule  whippet_quant:
    input:
        "data/fastq/{sample}.fastq.gz",
        "Whippet/Index/whippet.jls"
    params:
        bin = config["whippet_bin_folder"],
        output = "Whippet/Quant/{sample}"
    output:
        temp("Whippet/SAM/{sample}.sam"),
        "Whippet/Quant/{sample}.gene.tpm.gz",
        "Whippet/Quant/{sample}.isoform.tpm.gz",
        "Whippet/Quant/{sample}.jnc.gz",
        "Whippet/Quant/{sample}.map.gz",
        "Whippet/Quant/{sample}.psi.gz"
    shell:
        "julia {params.bin}/whippet-quant.jl {input[0]}  -x {input[1]} -o {params.output} --sam > {output[0]}"

rule unzip_psi:
    input:
        "Whippet/Quant/{sample}.psi.gz"
    output:
        "Whippet/Quant/{sample}.psi"
    shell:
        "zcat {input} > {output}"

rule fix_quant:
    input:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage.PSI",
        "Whippet/Quant/{sample}.psi"
    output:
        "Whippet/Quant/{sample}.psi.ME"
    shell:
        "python3 src/Replace_PSI_whippet.py {input} > {output}"


rule samTobam:
    input:
        "Whippet/SAM/{sample}.sam"
    output:
        "Whippet/BAM/{sample}.sorted.bam"
    params:
        samtools = config["samtools"]
    shell:
        "{params.samtools} view -b  {input}  | {params.samtools} sort - -o {output}"


rule indexBam:
    input:
        "Whippet/BAM/{sample}.sorted.bam"
    output:
        "Whippet/BAM/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"





##### Whippet Diff ####


import csv
from collections import defaultdict

Tissue_clusters = defaultdict(list)

with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/Tissue.clusters.txt") as csvfile:

    reader = csv.DictReader(csvfile, delimiter="\t")

    for row in reader:

        Tissue_clusters[int(row['Tissue_clusters'])].append(row["FILE_NAME"])




COMPARE = {"NM1":[ (12, 1) , (7, 5, 4, 14, 13, 6, 10, 9, 3, 8, 11, 2)],   #NM1
        "N1":[ (12, 7, 1, 4, 13) , (5, 14, 6, 10, 9, 3, 8, 11, 2 )],    #N1
        "NM2":[ (12, 7, 1, 5, 13) , (4, 14, 6, 10, 9, 3, 8, 11, 2 )],    #NM2
        "N2":[ (12, 7, 1, 5, 14) , (4, 13, 6, 10, 9, 3, 8, 11, 2 )],    #N2 - N5
        "NM3":[ (12, 7, 1, 5, 14) , (4, 13, 6, 10, 9, 3, 8, 11, 2 )],    #NM3
        "N3":[ (12, 7, 1, 5, 4, 14, 13) , (6, 10, 9, 3, 8, 11, 2 )],    #N3 - N4
        "N4":[ (12, 7, 1, 5, 4, 14, 13) , (6, 10, 9, 3, 8, 11, 2 )],    #N3 - N4
        "N5":[ (12, 7, 1, 5, 14) , (4, 13, 6, 10, 9, 3, 8, 11, 2 )],    #N2 - N5
        "NN1":[ (12, 7, 1, 5, 4, 14, 6) , (13, 10, 9, 3, 8, 11, 2 )]}    #NN1


COMPARE_FILES = dict()
COMPARE_rev_c1 = dict()
COMPARE_rev_c2 = dict()



for ME_cluster, compare in COMPARE.items():

    c1, c2 = compare
    c1_FILES = []
    c2_FILES = []

    for c in c1:
        for f in Tissue_clusters[c]:
            c1_FILES.append(f)



    # for c in c2:
    #     for f in Tissue_clusters[c]:
    #         c2_FILES.append(f)

    #COMPARE_FILES[ME_cluster] = [c1_FILES, c2_FILES]

    COMPARE_FILES[ME_cluster] = [c1_FILES, c2]


#print(COMPARE_FILES)

#print(COMPARE.keys())

ME_clusters_com = []

for ME_cluster, file_conditions in COMPARE_FILES.items():

    c1, c2 = file_conditions

    PSI_c1 = expand("Whippet/Quant/{FILE}.psi.gz", FILE=c1 )
    #PSI_c2 = expand("Whippet/Quant/{FILE}.psi.gz", FILE=c2 )

    for c in c2:

        PSI_c2 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  Tissue_clusters[c] ]

        #print(PSI_c2)

        ME_clusters_com.append(ME_cluster + "_" + str(c))

        rule:
            input:
                PSI_c1 + PSI_c2
            output:
                "Whippet/Delta/" + ME_cluster + "_" + str(c) + ".diff.gz"
            params:
                bin = config["whippet_bin_folder"],
                a = ",".join(PSI_c1),
                b = ",".join(PSI_c2),
                o = "Whippet/Delta/" + ME_cluster + "_" + str(c)
            shell:
                "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"


    # Rules for fixed values


    PSI_c1 = expand("Whippet/Quant/{FILE}.psi.ME", FILE=c1 )
    #PSI_c2 = expand("Whippet/Quant/{FILE}.psi.gz", FILE=c2 )

    for c in c2:

        PSI_c2 = [ "Whippet/Quant/" + x + ".psi.ME" for x in  Tissue_clusters[c] ]

        #print(PSI_c2)

        ME_clusters_com.append(ME_cluster + "_" + str(c))

        rule:
            input:
                PSI_c1 + PSI_c2
            output:
                "Whippet/Delta/" + ME_cluster + "_" + str(c) + ".ME" + ".diff.gz"
            params:
                bin = config["whippet_bin_folder"],
                a = ",".join(PSI_c1),
                b = ",".join(PSI_c2),
                o = "Whippet/Delta/" + ME_cluster + "_" + str(c) + ".ME"
            shell:
                "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"


#print(ME_clusters_com)



##### Single cell ####

import random
import csv
from collections import defaultdict

def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

present_quant_files = set([])
quant_files = [os.path.basename(f).split(".")[0] for f in glob.glob('Whippet/Quant/*.psi.gz')]

for f in quant_files:
    present_quant_files.add(f)


# rule all:
#     input:
#         expand("BAM/{sample}.sorted.bam.bai", sample=quant_files)



#cluster_compare = [ ('GABA-ergic Neuron', 'Glutamatergic Neuron') ]


cluster_files = defaultdict(list)


with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/Tasic_clustering.txt") as Tasic:

    Tasic_clustering = csv.DictReader(Tasic, delimiter="\t")

    for row in Tasic_clustering:

        if row["Run_s"] in present_quant_files:

            cluster_files[row["broad_type"]].append(row["Run_s"])





# 'GABA-ergic Neuron'
# 'Glutamatergic Neuron'

for compare_name, c in cluster_compare.items():


    g1, g2 = c

    c1_names = []
    for c1 in g1:

        c1_names += cluster_files[c1]

    c2_names = []
    for c2 in g2:
        c2_names += cluster_files[c2]


    np_A, np_B = cluster_compare_np[compare_name]



    rule :
        input:
            expand('Whippet/BAM/{sample}.sorted.bam', sample=c1_names)
        output:
            temp("Whippet/BAM/Merge/" + compare_name + ".A.bam")
        shell:
            "samtools merge  {output} {input}"

    rule :
        input:
            "Whippet/BAM/Merge/" + compare_name + ".A.bam"
        output:
            "Whippet/BAM/Merge/" + compare_name + ".A.sort.bam"
        shell:
            'samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}'


    rule :
        input:
            expand('Whippet/BAM/{sample}.sorted.bam', sample=c2_names)
        output:
            temp("Whippet/BAM/Merge/" + compare_name + ".B.bam")
        shell:
            "samtools merge  {output} {input}"

    rule :
        input:
            "Whippet/BAM/Merge/" + compare_name + ".B.bam"
        output:
            "Whippet/BAM/Merge/" + compare_name + ".B.sort.bam"
        shell:
            'samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}'




    for r in range(repeats):


        c1_pools = partition(c1_names, np_A)
        c2_pools = partition(c2_names, np_B)

        p = 0

        target_pool_psi_A = []
        target_pool_psi_B = []

        delta_name = "Whippet/Delta/Tasic/" + compare_name +  "_pool_" +  str(r+1)

        #for pc1, pc2 in zip(c1_pools, c2_pools):

        for pc1 in c1_pools:

            p += 1

            FASTQ_c1 = [ "data/fastq/" + x + ".fastq.gz" for x in  pc1 ]
            PSI_c1 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc1 ]

            pool_ID = "pool_" +str(r + 1) + "_"  + str(p)


            target_pool_psi_A.append("Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".psi.gz")




            rule:
                input:
                    fastq = FASTQ_c1,
                    index = "Whippet/Index/whippet.jls"
                output:
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".gene.tpm.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".isoform.tpm.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".jnc.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".map.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".psi.gz"
                params:
                    bin = config["whippet_bin_folder"],
                    output = "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID
                shell:
                    "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) -x {input.index} --force-gz -o {params.output}"


        for pc2 in c2_pools:

            p += 1

            FASTQ_c2 = [ "data/fastq/" + x + ".fastq.gz" for x in  pc2 ]
            PSI_c2 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc2 ]

            pool_ID = "pool_" +str(r + 1) + "_"  + str(p)


            target_pool_psi_B.append("Whippet/Quant/Tasic/" + compare_name + "_B_" + pool_ID + ".psi.gz")



            rule:
                input:
                    fastq = FASTQ_c2,
                    index = "Whippet/Index/whippet.jls"
                output:
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".gene.tpm.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".isoform.tpm.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".jnc.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".map.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".psi.gz"
                params:
                    bin = config["whippet_bin_folder"],
                    output = "Whippet/Quant/Tasic/" + compare_name + "_B_" + pool_ID
                shell:
                    "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) -x {input.index} --force-gz -o {params.output}"



        rule:
            input:
                target_pool_psi_A + target_pool_psi_B
            output:
                delta_name + ".diff.gz"

            params:
                bin = config["whippet_bin_folder"],
                a = ",".join(target_pool_psi_A),
                b = ",".join(target_pool_psi_B),
                o = delta_name
            shell:
                "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"




    # rule compare:
    #     input:
    #         PSI_c1 + PSI_c2
    #     output:
    #         "Whippet/Delta/Single_cell/" + compare_name + ".diff.gz"
    #     params:
    #         bin = config["whippet_bin_folder"],
    #         a = ",".join(PSI_c1),
    #         b = ",".join(PSI_c2),
    #         o = "Whippet/Delta/Single_cell/" + compare_name
    #     shell:
    #         "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"


#"Whippet/Delta/GABA-ergic_Neuron_vs_Glutamatergic_Neuron.diff.gz


#defaultdict(<class 'int'>, {'GABA-ergic Neuron': 739, 'Glutamatergic Neuron': 761, 'Endothelial Cell': 29, 'Astrocyte': 43, 'Microglia': 22, 'Oligodendrocyte': 38, 'Oligodendrocyte Precursor Cell': 22})

#
# rule all:
#    input:
#     expand("Whippet/Delta/{ME_cluster}.diff.gz", ME_cluster=ME_clusters_com)



#snakemake -j 30000  --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -k


#snakemake -j 30000 Report/out_filtered_ME.txt --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -k

#        "R -e" + ' rmarkdown::render("src/final_filters.Rmd", params = list(ME_table={input[0]}, ME_coverage={input[1]}, out_filtered_ME={output[0]}, out_low_scored_ME={output[1]}, out_shorter_than_3_ME={output[2]}), output_file={output[3]})'



# snakemake --dag Round2/TOTAL.filter1.ME_SJ_coverage | dot -Tsvg > Micro-Exonator.svg

# snakemake --dag Report/report.html | dot -Tsvg > Micro-Exonator.svg


#R -e  'rmarkdown::render("src/final_filters.Rmd", params = list(ME_table="/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric", ME_coverage="/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Round2/TOTAL.filter1.ME_SJ_coverage", out_filtered_ME= "/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_filtered_ME.txt", out_low_scored_ME="/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_low_scored_ME.txt", out_shorter_than_3_ME="/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_shorter_than_3_ME.txt", min_number_files_detected=1), output_file="/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/report.html")'


#R -e  'rmarkdown::render("src/final_filters.Rmd",params = list(ME_table="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric", ME_coverage="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Round2/TOTAL.filter1.ME_SJ_coverage", out_filtered_ME="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_filtered_ME.txt", out_low_scored_ME="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_low_scored_ME.txt", out_shorter_than_3_ME="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/out_shorter_than_3_ME.txt", min_number_files_detected=1), output_file="lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator/Report/report.html")'
