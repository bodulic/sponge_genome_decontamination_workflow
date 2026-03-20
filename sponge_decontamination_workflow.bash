#!/bin/bash
#Master script for the analysis performed in the paper "Hidden Contaminants in Sponge Genomes: Large-Scale Decontamination of 30 Public Assemblies" (Bodulić and Vlahoviček, 2026)

#Downloading analysed sponge genome assemblies from NCBI Genomes (RefSeq + GenBank)
echo "Downloading analysed sponge genome assemblies"
tail -n +2 analysed_sponge_genomes.tsv | while IFS=$'\t' read -r species accession
do
 datasets download genome accession "${accession}" --filename tmp.zip --include genome
 unzip -q tmp.zip -d tmp_dir
 genome_fasta="$(find tmp_dir -name "*.fna" | head -n 1)"
 mv "${genome_fasta}" "${species}.fa"
 rm -rf tmp.zip tmp_dir
done
mkdir analysed_sponge_genomes && mv *.fa analysed_sponge_genomes

#Validation spike-in dataset creation
#Downloading sponge chromosomes for the validation dataset
echo "Downloading sponge chromosomes for the validation dataset"
curl -fL http://hex.bioinfo.hr/~kbodulic/genome_decon_validation_dataset/validation_sponge_chromosomes.tar.gz -o validation_sponge_chromosomes.tar.gz
tar -xzf validation_sponge_chromosomes.tar.gz && rm validation_sponge_chromosomes.tar.gz

#Downloading contaminant scaffolds for the validation dataset
echo "Downloading contaminant scaffolds for the validation dataset"
curl -fL http://hex.bioinfo.hr/~kbodulic/genome_decon_validation_dataset/validation_contaminant_scaffolds.tar.gz -o validation_contaminant_scaffolds.tar.gz
tar -xzf validation_contaminant_scaffolds.tar.gz && rm validation_contaminant_scaffolds.tar.gz

#Generating the validation dataset
echo "Generating the validation dataset"
generate_validation_dataset.R 1000 "1000,2000,4000,10000,20000,30000,50000,100000,200000,500000,1000000,5000000,10000000,100000000" 14500 100 5 1

#Downloading and builiding prerequisite databases for the decontamination pipekine
mkdir databases && cd databases
echo "Downloading NCBI NR database"
curl -fL "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" -o nr.gz

echo "Downloading NCBI taxonomy mapping files"
curl -fL "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz" -o prot.accession2taxid.gz
gunzip prot.accession2taxid.gz

curl -fL "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" -o taxdump.tar.gz
tar -xzf taxdump.tar.gz

echo "Creating DIAMOND NR database"
diamond makedb --in nr.gz --taxonmap prot.accession2taxid --taxonnodes nodes.dmp --taxonnames names.dmp --threads 20 --db nr

echo "Downloading MEGAN taxonomy mapping database"
curl -fL "https://software-ab.cs.uni-tuebingen.de/download/megan7/megan-nr-r2.zip" -o megan-nr-r2.zip
unzip megan-nr-r2.zip

echo "Downloading NCBI core NT database Kraken2 index"
curl -fL "https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz" -o kraken2_core_nt.tar.gz
mkdir kraken2_core_nt
tar -xzf kraken2_core_nt.tar.gz -C kraken2_core_nt

echo "Building Taxonomizr database"
Rscript -e "taxonomizr::prepareDatabase(types='prot')"
cd ..

#Running the decontamination pipeline on the validation dataset
echo "Running the decontamination pipeline on the validation dataset"
decontamination_pipeline -n databases/nr.dmnd -k databases/kraken2_core_nt -T databases/nameNode.sqlite -G 20 -t 20 validation_dataset.fa Metazoa
cd validation_dataset.fa_decontamination_dir && mv scaffold_contamination_info.tsv scaffold_contamination_info_validation.tsv && mv contaminant_scaffolds_taxonomy.tsv contaminant_scaffolds_taxonomy_validation.tsv && ln -s ../taxonomy_info_validation.csv . && cd ..

#Running the decontamination pipeline on the analysed sponge genome assemblies
echo "Running the decontamination pipeline on the anaylsed sponge genome assemblies"
cd analysed_sponge_genomes
for file in $(ls)
do
 decontamination_pipeline -n ../databases/nr.dmnd -k ../databases/kraken2_core_nt -T ../databases/nameNode.sqlite -t 20 "${file}" Metazoa
 cd "${file}_decontamination_dir" && mv scaffold_contamination_info.tsv "${file}_scaffold_contamination_info.tsv" && mv contaminant_scaffolds_taxonomy.tsv "${file}_contaminant_scaffolds_taxonomy.tsv" && mv mc_umap_dt.tsv "${file}_mc_umap_dt.tsv" && cd ..
done

#Merging sponge decontamination results
mkdir sponge_decontamination_results
find . -type f -name "*_scaffold_contamination_info.tsv" | xargs -I{} ln {} sponge_decontamination_results/
find . -type f -name "*_contaminant_scaffolds_taxonomy.tsv" | xargs -I{} ln {} sponge_decontamination_results/
find . -type f -name "*_mc_umap_dt.tsv" | xargs -I{} ln {} sponge_decontamination_results/
cd sponge_decontamination_results

awk -F'\t' 'BEGIN{OFS="\t"} FNR==1{if(NR==1) print $0,"filename"; next} {print $0,FILENAME}' *_scaffold_contamination_info.tsv > scaffold_contamination_info_sponges.tsv
awk -F'\t' 'BEGIN{OFS="\t"} FNR==1{if(NR==1) print $0,"filename"; next} {print $0,FILENAME}' *_contaminant_scaffolds_taxonomy.tsv > contaminant_scaffolds_taxonomy_sponges.tsv
awk -F'\t' 'BEGIN{OFS="\t"} FNR==1{if(NR==1) print $0,"filename"; next} {print $0,FILENAME}' *_mc_umap_dt.tsv > mc_umap_dt_sponges.tsv
cd ..

#Running BUSCO on the original sponge genome assemblies
for file in $(ls | grep .fa$)
do
 busco -i "${file}" -m genome -l metazoa -c 20 -o "${file}_busco_original_genome"
 cd "${file}_busco_original_genome" && grep "Complete BUSCOs (C)" short_summary.specific.metazoa*.txt | awk '{print $1}' > "${file}_busco_original_genome_complete" && cd ..
done

#Running BUSCO on the decontaminated sponge genome assemblies
for dir in $(ls | grep _decontamination_dir)
do
 file="${dir%_decontamination_dir}"
 cd "${dir}"
 busco -i decontaminated_genome.fa -m genome -l metazoa -c 20 -o "${file}_busco_decontaminated_genome"
 cd "${file}_busco_decontaminated_genome" && grep "Complete BUSCOs (C)" short_summary.specific.metazoa*.txt | awk '{print $1}' > "${file}_busco_decontaminated_genome_complete" && cd ..
 cd ..
done

#Merging BUSCO results
mkdir busco_results
find . -type f -name "*busco_original_genome_complete*" | xargs -I{} ln -sr "{}" busco_results/
find . -type f -name "*busco_decontaminated_genome_complete*" | xargs -I{} ln -sr "{}" busco_results/
cd busco_results

for file in *busco_original_genome_complete
do
 printf "%s\t%s\n" "$(cat "${file}")" "${file}"
done > busco_original_genome_complete_merged.tsv

for file in *busco_decontaminated_genome_complete
do
 printf "%s\t%s\n" "$(cat "${file}")" "${file}"
done > busco_decontaminated_genome_complete_merged.tsv
cd ../../

#Analysing validation dataset decontamination results (Figure 2)
cd validation_dataset.fa_decontamination_dir
generate_fig2_elements.R
cd ..

#Analysing decontaminated sponge genome assemblies (Figure 3)
cd analysed_sponge_genomes
mkdir decontamination_analysis && cd decontamination_analysis
ln -s ../sponge_decontamination_results/scaffold_contamination_info_sponges.tsv
ln -s ../sponge_decontamination_results/contaminant_scaffolds_taxonomy_sponges.tsv
ln -s ../sponge_decontamination_results/mc_umap_dt_sponges.tsv
ln -s ../busco_results/busco_original_genome_complete_merged.tsv
ln -s ../busco_results/busco_decontaminated_genome_complete_merged.tsv
generate_fig3_elements.R
cd ../../

exit 0
