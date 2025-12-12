wget http://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz

gunzip Danio_rerio.GRCz11.112.gtf.gz


wget http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz

gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz

#Filter the protein coding regions

cellranger mkgtf \
    Danio_rerio.GRCz11.112.gtf \
    Danio_rerio.GRCz11.112.filtered.gtf \
    --attribute=gene_biotype:protein_coding


cellranger mkref \
--genome=danio_rerio_genome_112 \
--fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa \
--genes=Danio_rerio.GRCz11.112.filtered.gtf

#make a count directory to store the output in

mkdir count

#configure the count function to run from the project directory

cellranger count --id=output \
                --description=Alison_Voon_Magor_Danio_Rerio_3prime_gene_expression \
                --fastqs=count/fastqs \
                --sample=SI-TT-F11 \
                --transcriptome=reference/danio_rerio_genome_112 \
                --create-bam=false

