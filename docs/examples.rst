==============
Real Use Cases
==============

Estimating allele-specific expression from a F1 sample
----------------------------------------------------------------------------------------

To estimate allele-specific expression from RNA-seq data from a F1 hybrid, we need to have the transcriptomes of the two parental strains of F1. In this example, we will be using F1 hybrids from two inbred mouse strains B6 and CAST.

Building strain specific genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us use **g2gtools** (https://github.com/churchill-lab/g2gtools) to create parental genomes using known SNPs and Indels and create strain-specific annotation file in GTF format using g2gtools. Since B6 strain is almost the same as the strain of mouse reference genome, we will be creating CAST genome using the Sanger SNP and Indel data (ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/).::

    ### reference genome in fasta format
    REF_GENOME=path_to_mouse_reference_genome(B6)

    ### vcf file for indels
    ### for e.g. :mgp_v5_merged_indels_dbSNP142_CAST.vcf.gz
    INDELS_VCF=path_to_Sanger_vcf_file_for_indels

    ### vcf file for snps (gzipped) for e.g. mgp_v5_merged_snps_dbSNP142_CAST.vcf.gz
    SNPS_VCF=path_to_Sanger_vcf_file_for_snps

    ### gene annotation file in gtf format downloaded from ensembl
    ### (ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus)
    GTF=path_to_ensembl_annotation_in_gtf_format_Mus_musculus.GRCm38.84.gtf 
    
    ### iname of the strain that we want to use SNPs and Indels
    ### (it should match a column name in the Sanger vcf file), e.g., CAST_EiJ
    STRAIN=CAST_EiJ

    ### Create a chain file for mapping bases between two genomes. In this case, between reference and some other strain, like CAST_EiJ:
    g2gtools vcf2chain -f ${REF} -i ${INDELS_VCF} -s ${STRAIN} -o ${STRAIN}/REF-to-${STRAIN}.chain
    
    ### patch SNPs 
    g2gtools patch -i ${REF} -s ${STRAIN} -v ${SNPS_VCF} -o ${STRAIN}/${STRAIN}.patched.fa
    
    ### use the chain file and patched fa file to create the new strain's genome fasta file
    g2gtools transform -i ${STRAIN}/${STRAIN}.patched.fa -c ${STRAIN}/REF-to-${STRAIN}.chain -o ${STRAIN}/${STRAIN}.fa

    ### create new strain's annotation file
    g2gtools convert -c ${STRAIN}/REF-to-${STRAIN}.chain -i ${GTF} -f gtf -o ${STRAIN}/${STRAIN}.gtf
   
    ### We can also extract the transcripts from the new strain's genome (This step is not needed for F1 ase analysis)
    g2gtools gtf2db -i ${STRAIN}/${STRAIN}.gtf -o ${STRAIN}/${STRAIN}.gtf.db
    g2gtools extract --transcripts -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.transcripts.fa

Building Diploid Transcriptome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we have, two parental genomes (B6 and CAST) and their annotation files in GTF format. We can use **emase**. Let us first, use **prepare-emase** command in emase, to create diploid (F1) transcritome using the strain-specific transcriptomes::
    
    ### B6 genome fasta file
    GENOME1=path_to_genome_fasta_file_of_parent1
    
    ### CAST genome fasta file
    GENOME2=path_to_genome_fasta_file_of_parent2
    
    ### B6 annotation file in GTF format
    GTF1=path_to_annotation_file_in_GTF_format_for_parent1
    
    ### CAST annotation file in GTF format
    GTF1=path_to_annotation_file_in_GTF_format_for_parent2
    
    ### identifier of parent1 (allele1 in F1)
    SUFFIX1=B
    
    ### identifier of parent2 (allele2 in F1)
    SUFFIX2=C
    
    EMASE_DIR=path_to_output
    
    ### create diploid transcriptome of the F1 using prepare-emase
    ### -m option also outputs all transcripts for every gene
    prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -s ${SUFFIX1},${SUFFIX2} -o ${EMASE_DIR} -m 

This will create the following files in the ${OUT_DIR}::
  
    ${EMASE_DIR}/emase.pooled.transcripts.fa
    ${EMASE_DIR}/emase.pooled.transcripts.info
    ${EMASE_DIR}/emase.gene2transcripts.tsv
    ${EMASE_DIR}/bowtie.transcripts.1.ebwt
    ${EMASE_DIR}/bowtie.transcripts.2.ebwt
    ${EMASE_DIR}/bowtie.transcripts.3.ebwt
    ${EMASE_DIR}/bowtie.transcripts.4.ebwt
    ${EMASE_DIR}/bowtie.transcripts.rev.1.ebwt
    ${EMASE_DIR}/bowtie.transcripts.rev.2.ebwt

**emase.pooled.transcripts.fa** contains all the transcripts in two parental genomes of F1 hybrid and each allele has distinct identifier.
For example, the two alleles of the transcript "ENSMUST00000000001" will be represented as::

    >ENSMUST00000000001_B
    CACACATCCGGTTCTTCCGGGAGCTAGGGGAGCTGACGGAGAAGGCCACCGCCCAGCAGA
    AGACCCGTCTCCGCCGGTGTGTGGCGATTCCCGCGGTGTGTGTGAGTGAGCCCGGGCCCG
    >ENSMUST00000000001_C
    CACACATCCGGTTCTTCCGGGAGCTAGGGGAGCTGACGGAGAAGGCCACCGCCCAGCAGA
    AGACCCGTCCGCCGGTGTGTGGCGATTCCCGCGGTGTGTGTGAGTGAGCCCGGGCCCG

**emase.pooled.transcripts.info** contains all the transcripts in two parental genomes of F1 hybrid and their lengths.
For example, the two alleles of the transcript "ENSMUST00000000001" will be represented as::
    
    ENSMUST00000000001_B    120
    ENSMUST00000000001_C    118 

Quantfying ASE with Single End reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Aligning RNA-seq reads to the diploid transcriptome using bowtie1::

    bowtie -q -a --best --strata --sam -v 3 ${EMASE_DIR}/bowtie.transcriptome ${FASTQ} \
           | samtools view -bS - > ${BAM_FILE}


# Converting bam file to alignment profile in *emase*'s h5 format::

    bam-to-emase -a ${BAM_FILE} \
             -i ${EMASE_DIR}/emase.pooled.transcripts.info \
             -s B,C \
             -o emase.alignment.profile.h5


# Run EMASE

Now we are ready to run EMASE::

    run-emase -i emase.alignment.profile.h5 \
          -g emase.gene2transcripts.tsv \
          -L ${EMASE_DIR}/emase.pooled.transcripts.info \
          -M ${MODEL} \
          -o ${SAMPLE_DIR}/emase \
          -r ${READ_LENGTH} \
          -c

