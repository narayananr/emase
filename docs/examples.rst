==============
Real Use Cases
==============

Estimating allele-specific expression from human RNA-seq data
-------------------------------------------------------------

**1. Process reference data**

We need to first extract transcript information from the reference genome and gene annotation. Most importantly,
EMASE requires the list of transcript ID's and which gene each transcript belong to.::

    prepare-emase -G ${REF_FASTA} -g ${REF_GTF} -o ${REF_DIR} -m --no-bowtie-index

'prepare-emase' generates the following files for the reference genome::

    ${REF_DIR}/emase.transcriptome.fa
    ${REF_DIR}/emase.transcriptome.info
    ${REF_DIR}/emase.gene2transcripts.tsv

**2. Build an individualized genome**

We assume there is a vcf file that contains phased variant information for every sample of your population. Unless we
know which allele is M(aternal) or P(aternal), we are going to distinguish two alleles with suffices, L(eft) and
R(ight). We also recommend to use different ${SAMPLE_DIR} for each sample::

    python build_new_sequence_from_vcfs.py -r ${REFERENCE_FASTA} \
                                           -i ${INDEL_VCF} \
                                           -s ${SNP_VCF} \
                                           -d ${SAMPLE_DIR} \
                                           -o L.fa \
                                           ${SAMPLE_HAP1_ID_IN_VCF}
    python build_new_sequence_from_vcfs.py -r ${REFERENCE_FASTA} \
                                           -i ${INDEL_VCF} \
                                           -s ${SNP_VCF} \
                                           -d ${SAMPLE_DIR} \
                                           -o R.fa \
                                           ${SAMPLE_HAP2_ID_IN_VCF}

**3. Individualize gene annotation**

We want to incorporate individual variation into the gene annotation too so we can build individualized transcriptome in
the following step::

    python adjust_annotations.py -s 4 -e 5 -c 1 -t 9 \
                                 -d ${SAMPLE_DIR} \
                                 -o L.gtf \
                                 -C ${SAMPLE_HAP1_ID_IN_VCF}_comments.txt
                                 ${REFERENCE_GTF} \
                                 ${SAMPLE_HAP1_ID_IN_VCF}
    python adjust_annotations.py -s 4 -e 5 -c 1 -t 9 \
                                 -d ${SAMPLE_DIR} \
                                 -o R.gtf \
                                 -C ${SAMPLE_HAP1_ID_IN_VCF}_comments.txt
                                 ${REFERENCE_GTF} \
                                 ${SAMPLE_HAP2_ID_IN_VCF}

**4. Create a personalized diploid transcriptome**

From L.fa, R.fa, and the corresponding gtf files, we are going to create diploid transcriptome and other
information that EMASE requires. We assume bowtie v1.0.0 or newer is available.::

    prepare-emase -G ${SAMPLE_DIR}/L.fa,${SAMPLE_DIR}/R.fa -s L,R -o ${SAMPLE_DIR}

This will generate the following files::

    ${SAMPLE_DIR}/emase.pooled.transcriptome.fa
    ${SAMPLE_DIR}/emase.pooled.transcriptome.info
    ${SAMPLE_DIR}/bowtie.transcriptome.1.ebwt
    ${SAMPLE_DIR}/bowtie.transcriptome.2.ebwt
    ${SAMPLE_DIR}/bowtie.transcriptome.3.ebwt
    ${SAMPLE_DIR}/bowtie.transcriptome.4.ebwt
    ${SAMPLE_DIR}/bowtie.transcriptome.rev.1.ebwt
    ${SAMPLE_DIR}/bowtie.transcriptome.rev.2.ebwt

**5. Align RNA-seq reads against the diploid transcriptome**

Although EMASE is a flexible framework for many other alignment strategies, the current version of EMASE was most
intensely tested with bowtie1 transcriptome alignments with the following parameters::

    bowtie -q -a --best --strata --sam -v 3 ${SAMPLE_DIR}/bowtie.transcriptome ${FASTQ_FILE} \
        | samtools view -bS - > ${SAMPLE_DIR}/bowtie.transcriptome.bam

**6. Convert bam file into the emase format**

EMASE runs on an alignment profile of three-dimensional incidence matrix. We convert an alignment file (bam) to
the EMASE format using the following script::

    bam-to-emase -a ${SAMPLE_DIR}/bowtie.transcriptome.bam \
                 -i ${REF_DIR}/emase.transcriptome.info \
                 -s L,R \
                 -o ${SAMPLE_DIR}/bowtie.transcriptome.h5

**7. Run EMASE**

Now we are ready to run EMASE::

    run-emase -i ${SAMPLE_DIR}/bowtie.transcriptome.h5 \
              -g ${REF_DIR}/emase.gene2transcripts.tsv \
              -L ${SAMPLE_DIR}/emase.pooled.transcriptome.info \
              -M ${MODEL} \
              -o ${SAMPLE_DIR}/emase \
              -r ${READ_LENGTH} \
              -c

'run-emase' outputs the following files as a result::

    ${SAMPLE_DIR}/emase.isoforms.effective_read_counts
    ${SAMPLE_DIR}/emase.isoforms.tpm
    ${SAMPLE_DIR}/emase.isoforms.alignment_counts
    ${SAMPLE_DIR}/emase.genes.effective_read_counts
    ${SAMPLE_DIR}/emase.genes.tpm
    ${SAMPLE_DIR}/emase.genes.alignment_counts


Estimating allele-specific binding from ChIP-seq data
-----------------------------------------------------

We assume you have a set of individualized genome and annotation files, in this example, S1 and S2, created by Seqnature
package. We also assume you have a bed file that specifies genomic regions of your interest. First, you need to convert
your bed file into a simple gtf format::

    bed-to-gtf -i targets.bed -o targets.gtf

The targets.gtf files should be modified according to the strains of our interest::

    python adjust_annotations.py -s 4 -e 5 -c 1 -t 9 -o S1.gtf -C S1_comments.txt targets.gtf S1
    python adjust_annotations.py -s 4 -e 5 -c 1 -t 9 -o S2.gtf -C S2_comments.txt targets.gtf S2

Finally, run::

    prepare-emase -G S1.fa,S2.fa -g S1.gtf,S2.gtf -s S1,S2 -o S1xS2

This will store the following files in the folder 'S1xS2'::

    S1xS2/emase.pooled.transcriptome.fa
    S1xS2/emase.pooled.transcriptome.info
    S1xS2/bowtie.transcriptome.1.ebwt
    S1xS2/bowtie.transcriptome.2.ebwt
    S1xS2/bowtie.transcriptome.3.ebwt
    S1xS2/bowtie.transcriptome.4.ebwt
    S1xS2/bowtie.transcriptome.rev.1.ebwt
    S1xS2/bowtie.transcriptome.rev.2.ebwt

Now you can align your RNA-seq reads against the pooled bowtie index of target region::

    bowtie -q -a --best --strata --sam -v 3 S1xS2/bowtie.transcriptome ${FASTQ_FILE} \
        | samtools view -bS - > S1xS2/bowtie.transcriptome.bam

Next, we convert the alignment file into a format that EMASE use for running EM algorithm::

    bam-to-emase -a S1xS2/bowtie.transcriptome.bam \
                 -i S1xS2/emase.transcriptome.info \
                 -s S1,S2 \
                 -o S1xS2/bowtie.transcriptome.h5


It is now ready to run emase. We assume the read length is 100bp::

    run-emase -i bowtie.transcriptome.h5 -L S1xS2/emase.pooled.transcriptome.info -M 4 -c


Estimating allele-specific expression from a F1 sample
----------------------------------------------------------------------------------------
To estimate allele-specific expression from RNA-seq data from a F1 hybrid, we need to have the transcriptomes of the two parental strains of F1. In this example, we will be using F1 hybrids from two inbred mouse strains B6 and CAST.

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
    g2gtools vcf2chain -f ${REF} -i ${VCF_INDELS} -s ${STRAIN} -o ${STRAIN}/REF-to-${STRAIN}.chain
    
    ### patch SNPs 
    g2gtools patch -i ${REF} -s ${STRAIN} -v ${VCF_SNPS} -o ${STRAIN}/${STRAIN}.patched.fa
    
    ### use the chain file and patched fa file to create the new strain's genome fasta file
    g2gtools transform -i ${STRAIN}/${STRAIN}.patched.fa -c ${STRAIN}/REF-to-${STRAIN}.chain -o ${STRAIN}/${STRAIN}.fa

    ### create new strain's annotation file
    g2gtools convert -c ${STRAIN}/REF-to-${STRAIN}.chain -i ${GTF} -f gtf -o ${STRAIN}/${STRAIN}.gtf
   
    ### we can also extract the transcripts from the new strain's genome (This step is not needed for F1 ase analysis)
    g2gtools gtf2db -i ${STRAIN}/${STRAIN}.gtf -o ${STRAIN}/${STRAIN}.gtf.db
    g2gtools extract --transcripts -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.transcripts.fa

Now we have, two parental genomes (B6 and CAST) and their annotation files in GTF format. 

Now, we can use **emase**. Let us first, use **prepare-emase** command in emase, to create diploid (F1) transcritome using the strain-specific transcriptomes::
    
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
    
    OUT_DIR=path_to_output
    
    ### create diploid transcriptome of the F1 using prepare-emase
    ### -m option also outputs gene to transcript information
    prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -s ${SUFFIX1},${SUFFIX2} -o ${OUT_DIR} -m

This will create two files in the ${OUT_DIR}::
  
    emase.pooled.transcripts.fa
    emase.pooled.transcripts.info

emase.pooled.transcripts.fa contains all the transcripts in two parental genomes of F1 hybrid and each allele has distinct identifier.
For example, the two alleles of the transcript "ENSMUST00000000001" will be represented as::

    >ENSMUST00000000001_B
    CACACATCCGGTTCTTCCGGGAGCTAGGGGAGCTGACGGAGAAGGCCACCGCCCAGCAGA
    AGACCCGTCTCCGCCGGTGTGTGGCGATTCCCGCGGTGTGTGTGAGTGAGCCCGGGCCCG
    >ENSMUST00000000001_C
    CACACATCCGGTTCTTCCGGGAGCTAGGGGAGCTGACGGAGAAGGCCACCGCCCAGCAGA
    AGACCCGTCCGCCGGTGTGTGGCGATTCCCGCGGTGTGTGTGAGTGAGCCCGGGCCCG

emase.pooled.transcripts.info contains all the transcripts in two parental genomes of F1 hybrid and their lengths.
For example, the two alleles of the transcript "ENSMUST00000000001" will be represented as::
    
    ENSMUST00000000001_B    120
    ENSMUST00000000001_C    118 




