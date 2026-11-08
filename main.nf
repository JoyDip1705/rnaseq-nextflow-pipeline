#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.read = "/mnt/d/Nextflow_pipeline/raw_data/SRR*_{1,2}.fastq.gz"
params.ref = "/mnt/d/Nextflow_pipeline/ref_genome/hg38.fna"
params.adapters = "/mnt/d/Nextflow_pipeline/adapters/TruSeq3-PE.fa"
params.output = "/mnt/d/Nextflow_pipeline/result"
params.annot = "/mnt/d/Nextflow_pipeline/ref_genome/hg38.gtf"
params.meta = "/mnt/d/Nextflow_pipeline/sample_metadata.csv"


process QUALITY {
   tag "Quality of the sample: ${sample_id}"
   publishDir "${params.output}/qc_report", mode: 'copy'

   input:
   tuple val(sample_id), path(data)

   output:
   path "*.{html,zip}", emit: html_zip

   script:
   """
   fastqc -t ${task.cpus} ${data}
   """
}

                                                     
process TRIMMING {
   tag "Trimming of the sample: ${sample_id}"
   publishDir "${params.output}/trimming_output", mode: 'copy'

   input:
   tuple val(sample_id), path(data)
   path(adapter)

   output:
   tuple val(sample_id), path("${sample_id}_R1_paired.fastq.gz"), path("${sample_id}_R2_paired.fastq.gz"), emit: trim_files

   script:
   """
   trimmomatic PE -threads ${task.cpus} -phred33 ${data[0]} ${data[1]} \
   ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
   ${sample_id}_R2_paired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \
   ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
   """
}

process INDEXING {
   tag "Indexing of the Genome: Human"
   publishDir "${params.output}/indexing_output", mode: 'copy'

   input:
   path(genome_file)

   output:
   path "genome_index.*", emit: index_file

   script:
   """
   hisat2-build -p ${task.cpus} ${genome_file} genome_index
   """

}

process MAPPING_SORT {
   tag "Mapping of the sample: ${sample_id}"
   publishDir "${params.output}/mapping_output", mode: 'copy'

   input:
   tuple val(sample_id), path(data_R1), path(data_R2)
   path(index_file)

   output:
   tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam_bai

   script:
   """
   hisat2 -x genome_index -p ${task.cpus} -1 ${data_R1} -2 ${data_R2} | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam -

   samtools index ${sample_id}_sorted.bam
   """

}

process QUANTIFICATION {
   tag "Quantification of the sample: ${sample_id}"
   publishDir "${params.output}/quantification_output", mode: 'copy'

   input:
   tuple val(sample_id), path(bam), path(bai)
   path(annot_file)
 
   output:
   tuple val(sample_id), path("${sample_id}_counts.txt"), emit: count_file
   path "${sample_id}_counts.txt.summary", emit: summary_file

   script:
   """
   featureCounts -p -B -C  -T ${task.cpus} -a ${annot_file} -o ${sample_id}_counts.txt ${bam}
   """

}

process MERGING_COUNT {
   tag "Merging all count files"
   publishDir "${params.output}/all_merged_counts", mode: 'copy'

   input:
   path(count_files)

   output:
   path "raw_counts_matrix.csv", emit: merge_files

   script:
   """
   merge_counts.py ${count_files}
   """
}

process DEG {
   tag "Identifying DEGs of late vs early luad samples"
   publishDir "${params.output}/DEGs_output", mode: 'copy'

   input:
   path(count)
   path(metadata)

   output:
   path "DESeq2_results.csv", emit: results
   path "significant_DEGs.csv", emit: sig_degs
   path "normalized_counts.csv", emit: norm_counts
   path "volcano_plot.png", emit: volcano
   path "MA_plot.png", emit: ma
   path "PCA_plot.png", emit: pca
   path "heatmap_top50.png", emit: heatmap
   path "DESeq2_object.rds", emit: rds

   script:
   """
   DESeq2.R ${count} ${metadata}
   """

}


workflow {
   ch_rawdata = Channel.fromFilePairs(params.read, checkIfExists: true)
   ch_adapter = Channel.fromPath(params.adapters, checkIfExists: true).first()
   ch_ref = Channel.fromPath(params.ref, checkIfExists: true)
   ch_annot = Channel.fromPath(params.annot, checkIfExists: true).first()
   ch_meta = Channel.fromPath(params.meta, checkIfExists: true).first()

   qc_out = QUALITY(ch_rawdata)

   trim_out = TRIMMING(ch_rawdata, ch_adapter)

   index_out = INDEXING(ch_ref)

   mapping_out = MAPPING_SORT(trim_out.trim_files, index_out.index_file.collect())

   quant_out = QUANTIFICATION(mapping_out.bam_bai, ch_annot)

   ch_merge_counts = quant_out.count_file.map{ sample_id, counts -> counts }.collect()

   merge_all = MERGING_COUNT(ch_merge_counts)

   deg_out = DEG(merge_all.merge_files, ch_meta)


}



