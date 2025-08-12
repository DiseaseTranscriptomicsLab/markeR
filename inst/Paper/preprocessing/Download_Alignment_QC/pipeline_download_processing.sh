#!/bin/bash
# to launch pipeline:
#     screen -S benchmarking
#     conda activate benchmarking
#     ./pipeline_benchmarking.sh 2>&1 | tee logPipeline.txt 
# to kill pipeline: 
#     screen -XS <id_screen> quit 

# CONDA ENVIRONMENT ---------------------------------------------------------------------------------- 
# python: v3.11.5
# kallisto: conda install bioconda::kallisto #(v0.44.0)
# fasterq-dump: conda install bioconda::sra-tools #(v2.11.0)
# fastqc: conda install bioconda::fastqc  #(v0.12.1)
# multiqc: conda install multiqc -c conda-forge  #(v1.14)

# ----------------------------------------- DEFINE VARIABLES -----------------------------------------

# User-defined variables (ONLY CHANGE THE SCRIPT HERE!!!!!) ------------------------------------------ 

working_directory="path/to/working/directory"

dataset_table_aux="path/to/text_file/sampleIDs.txt"
dataset_table="path/where/formatted/text_file/sampleIDs/will/be/stored.txt"
dataset_description="NameOfRun"

N=12 # Number of samples to process in parallel
threads=20


# TRANSCRIPTOME  
# FASTA file - RefSeq Release 109 (19/11/2021), downloaded in April 6th, 2022, from:
# https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
# Correspondence isoform <-> gene retreived from UCSC's genome browser for RefSeq release 109
 
transcriptome="path/to/transcriptome"

#########################
#########################
#########################


# Auxiliary variables ---------------------------------------------------------------------------------- 

tr -d '\r' <${dataset_table_aux} >${dataset_table} #format to unix 
kallisto_index="kallisto_index"

# ----------------------------------------- CREATE DIRECTORIES -----------------------------------------

mkdir ${working_directory}/FASTQ

mkdir ${working_directory}/fasterqdump_logs

mkdir ${working_directory}/kallisto

mkdir ${working_directory}/kallisto_index_log

mkdir ${working_directory}/FASTQC

mkdir ${working_directory}/kallisto_logs_samples

mkdir ${working_directory}/MULTIQC

mkdir ${working_directory}/temp

######################################################### PIPELINE ########################################################################
 

 echo "####################### STARTING PIPELINE... ####################### "

 echo "--> Starting time: " $(date +%F_%T)

# ----------------------------------------- Kallisto index -----------------------------------------

echo "Creating Kallisto index..."

kallisto index -i ${working_directory}/${kallisto_index} ${transcriptome} > ${working_directory}/kallisto_index_log/kallisto_index_log.txt 2>&1
 

# ----------------------------------------- Pipeline for each sample -----------------------------------------

for sample in $(awk '{print $1}' ${dataset_table})
do 

  (

  # ----------------------------------------- DOWNLOAD SAMPLE -----------------------------------------

  # --split-files splits the FASTQ reads into two files, when its paired-end
  # https://rnnh.github.io/bioinfo-notebook/docs/fasterq-dump.html

  echo "Downloading sample " ${sample} "..." 
  fasterq-dump ${sample} --temp ${working_directory}/temp --split-files --threads ${threads} --outdir ${working_directory}/FASTQ/ > ${working_directory}/fasterqdump_logs/log_fasterqdump_${sample}.txt 2>&1
  
  nb_files=$(ls ${working_directory}/FASTQ/${sample}* | wc -l)

  # if download failed, there's no point in doing the rest of the pipeline
  if [ $nb_files -eq 0 ]; then
    
    echo "Failed download for sample " ${sample}

    echo ${sample} >> ${working_directory}/failedSamples_${dataset_description}.txt

  else


      echo "Compressing samples " ${sample} "..."
      gzip -f ${working_directory}/FASTQ/${sample}*

        # ----------------------------------------- CHECK IF SINGLE OR PAIRED-END -----------------------------------------

         
      echo "  Number of files: " ${nb_files}

      if [ $nb_files -ge 2 ]; then
         paired_end=true
         # Remove unmated reads, so that paired-end has only two samples
         rm -f ${working_directory}/FASTQ/${sample}.*

         echo "  info: Sample " ${sample} " is paired-end"
      else
        paired_end=false
        echo "  info: Sample " ${sample} " is single-end"
      fi


      # ------------------------------------------------ PERFORM FASTQC ------------------------------------------------

    echo "Performing FASTQC for sample " ${sample} "..."

     # Check the value of paired_end and execute different code accordingly
      if [ "$paired_end" = true ]; then

        fastqc --quiet ${working_directory}/FASTQ/${sample}_1* -o ${working_directory}/FASTQC
        fastqc --quiet ${working_directory}/FASTQ/${sample}_2* -o ${working_directory}/FASTQC

      else
        
        fastqc --quiet ${working_directory}/FASTQ/${sample}* -o ${working_directory}/FASTQC

      fi


      # ------------------------------------------------ ALIGN WITH KALLISTO ------------------------------------------------

      echo "Kallisto alignment for sample " ${sample} "..."

       # Check the value of paired_end and execute different code accordingly
      if [ "$paired_end" = true ]; then

        kallisto quant -i ${working_directory}/${kallisto_index} -o ${working_directory}/kallisto/${sample} ${working_directory}/FASTQ/${sample}_1.* ${working_directory}/FASTQ/${sample}_2.* -t ${threads} --bias --single-overhang > ${working_directory}/kallisto_logs_samples/log_kallisto_quant_${sample}.txt 2>&1

      else
        
        kallisto quant -i ${working_directory}/${kallisto_index} -o ${working_directory}/kallisto/${sample} --single -l 200 -s 20 ${working_directory}/FASTQ/${sample}.* -t ${threads} --bias --single-overhang > ${working_directory}/kallisto_logs_samples/log_kallisto_quant_${sample}.txt 2>&1

      fi

      # ------------------------------------------------ REMOVE SAMPLE ------------------------------------------------

      echo "Removing sample " ${sample} "..."

      rm ${working_directory}/FASTQ/${sample}* 
     

  fi
  
  ) & # to parallelise 
 
     # allow to execute up to $N jobs in parallel
    if [ $(jobs -r -p | wc -l) -ge $N ]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
 

done


# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait



# ------------------------------------------------ PERFORM MULTIQC BY DATASET ------------------------------------------------

# perform multiqc for each dataset in separate, including kallisto logs

# Iterate over each dataset

for dataset in $(awk '{print $2}' ${dataset_table} | sort -u); do

    (

    echo "Performing multiQC for dataset " ${dataset} "..."

    # Extract sample IDs for the current dataset

    samples_in_batch=$(awk -v dataset="$dataset" '$2 == dataset {print $1}' "$dataset_table" )
    #echo ${samples_in_batch} 

    full_paths_fastqc=()
    full_paths_kallisto=()
        # Iterate over each sample ID
    for sample in $samples_in_batch; do
        full_paths_fastqc+=" $(ls ${working_directory}/FASTQC/${sample}*fastqc.zip )"
        full_paths_kallisto+=" $(ls ${working_directory}/kallisto_logs_samples/log_kallisto_quant_${sample}*.txt )"
    done

    # Run MultiQC on the samples in the current batch
    multiqc -q  --interactive --title ${dataset} --module fastqc ${full_paths_fastqc} --module kallisto ${full_paths_kallisto} -o ${working_directory}/MULTIQC/$dataset

    ) &

    # allow to execute up to $N jobs in parallel
    if [ $(jobs -r -p | wc -l) -ge $N ]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
 

done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

 
echo "####################### DONE! ####################### "
echo "--> Ending time: " $(date +%F_%T)
