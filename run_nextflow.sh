genome_path='/mnt/galaxy/galaxyworks/tool-data/T2T_chm13v2.0+bs_controls/bwameth_index/T2T_chm13v2.0+bs_controls/T2T_chm13v2.0+bs_controls.fa'
             # mnt/galaxy/data/genome/grch38_hs38d1_noalt+bs_controls/bwameth_index/grch38+bs_controls/grch38_core+bs_controls.fa'

. /mnt/home/aerijman/miniconda3/bin/activate
conda activate nextflow

nextflow main.nf  \
  --input_glob "test_data/*1.fastq.gz" \
  --genome "${genome_path}" \
  --email "aerijman@neb.com" \
  --max_input_reads 10000 \
  --flowcell "test_pipeline" \
  -with-report  "emseq_metadata_report.html" \
  -with-timeline "emseq_metadata_timeline.html" \
  -with-dag "emseq_metadata_dag.html" \
  -w "/mnt/hpc_scratch/" \
  --read_length 151 \
  -resume
