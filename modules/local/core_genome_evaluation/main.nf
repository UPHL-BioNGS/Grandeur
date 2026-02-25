process CORE_GENOME_EVALUATION {
  tag         "Evaluating core genome"
  label       "process_single"
  container   'staphb/pandas:2.3.0'

  input:
  tuple file(fasta), file(summary), file(script)

  output:
  path("core_genome_values.csv"), emit: evaluation
  path "core_genome_evaluation/core_genome_evaluation.csv", emit: for_multiqc
  path "logs/${task.process}/*.log"                       , emit: log_files

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p core_genome_evaluation logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    python3 ${script} | tee -a \$log_file

    num_samples=\$(wc -l core_genome_evaluation.csv | awk '{print \$1}' )
    num_core_genes=\$(cut -f 3 core_genome_evaluation.csv -d "," | tail -n 1 | cut -f 1 -d "." )
    per_core_genes=\$(cut -f 7 core_genome_evaluation.csv -d "," | tail -n 1 )
    echo "\$num_samples,\$num_core_genes,\$per_core_genes" > core_genome_values.csv
    cp core_genome_evaluation.csv core_genome_evaluation/core_genome_evaluation.csv
  """
}
