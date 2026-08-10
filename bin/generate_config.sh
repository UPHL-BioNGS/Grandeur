#!/bin/bash

# Define the output configuration file name
OUTPUT_FILE="edit_me.config"

# 1. Write the static header and parameters to the file
cat << 'EOF' > "$OUTPUT_FILE"
//# Instructions --------------------------------------------
//# We think it's great that you want to adjust some parameters.
//# This is very useful when running this on the cloud.
//# This is especially useful for the following processes:
//#     - params.kraken2_db to specify where the kraken2 database is on your system
//# Right now, most everything is commented out with '//'.
//# To adjust a parameter, remove the '//' from in front of a param and replace the value
//# with the one that suits your needs.
//# Don't hesitate to ask for assistance at https://github.com/UPHL-BioNGS/Grandeur/issues
//# if something doesn't work (also, please include a copy of your config file).

//# Docker Params -------------------------------------------
//docker.enabled                  = true
//docker.runOptions               = '-u $(id -u):$(id -g)'
//docker.sudo                     = false
//docker.temp                     = /tmp
//docker.remove                   = true
//docker.fixOwnership             = true
//docker.engineOptions            = ''
//docker.mountFlags               = ''

//# Singularity Params --------------------------------------
//singularity.enabled             = true
//singularity.autoMounts          = true
//singularity.runOptions          = ""
//process.stageInMode             = "link"
//singularity.engineOptions       = ''
//singularity.cacheDir            = ''

//# AWS Batch Params ----------------------------------------
//process.executor                = 'awsbatch'
//process.queue                   = ''
//aws.batch.cliPath               = '/home/ec2-user/miniconda/bin/aws'
//aws.region                      = 'us-east-1'
//workDir                         = 's3://'

//# Google Cloud Params -------------------------------------
//process.executor                = ''
//google.project                  = ''
//google.location                 = ''
//google.region                   = ''
//workDir                         = ''
//google.lifeSciences.bootDiskSize = 50.GB

//# HPC / Cluster Params ------------------------------------
//process.executor                = 'slurm'
//process.queue                   = 'general'
//process.clusterOptions          = '--account=my_lab_account'
//executor.queueSize              = 100
//executor.submitRateLimit        = '10 sec'

//# Nextflow Tower ------------------------------------------
//tower.accessToken               = ''
//tower.enabled                   = true

//# Nextflow Reports & Tracing ------------------------------
//report.enabled                  = true
//report.file                     = "${params.outdir}/pipeline_info/execution_report.html"
//timeline.enabled                = true
//timeline.file                   = "${params.outdir}/pipeline_info/execution_timeline.html"
//trace.enabled                   = true
//trace.file                      = "${params.outdir}/pipeline_info/execution_trace.txt"
//dag.enabled                     = true
//dag.file                        = "${params.outdir}/pipeline_info/pipeline_dag.svg"

//# Disk Management -----------------------------------------
//cleanup                         = true

//# How this file was obtained -------------------------------
//params.config_file              = false

//# Adjustable Workflow parameters ---------------------------
EOF

# 2. Extract workflow parameters dynamically from nextflow_schema.json
if [[ -f "nextflow_schema.json" ]]; then
    # Look for .properties, .definitions[].properties, and ."$defs"[].properties
    jq -r '
      [ .properties // {}, .["$defs"][]?.properties // {}, .definitions[]?.properties // {} ]
      | add
      | to_entries[]
      | .key as $k
      | .value as $v
      | (if $v | has("default") then $v.default else null end) as $d
      | (
          if $d == null then "null"
          elif $v.type == "string" then "'\''\($d)'\''"
          else ($d | tojson)
          end
        ) as $val
      | "//params.\($k)\t\($val)"
    ' nextflow_schema.json | while IFS=$'\t' read -r param val; do
        printf "%-31s = %s\n" "$param" "$val" >> "$OUTPUT_FILE"
    done
else
    echo "//# WARNING: nextflow_schema.json not found in the current directory." >> "$OUTPUT_FILE"
fi

# 3. Dynamically add the process block header from conf/base.config
echo "" >> "$OUTPUT_FILE"
if [[ -f "conf/base.config" ]]; then
    # Delete the last line ($d) and prepend // to comment it out
    sed '$d' conf/base.config | grep -v "//" | grep -A 100000000 "process {" | sed 's/^/\/\//'  >> "$OUTPUT_FILE"
else
    echo "//# WARNING: conf/base.config not found." >> "$OUTPUT_FILE"
    echo "//process {" >> "$OUTPUT_FILE"
fi

# 4. Iterate through all main.nf files and parse out the process details
for nf_file in modules/local/*/main.nf; do
    # Skip if no files are found
    if [[ ! -f "$nf_file" ]]; then
        continue
    fi

    # Extract the process name
    proc_name=$(grep -E '^process ' "$nf_file" | awk '{print $2}' | tr -d '{')
    
    # Extract configuration keys
    tag=$(grep -E '^[[:space:]]*tag[[:space:]]+' "$nf_file" | sed -E 's/^[[:space:]]*tag[[:space:]]+//')
    label=$(grep -E '^[[:space:]]*label[[:space:]]+' "$nf_file" | sed -E 's/^[[:space:]]*label[[:space:]]+//')
    container=$(grep -E '^[[:space:]]*container[[:space:]]+' "$nf_file" | sed -E 's/^[[:space:]]*container[[:space:]]+//')

    # Extract explicit time and errorStrategy if they exist in main.nf
    time_val=$(grep -m 1 -E '^[[:space:]]*time[[:space:]=]+' "$nf_file" | sed -E 's/^[[:space:]]*time[[:space:]=]+//')
    error_val=$(grep -m 1 -E '^[[:space:]]*errorStrategy[[:space:]=]+' "$nf_file" | sed -E 's/^[[:space:]]*errorStrategy[[:space:]=]+//')
    
    # Clean the label (strip quotes and spaces) to accurately map resources
    label_clean=$(echo "$label" | tr -d "'\"[:space:]")
    
    # Determine cpus and memory based on base.config mapping
    case "$label_clean" in
        process_single)      cpus="1";  memory="6.GB" ;;
        process_low)         cpus="2";  memory="12.GB" ;;
        process_medium)      cpus="6";  memory="36.GB" ;;
        process_high)        cpus="12"; memory="72.GB" ;;
        process_high_memory) cpus="1";  memory="200.GB" ;;
        process_long)        cpus="1";  memory="6.GB" ;; # Using defaults for CPU/mem
        *)                   cpus="1";  memory="6.GB" ;; # Fallback defaults
    esac

    # Determine final values to write (favoring main.nf if it existed)
    final_time="${time_val:-$default_time}"
    final_error="${error_val:-{ task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }}"

    # Extract default ext arguments and prefix from the script block
    args_val=$(grep -E 'def args[[:space:]]*=[[:space:]]*task\.ext\.args[[:space:]]*\?:' "$nf_file" | sed -E 's/.*\?:[[:space:]]*//')
    prefix_val=$(grep -E 'def prefix[[:space:]]*=[[:space:]]*task\.ext\.prefix[[:space:]]*\?:' "$nf_file" | sed -E 's/.*\?:[[:space:]]*//')
    
    # Append the parsed block to the output file (now including cpus and memory)
    cat << EOF >> "$OUTPUT_FILE"
//    withName: ${proc_name} {
//        tag           = ${tag}
//        label         = ${label}
//        cpus          = ${cpus}
//        memory        = ${memory}
//        publishDir = [
//            path: { "\${params.outdir}" },
//            mode: params.publish_dir_mode,
//            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
//        ]
//        container     = ${container}
//        errorStrategy = ${final_error}
//        time          = ${final_time}
EOF

    # Conditionally append args if they exist
    if [[ -n "$args_val" ]]; then
        formatted_args=$(echo "$args_val" | sed "s/'/\"/g")
        echo "//        ext.args      = ${formatted_args}" >> "$OUTPUT_FILE"
    fi

    # Conditionally append prefix if it exists
    if [[ -n "$prefix_val" ]]; then
        echo "//        ext.prefix    = ${prefix_val}" >> "$OUTPUT_FILE"
    fi

    # Close the individual withName block
    echo "//    }" >> "$OUTPUT_FILE"

done

# 5. Close the global process block
echo "//}" >> "$OUTPUT_FILE"

echo "Configuration file successfully generated at: $OUTPUT_FILE"
echo "Use cp  $OUTPUT_FILE conf/grandeur_template.config"