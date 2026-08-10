process CONCAT_REPORTS {
    tag           "${output_name}"
    container     'staphb/pandas:3.0.5'
    label         "process_single"

    input:
    tuple file(files), val(output_name), val(subdir), val(keep_header)

    output:
    path "*/${output_name}", emit: summary

    script:
    if ( keep_header )
        """
        mkdir -p ${subdir}
        head -n 1 ${files[0]} > ${subdir}/${output_name}
        tail -q -n +2 ${files} | sort >> ${subdir}/${output_name}
        """
    else
        """
        mkdir -p ${subdir}
        cat ${files} > ${subdir}/${output_name}
        """
}