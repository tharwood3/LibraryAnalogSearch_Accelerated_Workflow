#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.inputlibraries = "data/libraries"
params.inputspectra = "data/spectra"

// Parameters
params.topk = 1

params.fragment_tolerance = 0.05

params.library_min_cosine = 0.7
params.library_min_matched_peaks = 6

params.analog_max_shift = 400

// Blink Parameters
params.blink_ionization = "positive"
params.blink_binwidth = 0.01
params.blink_numcores = 4

TOOL_FOLDER = "$baseDir/bin"

process searchDataBLINK {
    //publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each file(input_library)
    each file(input_spectrum)

    output:
    file 'search_results/*.csv' optional true

    script:
    def randomFilename = UUID.randomUUID().toString()
    def input_spectrum_abs = input_spectrum.toRealPath()
    def input_library_abs = input_library.toRealPath()

    """
    mkdir search_results
    echo $workDir
    previous_cwd=\$(pwd)
    echo \$previous_cwd

    python $TOOL_FOLDER/blink_analog_prefilter.py \
    $input_spectrum_abs \
    $input_library_abs \
    \$previous_cwd/search_results/${randomFilename}.csv \
    $params.blink_ionization \
    --tolerance $params.fragment_tolerance \
    --bin_width $params.blink_binwidth \
    --min_score $params.library_min_cosine \
    --min_matches $params.library_min_matched_peaks \
    --max_shift $params.analog_max_shift \
    --num_cores $params.blink_numcores \
    """
}

process mergeResults {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "results/*"

    output:
    path 'merged_results.tsv'

    """
    python $TOOL_FOLDER/tsv_merger.py \
    results \
    merged_results.tsv \
    --topk $params.topk
    """
}

process getGNPSAnnotations {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "merged_results.tsv"

    output:
    path 'merged_results_with_gnps.tsv'

    """
    python $TOOL_FOLDER/getGNPS_library_annotations.py \
    merged_results.tsv \
    merged_results_with_gnps.tsv
    """

}

workflow {
    libraries = Channel.fromPath(params.inputlibraries + "/*.mgf" )
    spectra = Channel.fromPath(params.inputspectra + "/**" )
    
    search_results = searchDataBLINK(libraries, spectra)
    
    merged_results = mergeResults(search_results.collect())

    getGNPSAnnotations(merged_results)
}