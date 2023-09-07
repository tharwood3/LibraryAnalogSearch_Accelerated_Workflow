import os
import sys
import argparse
import logging
import multiprocessing as mp
from timeit import default_timer as timer
import pandas as pd
import numpy as np
from functools import partial
from blink import blink
from matchms import Spectrum
from matchms.similarity import ModifiedCosine

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

blink_to_gnps_map = {'query_filename':'SpectrumFile', 'query_id':'#Scan#', 'score':'MQScore', 'ref_id':'LibrarySpectrumID', 'ref_filename':'LibraryName', 
'query_charge':'Charge', 'query_precursor_mz':'SpecMZ', 'mz_ppm_error':'mzErrorPPM', 'mz_diff':'ParentMassDiff', 'query_rt':'p-value', 'matches':'LibSearchSharedPeaks'}

def compute_nl_sum(S12):
    sum_scores = S12['mdi'] + S12['nli']
    sum_counts = S12['mdc'] + S12['nlc']
    
    return {'mzi':sum_scores, 'mzc':sum_counts}

def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)

def create_blink_parser():

    parser = argparse.ArgumentParser(description='REM-BLINK Efficiently Performs Analog Searches Accounting For Multiple Mass Differences')

    parser.add_argument('query_file', action='store', type=file_path, metavar='Q', help='MGF or mzML files containing experimental spectra')
    parser.add_argument('reference_file', action='store', type=file_path, metavar='R', help='MGF or mzML files containing reference spectra')
    parser.add_argument('output_file', action='store', type=str, metavar='O', help='path to output file. Output file is a CSV')
    parser.add_argument('polarity', action='store', type=str, metavar='P', choices=['positive', 'negative'], help='ion mode of query spectra. Allowed inputs as "positive" and "negative"')

    #Discretize options
    discretize_options = parser.add_argument_group()
    discretize_options.add_argument('-t', '--tolerance', action='store', type=float, default=0.01, required=False,
                                    help='allowed dalton tolerance between matched MS/MS fragment peaks')
    
    discretize_options.add_argument('-b','--bin_width', action='store', type=float, default=0.001, required=False,
                                 help='width of bins in m/z. Larger bins will be faster at the expense of precision.')
    
    discretize_options.add_argument('-i','--intensity_power', action='store', type=float, default=0.5, required=False,
                                 help='exponent used to adjust intensities prior to scoring')
    
    discretize_options.add_argument('--trim', action='store_true', default=False, required=False,
                                    help='remove empty spectra when discretizing')
    
    discretize_options.add_argument('--dedup', action='store_true', default=False, required=False,
                                    help='deduplicate fragment ions within 2 times bin_width')

    #Filtering options
    filter_options = parser.add_argument_group()

    filter_options.add_argument('-s','--min_score', type=float, default=0.2, required=False,
                                 help='minimum cosine score to include in output.')
    filter_options.add_argument('-m','--min_matches', type=int, default=3, required=False,
                                 help='minimum matches to include in output.')
    # filter_options.add_argument('-o', '--override_matches', type=int, default=5, required=False,
    #                             help='number of matches to keep comparison regardless of score')
    
    #GNPS options
    gnps_options = parser.add_argument_group()

    gnps_options.add_argument('--max_shift', type=float, default=200.0, required=False,
                              help='maximum precursor m/z difference to use for scoring')
    gnps_options.add_argument('--num_cores', type=int, default=4, required=False,
                              help='Number of cores used during multiprocessing steps')

    return parser

def mms_score_match(query_pmzs, reference_pmzs,
                    query_spectra, reference_spectra, 
                    mod, args, pair):
    
    query_pmz = query_pmzs[pair[0]]
    lib_pmz = reference_pmzs[pair[1]]
    
    pmz_diff = abs(query_pmz - lib_pmz)
    
    if pmz_diff > args.max_shift:
        return None, None
    
    query_spec = query_spectra[pair[0]]
    lib_spec = reference_spectra[pair[1]]

    query_spectrum = Spectrum(query_spec[0], query_spec[1], metadata={'precursor_mz':query_pmz})
    lib_spectrum = Spectrum(lib_spec[0], lib_spec[1], metadata={'precursor_mz':lib_pmz})
    res = mod.pair(query_spectrum, lib_spectrum)
    
    score = res['score'].item()
    matches = res['matches'].item()
    
    if score < args.min_score or matches < args.min_matches:
        return None, None
    
    return matches, score

def main():

    parser = create_blink_parser()
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    start = timer()
    query_df = blink.open_msms_file(args.query_file)
    reference_df = blink.open_msms_file(args.reference_file)
    end = timer()

    if args.polarity == 'positive':
        if 'ionmode' in reference_df.columns:
            reference_df = reference_df[reference_df['ionmode'].str.lower() != 'negative']

    else:
        if 'ionmode' in reference_df.columns:
            reference_df = reference_df[reference_df['ionmode'].str.lower() != 'positive']

    logging.info('Input files read time: {} seconds, {} spectra'.format(end-start, query_df.shape[0]+reference_df.shape[0]))

    query_spectra = query_df.spectrum.tolist()
    reference_spectra = reference_df.spectrum.tolist()

    query_pmzs = query_df.precursor_mz.tolist()
    reference_pmzs = reference_df.precursor_mz.tolist()

    start = timer()
    blink_prefilter_tol = args.tolerance * 3
    discretized_spectra = blink.discretize_spectra(query_spectra, reference_spectra, query_pmzs, reference_pmzs, network_score=True, 
                                             tolerance=blink_prefilter_tol, bin_width=args.bin_width, intensity_power=args.intensity_power,
                                             trim_empty=args.trim, remove_duplicates=args.dedup)
    end = timer()
    
    logging.info('Discretization time: {} seconds, {} spectra'.format(end-start, len(query_spectra)+len(reference_spectra)))
    
    start = timer()
    scores = blink.score_sparse_spectra(discretized_spectra)
    scores = compute_nl_sum(scores)
    end = timer()

    logging.info('Scoring time: {} seconds, {} comparisons'.format(end-start, len(query_spectra)*len(reference_spectra)))
    
    start = timer()
    filtered_hits = blink.filter_hits(scores, min_score=args.min_score, min_matches=args.min_matches, override_matches=None)
    end = timer()

    logging.info('Filtering time: {} seconds, {} comparisons'.format(end-start, len(query_spectra)*len(reference_spectra)))

    query_idxs, library_idxs = filtered_hits['mzi'].nonzero()
    filtered_pairs = np.array((query_idxs, library_idxs)).T

    start = timer()
    mod = ModifiedCosine(tolerance=args.tolerance)
    pool = mp.Pool(args.num_cores)
    score_match = partial(mms_score_match, query_pmzs, reference_pmzs,
                                   query_spectra, reference_spectra,
                                   mod, args)
    matchms_hits = pool.map(score_match, filtered_pairs)
    end = timer()

    logging.info('MatchMS time: {} seconds, {} comparisons'.format(end-start, filtered_pairs.shape[0]))

    output = pd.concat([pd.DataFrame(matchms_hits, columns=['matches', 'score']), pd.DataFrame(filtered_pairs, columns=['query', 'ref'])], axis=1)

    output = output.astype({'query':int, 'ref':int})

    output['query_filename'] = os.path.basename(args.query_file)
    output['ref_filename'] = os.path.basename(args.reference_file)

    output = pd.merge(output, query_df['precursor_mz'], left_on='query', right_index=True)
    output = pd.merge(output, query_df['charge'], left_on='query', right_index=True)
    output.rename(columns={'precursor_mz':'query_precursor_mz', 'charge':'query_charge'}, inplace=True)

    output = pd.merge(output, reference_df['precursor_mz'], left_on='ref', right_index=True)
    output = pd.merge(output, reference_df['charge'], left_on='ref', right_index=True)
    output.rename(columns={'precursor_mz':'ref_precursor_mz', 'charge':'ref_charge'}, inplace=True)

    #mzML Query
    if 'id' in query_df.columns:
        output = pd.merge(output, query_df['id'], left_on='query', right_index=True)
        output = pd.merge(output, (query_df['rt'] * 60), left_on='query', right_index=True)
        output.rename(columns={'id':'query_id', 'rt':'query_rt'}, inplace=True)

    #FBMN MGF Query
    elif 'scans' in query_df.columns:
        output = pd.merge(output, query_df['scans'], left_on='query', right_index=True)
        output = pd.merge(output, query_df['rtinseconds'], left_on='query', right_index=True)
        output.rename(columns={'scans':'query_id', 'rtinseconds':'query_rt'}, inplace=True)
    
    #Lib MGF Query
    elif 'spectrumid' in query_df.columns:
        output = pd.merge(output, query_df['spectrumid'], left_on='query', right_index=True)   
        output['query_rt'] = 0
        output.rename(columns={'spectrumid':'query_id'}, inplace=True)

    #mzML Ref
    if 'id' in reference_df.columns:
        output = pd.merge(output, reference_df['id'], left_on='ref', right_index=True)
        output.rename(columns={'id':'ref_id'}, inplace=True)
        
    #Lib MGF Ref
    elif 'spectrumid' in reference_df.columns:
        output = pd.merge(output, reference_df['spectrumid'], left_on='ref', right_index=True)
        output.rename(columns={'spectrumid':'ref_id'}, inplace=True)

    output['mz_diff'] = abs(output['query_precursor_mz'] - output['ref_precursor_mz'])
    output['mz_ppm_error'] = output['mz_diff'] / output['ref_precursor_mz'] * 1000000

    output.rename(columns=blink_to_gnps_map, inplace=True)
    output["FileScanUniqueID"] = output["SpectrumFile"].astype(str) + ":" + output["#Scan#"].astype(str)
    
    #boilerplate
    output["UnstrictEvelopeScore"] = 0

    start = timer()
    output[output['MQScore'].notna()].drop(columns=['query', 'ref']).to_csv(args.output_file)
    end = timer()

    logging.info('Output write time: {} seconds, {} rows'.format(end-start, output.shape[0]))

if __name__ == "__main__":
    main()