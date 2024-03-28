import peptacular as pt
from typing import List
import re
import pandas as pd
import sys
import argparse
import os

def parse_args(args: List[str]) -> dict:

    def writable_file_path(path):
        if os.path.exists(path) and not os.access(path, os.W_OK):
            raise argparse.ArgumentTypeError(f"File {path} is not writable.")
        return path

    def readable_file_path(path):
        if not os.path.exists(path) or not os.access(path, os.R_OK):
            raise argparse.ArgumentTypeError(f"File {path} does not exist or is not readable.")
        return path

    parser = argparse.ArgumentParser(description='Fragment a peptide sequence.')

    # Required arguments
    parser.add_argument('spectra',
                        type=readable_file_path,
                        help='The spectra file to read. Should be a plain text file with mz and intensity values '
                             'separated by a space.')
    parser.add_argument('sequence',
                        type=str,
                        help='The peptide sequence to fragment (Proforma2.0 format).')

    # Options for processing
    parser.add_argument('--fragment_types',
                        type=str,
                        nargs='+',
                        help='The fragment types to generate in the format: {charge}{ion_type}. '
                             'Supported special ions: [i, p] (i = Immonium, p = precursor)], '
                             'Supported terminal ions: [a, b, c, x, y, z], '
                             'Supported internal ions: [ax, ay, az, bx, by, bz, cx, cy, cz].',
                        default=['1b', '1y'])
    processing_group = parser.add_argument_group('Processing Options')
    processing_group.add_argument('--monoisotopic',
                                  action='store_true',
                                  help='Use monoisotopic mass.')
    processing_group.add_argument('--isotopes',
                                  type=int,
                                  nargs='+',
                                  help='The isotopes to consider, include 0 for base peak.',
                                  default=[0])
    processing_group.add_argument('--losses',
                                  type=int,
                                  nargs='+',
                                  help='The losses to consider, include 0 for no loss.',
                                  default=[0])

    # Options for error tolerance and peak assignment
    tolerance_group = parser.add_argument_group('Tolerance and Peak Assignment')
    tolerance_group.add_argument('--error_tolerance',
                                 type=float,
                                 help='The error tolerance to consider when matching peaks.',
                                 default=50)
    tolerance_group.add_argument('--error_tolerance_type',
                                 choices=['ppm', 'th'],
                                 help='The type of error tolerance to consider when matching peaks.',
                                 default='ppm')
    tolerance_group.add_argument('--peak_assignment',
                                 choices=['closest', 'largest'],
                                 help='The method to use when assigning peaks to fragments.',
                                 default='closest')

    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--fragment_output',
                              type=writable_file_path,
                              help='The output file to write the fragment data to.',
                              default='frags.csv')
    output_group.add_argument('--annotated_spectra_output',
                              type=writable_file_path,
                              help='The output file to write the annotated spectra data to.',
                              default='annots.csv')

    return vars(parser.parse_args(args))


def get_fragments(sequence: str, fragment_types: List[str], monoisotopic: bool, isotopes: List[int] | None,
                  losses: List[int] | None) -> List[pt.Fragment]:

    if isotopes is None:
        isotopes = [0]

    if losses is None:
        losses = [0]

    fragments = []

    fragmenter = pt.Fragmenter(sequence=sequence, monoisotopic=monoisotopic)
    for fragment_type in fragment_types:

        # Use regex to extract the charge and ion type. Assuming the ion type follows the charge in the string.
        charge = int(re.search(r'^\d+', fragment_type).group())
        ion_type = re.search(r'[a-z]+$', fragment_type, re.I).group()


        fragments.extend(fragmenter.fragment(ion_types=ion_type, charges=charge, isotopes=isotopes, losses=losses))

    return fragments

def get_annotations(mz_values: List[float], intensity_values: List[float], fragments: List[pt.Fragment],
                    mass_tolerance: float, mass_tolerance_type: str, peak_assignment:str) -> List[pt.FragmentMatch]:

    fragment_matches = pt.get_fragment_matches(fragments,
                                               mz_values,
                                               intensity_values,
                                               mass_tolerance,
                                               mass_tolerance_type,
                                               'largest' if peak_assignment == 'most intense' else 'closest')

    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)
    fragment_matches = {fm.mz: fm for fm in fragment_matches}  # keep the best fragment match for each mz

    return fragment_matches


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])


    mz_values, intensity_values = [], []
    with open(args['spectra'], 'r') as f:
        for line in f:
            mz, intensity = line.strip().split(' ')
            mz_values.append(float(mz))
            intensity_values.append(float(intensity))

    fragments = get_fragments(sequence=args['sequence'],
             fragment_types=args['fragment_types'],
             monoisotopic=args['monoisotopic'],
             isotopes=args['isotopes'],
             losses=args['losses'])

    frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])
    frag_df = frag_df.round(5)
    frag_df['start'] = frag_df['start'].astype(int)
    frag_df['end'] = frag_df['end'].astype(int)
    frag_df['isotope'] = frag_df['isotope'].astype(int)
    frag_df['charge'] = frag_df['charge'].astype(int)


    fragment_matches = get_annotations(mz_values,
                               intensity_values,
                               fragments,
                               args['error_tolerance'],
                               args['error_tolerance_type'],
                               args['peak_assignment'])


    if len(fragment_matches) == 0:

        empty_frag_match = pt.FragmentMatch(None, 0, 0)
        cols = list(empty_frag_match.to_dict().keys())

        data = []
        for mz, intensity in zip(mz_values, intensity_values):
            data.append({'mz': mz, 'intensity': intensity, 'matched': False})

        spectra_df = pd.DataFrame(data)

        # add other cols as nan
        for col in cols:
            if col not in spectra_df.columns:
                spectra_df[col] = None


    else:
        match_data = []
        data = []
        for mz, intensity in zip(mz_values, intensity_values):
            fm = fragment_matches.get(mz, None)
            if fm:
                match_data.append(fm.to_dict())
            else:
                data.append({'mz': mz, 'intensity': intensity})

        unassigned_df = pd.DataFrame(data)
        unassigned_df['matched'] = False

        match_df = pd.DataFrame(match_data)
        match_df['matched'] = True
        match_df['abs_error'] = match_df['error'].abs()
        match_df['abs_error_ppm'] = match_df['error_ppm'].abs()

        # merge dfs
        spectra_df = pd.concat([match_df, unassigned_df], ignore_index=True)


    spectra_df = spectra_df.sort_values(by='mz')
    spectra_df = spectra_df.round(5)
    spectra_df['start'] = spectra_df['start'].astype('Int16')
    spectra_df['end'] = spectra_df['end'].astype('Int16')
    spectra_df['isotope'] = spectra_df['isotope'].astype('Int16')
    spectra_df['charge'] = spectra_df['charge'].astype('Int16')


    # make macthed column be 3rd column
    cols = spectra_df.columns.tolist()
    cols.remove('matched')
    cols.insert(2, 'matched')
    spectra_df = spectra_df[cols]

    if args['fragment_output']:
        frag_df.to_csv(args['fragment_output'], index=False)


    if args['annotated_spectra_output']:
        spectra_df.to_csv(args['annotated_spectra_output'], index=False)
    else:
        print(spectra_df.to_string())


