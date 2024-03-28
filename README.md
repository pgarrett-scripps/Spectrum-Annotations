### How to Run:

```bash
python main.py spectra.ms2 [164.0700]-FDSFGDLSSASAIM[16]GNPK --error_tolerance 50 --error_tolerance_type th --isotopes 0 1 --losses 0 18
```

### Command Line Arguments:

```
usage: main.py [-h] [--fragment_types FRAGMENT_TYPES [FRAGMENT_TYPES ...]]
               [--monoisotopic] [--isotopes ISOTOPES [ISOTOPES ...]]
               [--losses LOSSES [LOSSES ...]]
               [--error_tolerance ERROR_TOLERANCE]
               [--error_tolerance_type {ppm,th}]
               [--peak_assignment {closest,largest}]
               [--fragment_output FRAGMENT_OUTPUT]
               [--annotated_spectra_output ANNOTATED_SPECTRA_OUTPUT]
               spectra sequence

Fragment a peptide sequence.

positional arguments:
  spectra               The spectra file to read. Shoul dbe a plain text file   
                        with mz and intensity values separated by a space.      
  sequence              The peptide sequence to fragment (Proforma2.0 format).  

options:
  -h, --help            show this help message and exit
  --fragment_types FRAGMENT_TYPES [FRAGMENT_TYPES ...]
                        The fragment types to generate in the format:
                        {charge}{ion_type}. Supported Ions: a, b, c, x, y, z,   
                        i, p. (i = Immonium, p = precursor)

Processing Options:
  --monoisotopic        Use monoisotopic mass.
  --isotopes ISOTOPES [ISOTOPES ...]
                        The isotopes to consider, include 0 for base peak.      
  --losses LOSSES [LOSSES ...]
                        The losses to consider, include 0 for no loss.

Tolerance and Peak Assignment:
  --error_tolerance ERROR_TOLERANCE
                        The error tolerance to consider when matching peaks.    
  --error_tolerance_type {ppm,th}
                        The type of error tolerance to consider when matching   
                        peaks.
  --peak_assignment {closest,largest}
                        The method to use when assigning peaks to fragments.    

Output Options:
  --fragment_output FRAGMENT_OUTPUT
                        The output file to write the fragment data to.
  --annotated_spectra_output ANNOTATED_SPECTRA_OUTPUT
                        The output file to write the annotated spectra data     
                        to.


```
"# Spectrum-Annotations" 
