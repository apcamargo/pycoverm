use coverm::{
    bam_generator::*, contig::*, coverage_takers::*, mosdepth_genome_coverage_estimators::*,
    FlagFilter,
};
use ndarray::Array;
use numpy::convert::ToPyArray;
use pyo3::{prelude::*, wrap_pyfunction};
use rust_htslib::{bam, bam::Read};
use std::collections::HashSet;
use std::str;

struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}

struct EstimatorsAndTaker {
    estimators: Vec<CoverageEstimator>,
    taker: CoverageTakerType,
}

/// is_bam_sorted(bam_file)
/// --
///
/// Checks whether a BAM file is sorted by coordinate.
///
/// Parameters
/// ----------
/// bam_file : str
///     Path to a BAM file.
///
/// Returns
/// -------
/// bool
///     Returns `True` if the BAM file is sorted by coordinate and `False`
///     otherwise.
#[pyfunction]
fn is_bam_sorted(bam_file: &str) -> PyResult<bool> {
    let bam = bam::Reader::from_path(bam_file).unwrap();
    let bytes = &bam::Header::from_template(bam.header()).to_bytes();
    Ok(bytes.windows(13).any(|window| window == b"SO:coordinate"))
}

/// get_coverages_from_bam(bam_list, contig_end_exclusion=75, min_identity=0.97,
/// trim_lower=0.0, trim_upper=0.0, contig_list=None, threads=1)
/// --
///
/// Computes contig mean coverages from sorted BAM files. All BAM files must be
/// mapped to the same reference.
/// Trimmed means will be computed if `trim_min` and/or `trim_max` are set to
/// values greater than 0.
///
/// Parameters
/// ----------
/// bam_list : list
///     Paths to input BAM files.
/// contig_end_exclusion : int, optional
///     Exclude bases at the ends of reference sequences from calculation.
///     Default is 75.
/// min_identity : float, optional
///     Exclude reads by overall identity to the reference sequences.
///     Default is 0.97.
/// trim_lower : float, optional
///     Fraction to trim from the lower tail of the coverage distribution.
///     Default is 0.0.
/// trim_upper : float, optional
///     Fraction to trim from the upper tail of the coverage distribution.
///     Default is 0.0.
/// contig_set : set, optional
///     If provided, only the coverages of the contigs within `contig_set` will
///     returned.
///     Default is None (return the coverages of all contigs).
/// threads : int, optional
///     Number of threads to use for coverage computation. Default is 1.
///
/// Returns
/// -------
/// tuple
///     A tuple whose fist element is a list of the contig names and the second
///     one is a numpy matrix of contig coverages in the input BAM files.
#[pyfunction(
    contig_end_exclusion = "75",
    min_identity = "0.97",
    trim_lower = "0.",
    trim_upper = "0.",
    contig_set = "None",
    threads = "1"
)]
fn get_coverages_from_bam(
    py: Python,
    bam_list: Vec<&str>,
    contig_end_exclusion: u64,
    min_identity: f32,
    trim_lower: f32,
    trim_upper: f32,
    contig_set: Option<HashSet<&str>>,
    threads: usize,
) -> (PyObject, PyObject) {
    let trim_upper = 1. - trim_upper;
    let min_fraction_covered_bases = 0.;
    let filter_params = FilterParameters {
        flag_filters: FlagFilter {
            include_improper_pairs: true,
            include_secondary: false,
            include_supplementary: false,
        },
        min_aligned_length_single: 0,
        min_percent_identity_single: min_identity,
        min_aligned_percent_single: 0.,
        min_aligned_length_pair: 0,
        min_percent_identity_pair: 0.,
        min_aligned_percent_pair: 0.,
    };
    let estimators = vec![CoverageEstimator::new_estimator_trimmed_mean(
        trim_lower,
        trim_upper,
        min_fraction_covered_bases,
        contig_end_exclusion,
    )];
    let taker = CoverageTakerType::new_cached_single_float_coverage_taker(estimators.len());
    let mut estimators_and_taker = EstimatorsAndTaker { estimators, taker };
    let bam_readers = generate_filtered_bam_readers_from_bam_files(
        bam_list.clone(),
        filter_params.flag_filters.clone(),
        filter_params.min_aligned_length_single,
        filter_params.min_percent_identity_single,
        filter_params.min_aligned_percent_single,
        filter_params.min_aligned_length_pair,
        filter_params.min_percent_identity_pair,
        filter_params.min_aligned_percent_pair,
    );

    // Get the headers from the first file, and verify all files have the same header.
    let mut headers: Vec<String> = match verify_same_bam_headers(&bam_list, &bam_readers) {
        // If no files: Return early with empty list
        None => {
            assert!(bam_list.is_empty());
            return default_return_value(py, 0);
        }
        // Else: Copy the headers, since they must be moved into `contig_coverage` below.
        Some((_, headers)) => headers
            .iter()
            .map(|&v| String::from_utf8_lossy(v).into_owned())
            .collect(),
    };

    // The map lets us filter out contigs not in the contig map.
    let map = contig_set.map(|names| index_map(&headers, &names));
    if let Some(ref map) = map {
        // If the empty set is passed, return early with empty list
        if map.is_empty() {
            return default_return_value(py, bam_list.len());
        }
        let mut itr = map.iter();
        headers.retain(|_| *itr.next().unwrap() != u32::MAX)
    }

    // Compute coverage
    contig_coverage(
        bam_readers,
        &mut estimators_and_taker.taker,
        &mut estimators_and_taker.estimators,
        true,
        filter_params.flag_filters,
        threads,
    );

    // Fill in matrix
    let mut matrix = Array::zeros((headers.len(), bam_list.len()));
    match &estimators_and_taker.taker {
        CoverageTakerType::CachedSingleFloatCoverageTaker {
            stoit_names: _,
            entry_names: _,
            coverages,
            current_stoit_index: _,
            current_entry_index: _,
            num_coverages: _,
        } => {
            for (col, input_bam) in coverages.iter().enumerate() {
                for coverage_entry in input_bam {
                    let mut row = coverage_entry.entry_index;
                    if let Some(ref map) = map {
                        let i = map[coverage_entry.entry_index];
                        if i == u32::MAX {
                            continue;
                        } else {
                            row = i as usize
                        }
                    }
                    matrix[[row, col]] = coverage_entry.coverage
                }
            }
        }
        _ => unreachable!(),
    }
    (
        headers.into_py(py),
        matrix
            .to_pyarray(Python::acquire_gil().python())
            .into_py(py),
    )
}

fn default_return_value(py: Python, n_files: usize) -> (Py<PyAny>, Py<PyAny>) {
    (
        Vec::<String>::new().into_py(py),
        Array::from_elem((0, n_files), 0f32)
            .to_pyarray(Python::acquire_gil().python())
            .into_py(py),
    )
}

/// Map from index of original reference to index in new reference, only
/// keeping those headers that are in names.
/// Headers not in names will map to u32::MAX.
/// Hence, header ["a", "b", "c"] and names ["c", "a"] will give
/// [0, u32::MAX, 1]
fn index_map(headers: &[String], names: &HashSet<&str>) -> Vec<u32> {
    // Make sure there are no names in names which are not in headers.
    let mut n_seen: usize = 0;
    let mut result = vec![u32::MAX; headers.len()];
    let mut new_index: u32 = 0;
    for (i, header) in headers.iter().enumerate() {
        if names.contains(&header.as_str()) {
            result[i] = new_index;
            new_index += 1;
            n_seen += 1
        }
    }
    if n_seen != names.len() {
        panic!("Some names in `contig_set` were not found in BAM headers")
    }
    result
}

fn verify_same_bam_headers<'a, 'b>(
    bam_paths: &[&'a str],
    readers: &'b [FilteredBamReader],
) -> Option<(&'a str, Vec<&'b [u8]>)> {
    bam_paths
        .iter()
        .zip(readers.iter())
        .fold(None, |state, (filename, reader)| {
            let new_names = reader.header().target_names();
            if let Some((old_filename, old_names)) = state {
                if old_names != new_names {
                    panic!(
                        "Headers of BAM file {} does not match those of BAM file {}.",
                        old_filename, filename
                    )
                }
            }
            Some((*filename, new_names))
        })
}

#[pymodule]
fn pycoverm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(is_bam_sorted))?;
    m.add_wrapped(wrap_pyfunction!(get_coverages_from_bam))?;
    Ok(())
}
