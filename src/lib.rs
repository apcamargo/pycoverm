use coverm::{
    bam_generator::*, contig::*, coverage_takers::*, mosdepth_genome_coverage_estimators::*,
    FlagFilter,
};
use ndarray::Array2;
use numpy::convert::ToPyArray;
use pyo3::{prelude::*, wrap_pyfunction};
use rust_htslib::{bam, bam::Read};
use std::collections::{HashMap, HashSet};
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
/// bam_file : list
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
    // Read the first 64 bytes of the BAM header
    let header = &bam::Header::from_template(bam.header()).to_bytes()[0..64];
    let header = str::from_utf8(header).unwrap();
    if header.contains("SO:coordinate") {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// get_coverages_from_bam(bam_list, contig_end_exclusion=75, min_identity=0.97, threads=1)
/// --
///
/// Computes contig mean coverages from sorted BAM files. Trimmed means will be
/// computed if `trim_min` and/or `trim_max` are set to values greater than 0.
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
    let mut estimators_and_taker = EstimatorsAndTaker {
        estimators: estimators,
        taker: taker,
    };
    let bam_readers = generate_filtered_bam_readers_from_bam_files(
        bam_list,
        filter_params.flag_filters.clone(),
        filter_params.min_aligned_length_single,
        filter_params.min_percent_identity_single,
        filter_params.min_aligned_percent_single,
        filter_params.min_aligned_length_pair,
        filter_params.min_percent_identity_pair,
        filter_params.min_aligned_percent_pair,
    );
    contig_coverage(
        bam_readers,
        &mut estimators_and_taker.taker,
        &mut estimators_and_taker.estimators,
        true,
        filter_params.flag_filters,
        threads,
    );

    let mut contig_coverages = HashMap::new();
    match &estimators_and_taker.taker {
        CoverageTakerType::CachedSingleFloatCoverageTaker {
            stoit_names: _,
            entry_names,
            coverages,
            current_stoit_index: _,
            current_entry_index: _,
            num_coverages: _,
        } => {
            for input_bam in coverages {
                for coverage_entry in input_bam {
                    if let Some(contig_name) = &entry_names[coverage_entry.entry_index] {
                        let contig_coverage_vector =
                            contig_coverages.entry(contig_name).or_insert(vec![]);
                        contig_coverage_vector.push(coverage_entry.coverage);
                    }
                }
            }
        }
        _ => unreachable!(),
    }

    // If `contig_set` is `None`, the coverage of every contig found within the input BAM files
    // will be stored. If `contig_set` is a set containing contig names, only the coverages of those
    // contigs will be stored.

    let (filter_contigs, contig_set) = match contig_set {
        Some(_) => (true, contig_set.unwrap()),
        None => (false, HashSet::new()),
    };

    let mut coverage_vector: std::vec::Vec<f32> = Vec::new();
    let mut contig_names_vector = Vec::new();

    if filter_contigs {
        for (contig_name, contig_coverage_vector) in contig_coverages.iter() {
            if contig_set.contains(contig_name.as_str()) {
                coverage_vector.extend(contig_coverage_vector);
                contig_names_vector.push(*contig_name);
            }
        }
    } else {
        for (contig_name, contig_coverage_vector) in contig_coverages.iter() {
            coverage_vector.extend(contig_coverage_vector);
            contig_names_vector.push(*contig_name);
        }
    }

    let coverage_vector = Array2::from_shape_vec(
        (
            contig_names_vector.len(),
            coverage_vector.len() / contig_names_vector.len(),
        ),
        coverage_vector,
    )
    .unwrap()
    .to_pyarray(Python::acquire_gil().python())
    .into_py(py);
    let contig_names_vector = contig_names_vector.into_py(py);

    (contig_names_vector, coverage_vector)
}

#[pymodule]
fn pycoverm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(is_bam_sorted))?;
    m.add_wrapped(wrap_pyfunction!(get_coverages_from_bam))?;
    Ok(())
}
