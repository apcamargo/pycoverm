# pyCoverM

This package is still in experimental stages and aims to be a simple Python interface to [CoverM](https://github.com/wwood/CoverM)'s fast coverage estimation functions. Currently, pyCoverM provides two functions: `is_bam_sorted`, which checks whether a BAM file is sorted by coordinate, and `get_coverages_from_bam`, that computes average contig coverages from sorted BAM files.

## Installation

```
pip install pycoverm
```

## Usage

```rust
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
```

```rust
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
```
