<h1 align="center">pyCoverM</h1>

pyCoverM is a Python library that provides a simple interface to [CoverM](https://github.com/wwood/CoverM) for fast coverage estimation.

## Installation

pyCoverM is available through PyPI or Conda.

### PyPI installation

```sh
pip install pycoverm
```

### Conda installation

The Conda package can be installed though [Pixi](https://pixi.sh/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/)/[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html).

```sh
# Pixi
pixi init pycoverm_project
cd pycoverm_project
pixi project channel add bioconda
pixi add pycoverm

# Mamba (just replace 'mamba' with 'conda' if you have Conda installed)
mamba create -n pycoverm_env -c conda-forge -c bioconda pycoverm
mamba activate pycoverm_env
```

### Quick start

pyCoverM provides two functions:
1. `is_bam_sorted` - Verifies if a BAM file is sorted by coordinate.
2. `get_coverages_from_bam` - Computes average contig coverages from sorted BAM files.


```py
>>> import pycoverm
>>> TEST_BAM = "tests/test_data.bam"
>>> pycoverm.is_bam_sorted(TEST_BAM)
```

    True

```py
>>> coverages = pycoverm.get_coverages_from_bam([TEST_BAM])
>>> coverages[0]
```

    ['contig_7847997', 'contig_11746202', 'contig_9129108', …, 'contig_2917594']

```py
>>> coverages[1]
```

    array([[0.        ],
           [0.526652  ],
           [0.08541025],
           …           ,
           [0.00907206]], dtype=float32)


> [!NOTE]
> If multiple BAM files are provided, the resulting NumPy array will contain one column for each BAM file, with each column corresponding to the coverage values from a specific BAM file.

## API

```rs
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

```rs
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
///     A list of paths to input BAM files.
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