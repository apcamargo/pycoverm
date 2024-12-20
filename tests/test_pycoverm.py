import numpy as np
import pycoverm

TEST_BAM = "tests/test_data.bam"


def test_is_bam_sorted():
    result = pycoverm.is_bam_sorted(TEST_BAM)
    assert result is True


def test_get_coverages_from_single_bam():
    result = pycoverm.get_coverages_from_bam([TEST_BAM])
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], list)
    assert isinstance(result[1], np.ndarray)
    assert result[1].shape == (547, 1)
    np.testing.assert_almost_equal(result[1].sum(), 30.173173904418945)


def test_get_coverages_from_multiple_bams():
    result = pycoverm.get_coverages_from_bam([TEST_BAM, TEST_BAM, TEST_BAM])
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], list)
    assert isinstance(result[1], np.ndarray)
    assert result[1].shape == (547, 3)
    np.testing.assert_almost_equal(result[1].sum(), 90.51951599121094)
