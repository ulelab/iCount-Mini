# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from pybedtools import create_interval_from_list

from iCount.analysis import rnamaps
from iCount.files import remove_extension
from iCount.tests.utils import get_temp_file_name, make_file_from_list, make_list_from_file


class TestGetGeneName(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        seg_level0 = create_interval_from_list(['1', '.', 'CDS', '1', '2', '.', '+', '.', 'gene_name "G0";'])
        seg_level1 = create_interval_from_list(['1', '.', 'UTR3', '1', '2', '.', '+', '.', 'gene_name "G1";'])
        seg_level4 = create_interval_from_list(['1', '.', 'intron', '1', '2', '.', '+', '.', 'gene_name "G2";'])

        self.assertEqual('G0', rnamaps.get_gene_name(seg_level0, seg_level1))
        self.assertEqual('G0', rnamaps.get_gene_name(seg_level1, seg_level0))
        self.assertEqual('G1', rnamaps.get_gene_name(seg_level1, seg_level4))
        with self.assertRaises(ValueError):
            self.assertEqual('B', rnamaps.get_gene_name(seg_level0, seg_level0))


class TestMakeLandmarksFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '+', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [
            ['chr1', '200', '201', 'A', '.', '+'],
        ])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [
            ['chr1', '199', '200', 'A', '.', '-'],
        ])

    def test_limits_upstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '400', '.', '+', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '350', '.', '-', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [])

    def test_limits_downstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '350', '.', '+', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '-', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', 'gene_name "A";'],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [])


class TestComputeDistances(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        self.landmarks = make_file_from_list([
            ['chr1', '10', '11', 'G1', '.', '+'],
            ['chr1', '20', '21', 'G1', '.', '+'],
            ['chr1', '20', '21', 'G2', '.', '-'],
            ['chr2', '10', '11', 'G3', '.', '+'],
        ])

    def test_basic(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 3)
        self.assertEqual(distances, {
            'chr1__+__10__G1': {2: 3},
        })

    def test_scores_sum(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr1', '12', '13', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1__+__10__G1': {2: 4},
        })

    def test_strands_not_mixed(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr1', '12', '13', '.', '1', '-'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1__+__10__G1': {2: 3},
            'chr1__-__20__G2': {8: 1},
        })

    def test_chroms_not_mixed(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr2', '12', '13', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1__+__10__G1': {2: 3},
            'chr2__+__10__G3': {2: 1},
        })

    def test_size_limit(self):
        xlinks = make_file_from_list([
            ['chr1', '22', '23', '.', '3', '+'],
            ['chr2', '222', '223', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1__+__20__G1': {2: 3},
        })

    def test_no_landmark(self):
        """Landmark is missing on this chromosome / stramd."""
        xlinks = make_file_from_list([
            ['chrX', '22', '23', '.', '3', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, 'exon-intron')
        self.assertEqual(total_cdna, 3)
        self.assertEqual(distances, {})


class TestMakeFullResultsFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.fname = get_temp_file_name()

    def test_make_full_results_file(self):
        distances = {
            'chr1__+__10__G1': {2: 3},
            'chr2__+__20__G2': {-3: 3, 1: 5},
        }
        rnamaps.make_results_raw_file(distances, self.fname, total_cdna=11, maptype='exon-intron')

        result = make_list_from_file(self.fname, fields_separator='\t')
        self.assertEqual(result[0], ['total_cdna:11'])
        self.assertEqual(result[1], ['.'] + list(map(str, range(-50, 151))))
        self.assertEqual(result[2][0], 'chr1__+__10__G1')
        self.assertEqual(result[2][45:55], ['0', '0', '0', '0', '0', '0', '0', '0', '3', '0'])
        self.assertEqual(result[3][45:55], ['0', '0', '0', '3', '0', '0', '0', '5', '0', '0'])


class TestRun(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.outdir = get_temp_file_name()

    def test_run(self):
        regions = make_file_from_list(sort=True, data=[
            ['chr1', '.', 'intergenic', '1', '210', '.', '+', '.', 'gene_name "None";'],
            ['chr1', '.', 'UTR5', '211', '270', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'CDS', '271', '330', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intron', '331', '490', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'CDS', '491', '550', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'UTR3', '551', '760', '.', '+', '.', 'gene_name "A";'],
            ['chr1', '.', 'intergenic', '761', '1100', '.', '+', '.', 'gene_name "None";'],

            ['chr1', '.', 'intergenic', '1', '300', '.', '-', '.', 'gene_name "None";'],
            ['chr1', '.', 'ncRNA', '301', '500', '.', '-', '.', 'gene_name "B";'],
            ['chr1', '.', 'intron', '501', '600', '.', '-', '.', 'gene_name "B";'],
            ['chr1', '.', 'ncRNA', '601', '750', '.', '-', '.', 'gene_name "B";'],
            ['chr1', '.', 'intergenic', '751', '1000', '.', '-', '.', 'gene_name "None";'],
        ])

        sites = make_file_from_list([
            ['chr1', '220', '221', '.', '1', '+'],
            ['chr1', '350', '351', '.', '1', '+'],
            ['chr1', '350', '351', '.', '1', '-'],
            ['chr1', '550', '551', '.', '1', '+'],
            ['chr1', '740', '741', '.', '1', '+'],
            ['chr1', '750', '751', '.', '1', '-'],
        ])

        rnamaps.run(sites, regions, outdir=self.outdir)

        self.assertTrue(os.path.isdir(self.outdir))

        sites_name = remove_extension(sites, ['.bed', '.bed.gz'])
        for maptype in rnamaps.MAP_TYPES:
            basename = os.path.join(self.outdir, '{}_{}'.format(sites_name, maptype))
            # for extension in ['.tsv', '.png', '_plot_data.txt']:
            for extension in ['.tsv', '.png']:
                fname = basename + extension
                self.assertTrue(os.path.isfile(fname))
                self.assertGreater(os.path.getsize(fname), 1)


if __name__ == '__main__':
    unittest.main()
