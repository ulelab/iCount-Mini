""".. Line to protect from pydocstyle D205, D400.

RNA maps
--------

Perform RNA-maps analysis.
"""
import csv
import logging
import os

import pandas as pd
from pybedtools import BedTool, create_interval_from_list

import iCount  # pylint: disable=wrong-import-position

LOGGER = logging.getLogger(__name__)

MAP_TYPES = {
    'exon-intron': {
        'upstream-type': ['CDS', 'UTR5'],
        'upstream-size-limit': 50,
        'downstream-type': ['intron'],
        'downstream-size-limit': 150,
        'downstream-plot-width': 100,
    },
    'intron-exon': {
        'upstream-type': ['intron'],
        'upstream-size-limit': 150,
        'upstream-plot-width': 100,
        'downstream-type': ['CDS', 'UTR3'],
        'downstream-size-limit': 50,
    },
    'gene-start': {
        'upstream-type': ['intergenic'],
        'upstream-size-limit': 200,
        'downstream-type': ['UTR5', 'CDS'],
        'downstream-size-limit': 50,
    },
    'gene-end': {
        'upstream-type': ['CDS', 'UTR3'],
        'upstream-size-limit': 200,
        'downstream-type': ['intergenic'],
        'downstream-size-limit': 200,
    },
    'translation-start': {
        'upstream-type': ['UTR5'],
        'upstream-size-limit': 50,
        'downstream-type': ['CDS'],
        'downstream-size-limit': 50,
    },
    'translation-end': {
        'upstream-type': ['CDS'],
        'upstream-size-limit': 50,
        'downstream-type': ['UTR3'],
        'downstream-size-limit': 200,
    },
    'noncoding-gene-start': {
        'upstream-type': ['intergenic'],
        'upstream-size-limit': 200,
        'downstream-type': ['ncRNA'],
        'downstream-size-limit': 100,
    },
    'noncoding-gene-end': {
        'upstream-type': ['ncRNA'],
        'upstream-size-limit': 100,
        'downstream-type': ['intergenic'],
        'downstream-size-limit': 200,
    },
}


def get_gene_name(seg1, seg2):
    """Determine gene name to use in landmarks file."""
    level1 = iCount.genomes.segment.TYPE_HIERARCHY.index(seg1[2])
    level2 = iCount.genomes.segment.TYPE_HIERARCHY.index(seg2[2])

    if level1 < level2:
        return seg1.attrs.get('gene_name', '.')
    elif level2 < level1:
        return seg2.attrs.get('gene_name', '.')
    else:
        raise ValueError('Not possible')


def make_landmarks_file(regions, maptype):
    """Create landmarks file from regions file."""
    up_types = MAP_TYPES[maptype]['upstream-type']
    up_limit = MAP_TYPES[maptype].get('upstream-size-limit', 0)
    down_types = MAP_TYPES[maptype]['downstream-type']
    down_limit = MAP_TYPES[maptype].get('downstream-size-limit', 0)

    intervals = []
    for strand in ['+', '-']:
        seg1 = None
        for seg2 in BedTool(regions).filter(lambda x: x.strand == strand):  # pylint: disable=cell-var-from-loop
            if seg1 and seg1.chrom == seg2.chrom:

                # pylint: disable=unsubscriptable-object
                if (strand == '+' and
                        seg1[2] in up_types and
                        seg2[2] in down_types and
                        len(seg1) > up_limit and
                        len(seg2) > down_limit):
                    # BED6 iz zero based and does not include stop nuclotide
                    # Landmark position is the first nucleotide in the downstream feature
                    # + strand: ---up---||---down---
                    #           01234567890123456789
                    # - strand: ---down---||---up---
                    gene_name = get_gene_name(seg1, seg2)
                    intervals.append([seg1.chrom, seg1.stop, seg1.stop + 1, gene_name, '.', strand])

                elif (strand == '-' and
                      seg1[2] in down_types and
                      seg2[2] in up_types and
                      len(seg1) > down_limit and
                      len(seg2) > up_limit):
                    gene_name = get_gene_name(seg1, seg2)
                    intervals.append([seg1.chrom, seg1.stop - 1, seg1.stop, gene_name, '.', strand])
                # pylint: enable=unsubscriptable-object

            seg1 = seg2

    landmarks = BedTool(create_interval_from_list(list_) for list_ in intervals).sort().saveas()
    return landmarks.fn


def compute_distances(landmarks, sites, maptype):
    """Compute distances between each xlink and it's closest landmark."""
    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    closest = BedTool(sites).closest(  # pylint: disable=assignment-from-no-return
        landmarks,
        s=True,
        t='first',
        D='a',
        nonamecheck=True,
    )
    # pylint: enable=too-many-function-args,unexpected-keyword-arg

    distances = {}
    total_cdna = 0
    for mseg in closest:
        # "mseg" means merged segment, since it consists of 3 parts:
        # sites segment (BED6) + landmark segment (BED6) + distance
        if len(mseg.fields) != 13:
            LOGGER.warning('Segment length shoud be 13, not %s. Segment fields: %s', str(len(mseg)), str(mseg.fields))

        score = int(mseg[4])
        chrom = mseg[6]
        pos = mseg[7]
        gene_name = mseg[9]
        strand = mseg[11]
        distance = -int(mseg[-1])
        total_cdna += score

        if not -MAP_TYPES[maptype]['upstream-size-limit'] <= distance <= MAP_TYPES[maptype]['downstream-size-limit']:
            continue
        if chrom == '.' or strand not in ['+', '-'] or '-' in pos:
            continue

        # loc = Landmark "ID" = landmark exact coordinates
        loc = '{}__{}__{}__{}'.format(chrom, strand, pos, gene_name)

        distances.setdefault(loc, {})[distance] = distances.get(loc, {}).get(distance, 0) + score

    return distances, total_cdna


def make_results_raw_file(distances, fname, total_cdna, maptype):
    """Write distances data to file."""
    up_limit = -MAP_TYPES[maptype]['upstream-size-limit']
    down_limit = MAP_TYPES[maptype]['downstream-size-limit']
    header = list(range(up_limit, down_limit + 1))

    with open(fname, 'wt') as handle:
        outfile = csv.writer(handle, delimiter='\t')
        outfile.writerow(['total_cdna:{}'.format(total_cdna)])

        outfile.writerow(['.'] + header)
        for loc, positions in sorted(distances.items()):
            outfile.writerow([loc] + [positions.get(pos, 0) for pos in header])


def make_results_summarised_file(outdir, fname):
    """Write "plot data" to file."""
    data = {}
    up_limit = -max([mtyp['upstream-size-limit'] for mtyp in MAP_TYPES.values()])
    down_limit = max([mtyp['downstream-size-limit'] for mtyp in MAP_TYPES.values()])
    header = list(range(up_limit, down_limit + 1))

    for basename in [file_ for file_ in os.listdir(outdir) if file_.endswith('.tsv')]:
        results_file = os.path.join(outdir, basename)
        maptype = iCount.plotting.plot_rnamap.guess_maptype(results_file)
        plot_data, _ = iCount.plotting.plot_rnamap.parse_results(results_file)

        # Make date of the same size, impute with 0 score on locations with no data.
        data[maptype] = [plot_data.get(pos, 0) for pos in header]
        assert len(header) == len(data[maptype])

    dframe = pd.DataFrame.from_dict(data, orient='index', columns=header)
    dframe.to_csv(fname, sep='\t')


def run(sites,
        regions,
        outdir=None,
        top_n=100,
        smoothing=1,
        nbins=50,
        binsize=None,
        colormap='Greys',
        imgfmt='png',
        ):
    """
    Compute distribution of cross-links relative to genomic landmarks.

    Parameters
    ----------
    sites : str
        Croslinks file (BED6 format). Should be sorted by coordinate.
    regions : str
        Regions file (regions.gtf.gz) that is produced by ``iCount segment``.
    outdir : str
        Output directory.
    top_n : int
        Plot heatmap for top_n best covered landmarks.
    smoothing : int
        Smoothing half-window. Average smoothing is used.
    nbins : int
        Number of bins. Either nbins or binsize can be defined, but not both.
    binsize : int
        Bin size. Either nbins or binsize can be defined, but not both.
    colormap : str
        Colormap to use. Any matplotlib colormap can be used.
    imgfmt : str
        OUtput image format.

    Returns
    -------
    None

    """
    iCount.logger.log_inputs(LOGGER)

    sites_name = iCount.files.remove_extension(sites, ['.bed', '.bed.gz'])
    if outdir is None:
        outdir = os.path.join(os.path.abspath(os.getcwd()), 'rnamaps_{}'.format(sites_name))
        LOGGER.info('Output directory not given, creating one at %s', outdir)
    os.makedirs(outdir, exist_ok=True)

    for maptype in MAP_TYPES:
        LOGGER.info('Creating landmarks for %s', maptype)
        landmarks = make_landmarks_file(regions, maptype)
        if not BedTool(landmarks).head(n=1, as_string=True):
            LOGGER.info('No landmarks for %s', maptype)
            continue

        LOGGER.info('Processing data for %s', maptype)
        distances, total_cdna = compute_distances(landmarks, sites, maptype)
        if not distances:
            LOGGER.warning('No distances for %s', maptype)
            continue

        LOGGER.info('Writing results to file for %s', maptype)
        # File with full set of results:
        results_raw_file = os.path.join(outdir, '{}_{}.tsv'.format(sites_name, maptype))
        make_results_raw_file(distances, results_raw_file, total_cdna, maptype)

        iCount.plotting.plot_combined.plot_combined(
            results_raw_file,
            os.path.join(outdir, '{}_{}.{}'.format(sites_name, maptype, imgfmt)),
            up_limit=MAP_TYPES[maptype].get('upstream-plot-width', MAP_TYPES[maptype]['upstream-size-limit']),
            down_limit=MAP_TYPES[maptype].get('downstream-plot-width', MAP_TYPES[maptype]['downstream-size-limit']),
            top_n=top_n,
            smoothing=smoothing,
            nbins=nbins,
            binsize=binsize,
            colormap=colormap,
        )

    # Single file with only RNA-maps distibution plot data.
    results_summarised_file = os.path.join(outdir, '{}_plot_data.tsv'.format(sites_name))
    make_results_summarised_file(outdir, results_summarised_file)

    LOGGER.info('Done.')
