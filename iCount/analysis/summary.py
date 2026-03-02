""".. Line to protect from pydocstyle D205, D400.

Cross-link site summary
-----------------------

Report count of cross-link events in each region type.
"""
import logging
import math
import re
import os
import tempfile

from pybedtools import BedTool

import iCount
from iCount.genomes.segment import summary_templates, sort_types_subtypes, TEMPLATE_TYPE, TEMPLATE_SUBTYPE, \
    TEMPLATE_GENE, SUMMARY_TYPE, SUMMARY_SUBTYPE, SUMMARY_GENE, SUMMARY_TRNA_ISOTYPE

LOGGER = logging.getLogger(__name__)


def _parse_isotype(gene_name):
    """Parse tRNA isotype from gene_name (e.g. 'tRNA-Ala-AGC-1-1' -> 'Ala')."""
    if not gene_name:
        return gene_name
    parts = gene_name.split('-')
    if len(parts) >= 2 and parts[0] in ('tRNA', 'tRX'):
        return parts[1]
    return gene_name


def _parse_overlay_annotations(overlay_annotations):
    """
    Parse overlay_annotations string into list of (gtf_path, name, group_by) tuples.

    Parameters
    ----------
    overlay_annotations : str
        Semicolon-separated specs in format ``gtf_path:name:group_by_attribute``.

    Returns
    -------
    list
        List of (gtf_path, name, group_by) tuples.

    """
    specs = []
    for spec in overlay_annotations.split(';'):
        parts = spec.strip().split(':')
        if len(parts) != 3:
            raise ValueError(
                'Invalid overlay_annotations spec "{}". Expected format: gtf_path:name:group_by_attribute'.format(spec))
        gtf_path, name, group_by = parts
        specs.append((gtf_path.strip(), name.strip(), group_by.strip()))
    return specs


def _build_site_type_lookup(sites, annotation):
    """Build a lookup mapping (chrom, start, stop, strand) -> region type."""
    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    overlaps = BedTool(sites).intersect(
        BedTool(annotation), sorted=True, s=True, wb=True, nonamecheck=True).saveas()
    # pylint: enable=too-many-function-args,unexpected-keyword-arg
    lookup = {}
    for seg in overlaps:
        key = (seg.chrom, seg.start, seg.stop, seg.strand)
        lookup[key] = seg[8]
    return lookup


def isotype_summary(overlaps, sum_cdna, out_dir):
    """
    Produce tRNA isotype-level summary cross-tabulated with the runner-up type.

    Scans the main site-vs-annotation overlaps for regions with tRNA biotype,
    parses the isotype from gene_name, and uses the runner_up attribute to
    show what other type underlies the tRNA (ncRNA) region.

    Auto-detected: only produces output if tRNA regions are found.

    Parameters
    ----------
    overlaps : pybedtools.BedTool
        Pre-computed site-vs-annotation intersection (with wb=True).
    sum_cdna : float
        Total cDNA count across all sites.
    out_dir : str
        Output directory.

    """
    from iCount.genomes.segment import SUBTYPE_GROUPS
    trna_biotypes = set(SUBTYPE_GROUPS.get('tRNA', []))

    isotype_counter = {}
    label_regions = {}
    all_regions = set()
    for segment in overlaps:
        biotype = re.match(r'.*biotype "(.*?)";', segment[-1])
        biotype = biotype.group(1) if biotype else ''
        biotype_list = biotype.split(',')

        if not any(bt in trna_biotypes for bt in biotype_list):
            continue

        score = float(segment.score)
        gene_name = re.match(r'.*gene_name "(.*?)";', segment[-1])
        gene_name = gene_name.group(1) if gene_name else ''
        isotype = _parse_isotype(gene_name)

        runner_up = re.match(r'.*runner_up "(.*?)";', segment[-1])
        runner_up = runner_up.group(1) if runner_up else 'NA'

        label = '{}:{}'.format(isotype, runner_up)
        isotype_counter[label] = isotype_counter.get(label, 0) + score

        region_key = (segment[6], segment[9], segment[10], segment[12])
        label_regions.setdefault(label, set()).add(region_key)
        all_regions.add(region_key)

    if not isotype_counter:
        LOGGER.info('No tRNA cross-links found, skipping isotype summary.')
        return

    total_trna_cdna = sum(isotype_counter.values())
    label_lengths = {lbl: sum(int(e) - int(s) for _, s, e, _ in regs) for lbl, regs in label_regions.items()}
    total_trna_length = sum(int(e) - int(s) for _, s, e, _ in all_regions)

    LOGGER.info('Writing tRNA isotype report...')
    with open(os.path.join(out_dir, SUMMARY_TRNA_ISOTYPE), 'wt') as out:
        header = ['Isotype:RunnerUp', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        out.write('\t'.join(map(str, [
            'TOTAL', total_trna_length, math.floor(total_trna_cdna),
            total_trna_cdna / sum_cdna * 100])) + '\n')
        for label, cdna in sorted(isotype_counter.items()):
            line = [label, label_lengths.get(label, 0), math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')


def overlay_summary(sites, overlay_gtf, regions, name, group_by, out_dir):
    """
    Produce a cross-tabulated summary of overlay annotation groups x region type.

    For each cross-link site overlapping an overlay feature, extract the group
    label from the specified attribute and look up the region type.

    Parameters
    ----------
    sites : str
        Cross-links file (BED6 format).
    overlay_gtf : str
        Overlay annotation GTF (e.g. TE or CRE annotations).
    regions : str
        Regions annotation file (regions.gtf.gz).
    name : str
        Name for the overlay (used in output filename).
    group_by : str
        GTF attribute to use as group label (e.g. 'gene_id', 'family_id').
    out_dir : str
        Output directory.

    """
    LOGGER.info('Building overlay summary for %s (group_by=%s)...', name, group_by)

    type_lookup = _build_site_type_lookup(sites, regions)

    # Use strand-specific intersection only if the overlay has strand info
    overlay_bt = BedTool(overlay_gtf)
    first_entry = next(iter(overlay_bt))
    stranded = first_entry.strand in ('+', '-')
    if stranded:
        LOGGER.info('Overlay %s is stranded, using strand-specific intersection.', name)
    else:
        LOGGER.info('Overlay %s is unstranded, using strand-agnostic intersection.', name)

    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    intersect_kwargs = dict(sorted=True, wb=True, nonamecheck=True)
    if stranded:
        intersect_kwargs['s'] = True
    overlay_overlaps = BedTool(sites).intersect(overlay_bt, **intersect_kwargs).saveas()
    # pylint: enable=too-many-function-args,unexpected-keyword-arg

    try:
        overlay_overlaps[0]
    except (IndexError, TypeError):
        LOGGER.warning('No intersections found for overlay %s. Skipping.', name)
        return

    sum_cdna = 0
    for seg in BedTool(sites):
        sum_cdna += float(seg.score)

    overlay_counter = {}
    label_features = {}
    all_features = set()
    for seg in overlay_overlaps:
        score = float(seg.score)
        key = (seg.chrom, seg.start, seg.stop, seg.strand)
        region_type = type_lookup.get(key, 'intergenic')

        group_label = re.match(r'.*{} "(.*?)";'.format(re.escape(group_by)), seg[-1])
        group_label = group_label.group(1) if group_label else 'unknown'

        label = '{}:{}'.format(group_label, region_type)
        overlay_counter[label] = overlay_counter.get(label, 0) + score

        feat_key = (seg[6], seg[9], seg[10], seg[12])
        label_features.setdefault(label, set()).add(feat_key)
        all_features.add(feat_key)

    total_overlay_cdna = sum(overlay_counter.values())
    label_lengths = {lbl: sum(int(e) - int(s) for _, s, e, _ in feats) for lbl, feats in label_features.items()}
    total_overlay_length = sum(int(e) - int(s) for _, s, e, _ in all_features)

    out_file = 'summary_{}.tsv'.format(name)
    LOGGER.info('Writing overlay report: %s', out_file)
    with open(os.path.join(out_dir, out_file), 'wt') as out:
        header = ['Group:Type', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        out.write('\t'.join(map(str, [
            'TOTAL', total_overlay_length, math.floor(total_overlay_cdna),
            total_overlay_cdna / sum_cdna * 100])) + '\n')
        for label, cdna in sorted(overlay_counter.items()):
            line = [label, label_lengths.get(label, 0), math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')


def summary_reports(annotation, sites, out_dir, templates_dir=None, overlay_annotations=None):
    """
    Make summary reports for a cross-link file.

    Parameters
    ----------
    annotation : str
        Annotation file (GTF format). It is recommended to use genome-level segmentation (e.g. regions.gtf.gz), that
        is produced by ``iCount segment`` command.
    sites : str
        Croslinks file (BED6 format). Should be sorted by coordinate.
    out_dir : str
        Output directory.
    templates_dir : str
        Directory containing templates for summary calculation. Made by ``iCount segment`` command. If this argument
        is not provided, summary templates are made on the fly.
    overlay_annotations : str
        Overlay annotation GTFs for cross-tabulated summaries. Semicolon-separated specs in format
        gtf_path:name:group_by_attribute (e.g. TE.gtf:TE:gene_id;CRE.gtf:CRE:gene_id).

    Returns
    -------
    iCount.Metrics
        iCount Metrics object.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()

    if templates_dir is None:
        templates_dir = tempfile.mkdtemp()
        summary_templates(annotation, templates_dir)

    LOGGER.info('Calculating intersection between cross-link and annotation...')
    # pylint: disable=too-many-function-args,unexpected-keyword-arg
    overlaps = BedTool(sites).intersect(
        BedTool(annotation),
        sorted=True,  # invokes memory efficient algorithm for large files
        s=True,  # only report hits in B that overlap A on the same strand
        wb=True,  # write the original entry in B for each overlap
        nonamecheck=True,  # Do not print warnings about name inconsistency to stdout
    ).saveas()
    # pylint: enable=too-many-function-args,unexpected-keyword-arg
    try:
        overlaps[0]  # will raise Error if overlaps is empty:
    except (IndexError, TypeError):
        raise ValueError('No intersections found. This may be caused by different naming of chromosomes in annotation'
                         'and cross-links file (example: "chr1" vs. "1")')

    type_counter, subtype_counter, gene_counter = {}, {}, {}
    LOGGER.info('Extracting summary data from intersection...')
    for segment in overlaps:
        score = float(segment.score)

        type_ = segment[8]
        type_counter[type_] = type_counter.get(type_, 0) + score

        biotype = re.match(r'.*biotype "(.*?)";', segment[-1])
        biotype = biotype.group(1) if biotype else ''
        biotypes = biotype.split(',')
        for biotype in biotypes:
            sbtyp = iCount.genomes.segment.make_subtype(type_, biotype)
            subtype_counter[sbtyp] = subtype_counter.get(sbtyp, 0) + score / len(biotypes)

        gene_id = re.match(r'.*gene_id "(.*?)";', segment[-1])
        gene_id = gene_id.group(1) if gene_id else None
        gene_counter[gene_id] = gene_counter.get(gene_id, 0) + score

    sum_cdna = 0
    for seg in BedTool(sites):
        sum_cdna += float(seg.score)

    def parse_template(template_file):
        """Parse template file."""
        template = {}
        with open(template_file, 'rt') as ifile:
            for line in ifile:
                line = line.strip().split('\t')
                template[line[0]] = line[1:]
        return template

    LOGGER.info('Writing type report...')
    type_template = parse_template(os.path.join(templates_dir, TEMPLATE_TYPE))
    with open(os.path.join(out_dir, SUMMARY_TYPE), 'wt') as out:
        header = ['Type', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for type_, cdna in sorted(type_counter.items(), key=lambda x: sort_types_subtypes(x[0])):
            line = [type_, type_template.get(type_, [-1])[0], math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Writing subtype report...')
    subtype_template = parse_template(os.path.join(templates_dir, TEMPLATE_SUBTYPE))
    with open(os.path.join(out_dir, SUMMARY_SUBTYPE), 'wt') as out:
        header = ['Subtype', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for stype, cdna in sorted(subtype_counter.items(), key=lambda x: sort_types_subtypes(x[0])):
            line = [stype, subtype_template.get(stype, [-1])[0], math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    LOGGER.info('Writing gene report...')
    gene_template = parse_template(os.path.join(templates_dir, TEMPLATE_GENE))
    with open(os.path.join(out_dir, SUMMARY_GENE), 'wt') as out:
        header = ['Gene name (Gene ID)', 'Length', 'cDNA #', 'cDNA %']
        out.write('\t'.join(header) + '\n')
        for gene_id, cdna in sorted(gene_counter.items()):
            gene_name, length = gene_template.get(gene_id, ['', -1])
            if gene_id == '.':
                gene_name = 'intergenic'
            line = ['{} ({})'.format(gene_name, gene_id), length, math.floor(cdna), cdna / sum_cdna * 100]
            out.write('\t'.join(map(str, line)) + '\n')

    # tRNA isotype summary (auto-detected from annotation biotypes)
    isotype_summary(overlaps, sum_cdna, out_dir)

    # Overlay summaries (TE, CRE, etc.)
    if overlay_annotations:
        overlay_specs = _parse_overlay_annotations(overlay_annotations)
        for gtf_path, name, group_by in overlay_specs:
            overlay_summary(sites, gtf_path, annotation, name, group_by, out_dir)

    LOGGER.info('Done.')
    return metrics
