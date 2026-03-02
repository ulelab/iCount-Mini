# iCount-Mini

This is a fork of iCount maintained by members of [Jernej Ule's](http://ulelab.info) group, focussing on the peak calling features of iCount.

Run commands using:
`iCount-Mini <command>`

## What's new in v4.0.0

### Updated type hierarchy

The default type hierarchy has been updated to split non-coding RNA exons into short ncRNAs (highest priority) and lncRNA exons (below UTRs):

`ncRNA > CDS > UTR3 > UTR5 > lncRNA > intron > intergenic`

This ensures short non-coding RNAs (miRNA, snRNA, snoRNA, rRNA, tRNA, etc.) take priority over protein-coding annotations at overlapping positions, while long non-coding RNA exons are ranked below UTRs. lncRNA introns remain in the standard intron category.

Additionally, `SUBTYPE_GROUPS` has been updated to match current ENSEMBL annotations (e.g. `lncRNA` biotype alongside the older `lincRNA`, `vault_RNA`, `protein_coding_CDS_not_defined`, `artifact`).

### Runner-up type attribute

Every region in `regions.gtf.gz` now carries a `runner_up` attribute recording the next-best type at overlapping positions. For example, if an ncRNA gene overlaps a CDS from another gene, the region type is `ncRNA` and `runner_up` is `CDS`. When there is no overlap, `runner_up` is `NA`. This is useful for understanding what lies underneath the winning annotation.

### tRNA annotation in segment

`iCount segment` now accepts an optional `--trna_annotation` parameter pointing to a BED file with tRNA gene coordinates (e.g. from GtRNAdb / tRNAscan-SE). These entries are loaded as ncRNA with biotype=tRNA and merged into the segmentation alongside the main GTF annotation.

Example:
```bash
iCount-Mini segment annotation.gtf segmentation.gtf genome.fai \
    --trna_annotation hg38-tRNAs.bed
```

### Summary output format

All summary TSV files include a `Length` column showing the total genomic length (bp) of contributing annotation regions for each row. Overlay and isotype summaries include a `TOTAL` row at the top showing the aggregate cDNA count and percentage for the entire category.

Cross-tabulated labels use a colon separator (e.g. `Ala:intron`, `LINE:intron`) for easy downstream parsing.

### tRNA isotype summary

When tRNA regions are present in the annotation, `iCount summary` automatically produces a tRNA isotype-level summary (`summary_tRNA_isotype.tsv`). This cross-tabulates tRNA isotypes (parsed from `gene_name`, e.g. `tRNA-Ala-AGC-1-1` yields isotype `Ala`) with the runner-up type from the annotation (CDS, intron, intergenic, etc.), formatted as `Isotype:RunnerUp`.

### Overlay annotation summaries

`iCount summary` now accepts `--overlay_annotations` to produce cross-tabulated summaries for additional annotation layers such as repetitive elements (TEs) or DNA cis-regulatory elements (CREs). Each overlay annotation is intersected with cross-link sites and combined with the region type to produce a `summary_{name}.tsv` file. Strand-specific intersection is auto-detected (used when overlay features have strand, omitted for unstranded annotations like CREs).

Example:
```bash
iCount-Mini summary regions.gtf.gz sites.bed out_dir/ \
    --overlay_annotations "TE.gtf:TE:gene_id;CRE.gtf:CRE:gene_id"
```

The format for overlay annotations is `gtf_path:name:group_by_attribute`, semicolon-separated. The `group_by_attribute` specifies which GTF attribute to use as the group label (e.g. `gene_id`, `family_id`). Output rows are cross-tabulated as `Group:Type` (e.g. `LINE:intron`, `SINE:UTR3`, `Alu:intergenic`).

### Performance improvements

Segmentation performance has been significantly improved through the following changes:

| Change | Impact |
|--------|--------|
| Removed redundant full-GTF parse that was only used for progress counting | ~50% reduction in GTF parsing time |
| Biotype classification now uses a dict lookup instead of linear scan | O(1) vs O(n) per biotype |
| Gene/transcript ID tracking uses sets instead of lists | O(1) vs O(n) membership checks (~60K genes, ~250K transcripts) |
| Chromosome filtering uses a set | O(1) vs O(n) per GTF line (~3.5M lines for human) |
| Pre-compiled regex for attribute parsing in the regions hot loop | Avoids repeated pattern compilation |

### Chromosome name auto-detection

When a tRNA BED file uses UCSC-style chromosome names (`chr1`, `chrM`, ...) but the annotation uses ENSEMBL-style names (`1`, `MT`, ...), the mismatch is auto-detected and names are converted using a shipped mapping table (`hg38_ucsc_to_ensembl.txt`). Entries on chromosomes with no ENSEMBL equivalent are silently skipped.

---

**Note on small differences of terminology between iCount-Mini and iCount**
+ In iCount-Mini, sigxls = iCount peaks and iCount-Mini peaks = iCount clusters. This is to bring the terminology more in line with the rest of the field.
+ In iCount-Mini RNA-maps have been renamed to 'metagene', to distinguish these plots which include only CLIP data from other RNA-maps which group crosslinks into categories dependent on orthogonal data, such as alternatively spliced exons.

**Note on peak calling with iCount-Mini**

Note that to call peaks with iCount-Mini you must run three commands:
1. Firstly you will need to run `iCount-Mini segment` to segment your gtf file into genomic regions.
2. You need to run `iCount-Mini sigxls` to call statistically significant crosslinks.
3. You need to run `iCount-Mini peaks` to merge your significant crosslinks into broader peak regions.

# iCount: protein-RNA interaction analysis

iCount is a Python module and associated command-line interface (CLI), which provides all the commands needed to process iCLIP data on protein-RNA interactions and generate:
 
+ demultiplexed and adapter-trimmed FASTQ files
+ BAM files with mapped iCLIP reads
+ identified protein-RNA cross-linked sites, saved to BED files
+ statistically significant cross-linked sites, saved to BED files
+ peaks of significant cross-linked sites, saved to BED files
+ grouping of individual replicate experiments
+ metagene generation showing the positional distribution of cross-linked sites relative to genomic landmarks
+ kmer enrichment analysis

You may start with the [tutorial](http://icount.readthedocs.io/en/latest/tutorial.html) or dive into the 
[documentation](http://icount.readthedocs.io/en/latest/index.html).

## iCount-Mini Authors

iCount-Mini is maintained by members of [Jernej Ule's](http://ulelab.info) group.

## iCount Authors

iCount is developed and supported by [Tomaž Curk](http://curk.info) from the [Bioinformatics Laboratory](http://biolab.si) at the [University of Ljubljana](http://www.uni-lj.si), [Faculty of Computer and Information Science](http://www.fri.uni-lj.si) and in collaboration with the laboratory of [Jernej Ule](http://ulelab.info).

The development started in late 2008 when Tomaž Curk and Gregor Rot wrote a first prototype of iCount.
In mid-2016, [Jure Zmrzlikar](https://github.com/JureZmrzlikar) from [Genialis](http://www.genialis.com) helped refactoring and improving the code, which is now available here.

## Development

To install a development version of iCount-Mini, use this command.
It's recommended to do this within a Python virtual environment.

```bash
pip install --upgrade -r requirements-rtd.txt -e .
```
