# iCount-Mini

This is a fork of iCount maintained by members of [Jernej Ule's](http://ulelab.info) group, focussing on the peak calling features of iCount.

Run commands using:
`iCount-Mini <command>`

**Note on small differences of terminology between iCount-Mini and iCount**
+ In iCount-Mini, sigxls = iCount peaks and iCount-Mini peaks = iCount clusters. This is to bring the terminology more in line with the rest of the field.
+ In iCount-Mini RNA-maps have been renamed to 'metagene', to distinguish these plots which include only CLIP data from other RNA-maps which group crosslinks into categories dependent on orthogonal data, such as alternatively spliced exons.

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
