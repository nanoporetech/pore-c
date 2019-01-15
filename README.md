# Pore-C tools


## Getting Started

### Creating a development environment

```
conda env create
source activate poreC
pip install -e .
pore_c --help
```

### Running tests

```
pytest tests
```


## The Pore-C tool suite

- generate-fragments
- cluster-reads
- map-to-fragments
- flatten-multiway
- hic-to-hicpro

### generate-fragments

This is redundant with (https://github.com/aidenlab/juicer/blob/master/misc/generate_site_positions.py) in the juicer universe and generates an identical file. The included version of the tools supports degenerate sequences and sequence groups coded as REGEXes and prints to stdout. 

Signature:
` Usage: pore_c generate-fragments [OPTIONS] REFERENCE_FASTA RESTRICTION_PATTERN`
`                                 OUTFILE`
`Try "pore_c generate-fragments --help" for help.`

Examples:

`pore_c generate-fragments reference.fa “AAGCTT” ref.HindIII.hicREF		#a single site`

`pore_c generate-fragments “(GAATTC |GCGGCCGC)” ref.mix.hicREF	#two sites: ecoRI + NotI`

`pore_c generate-fragments reference.fa “GCCNNNNNGGC” ref.BglI.hicREF		#BglI (with degeneracy)`

`pore_c generate-fragments reference.fa “RAATY” ref.ApoI.hicREF		#ApoI (with degeneracy)`

The resulting reference file consists of one space-separated line per chromosome. The first item of each line is the chromosome name, the subsequent items are restriction site positions. The final item is the length of the chromosome.



### clusterReads.py

This tools is applied to the bam file produced during alignment of poreC reads to a reference genome. It removes spurious mappings, such as is seen when part or all of a restriction fragment contains an Alu, causing that interval to map to thousands of copies of that Alu around the genome. The step generates two .bam files, a file name <yourfile>.split.XX.bam and <yourfile>.split.XX.removed.bam. The first file contains a most complete tilepath of mappings across each of the pore-C reads such that no two mappings overlaps by more than XX bp. The removed file contains, as the name suggests, all the mappings that were deemed uninformative.

Example:
`clusterReads.py 20180821_FAJ15224_DD.bwa.namesort.bam 20`

### mapToFrags.py

This tool assigns each mapping in the filtered bam file to a restriction fragment and generates a multiway contact .poreC file, which is similar to the text version of the juicer .hic file format, but expanded so that there are more than two contacts per entry. 

Example:
`mapToFrags.py --alignments 20180821_FAJ15224_DD.bwa.namesort.split.20.bam --fragments ref.EcoRI.hicREF > 20180821_FAJ15224_DD.bwa.namesort.split.20.poreC`

The file specification for a .poreC file is derived from the text version of the .hic file described by the Aiden lab. It contains the same tuples of fields (strand, chr,fragPosition, fragment ID), but because multiway contacts contain a notion of fragment order (whereas Hi-C entries only contain an association without orientation), the entries are not sorted as they are in a .hic file. The mapq values are appended to the end of the line. The fields for a 3-way contact would be as follows:

`readID str1 chr1 fragPos1 fragID1 str2 chr2 fragPos2 fragID2 str3 chr3 fragPos3 fragID3 mapq1 mapq2 mapq3`

and a 4 way contact:

`readID str1 chr1 fragPos1 fragID1 str2 chr2 fragPos2 fragID2 str3 chr3 fragPos3 fragID3 str4 chr4 fragPos4 fragID4 mapq1 mapq2 mapq3 mapq4`

If the space separated length of the line is L, the number of contacts in a multiway contact is therefore (L-1) / 5

**Nota bene**: If you would like to combine multiple runs into a single data set, this can be done with a simple `cat` operation of all the .poreC files into a single .poreC file. I think that concatenating .hic.txt files can also be done (after flattening, for example), however this requires multiple rounds of sorting in order to retain compatibility with `juicer_tools.jar pre`. 


### poreCflatten.py

This tool takes in a multiway contact .poreC file and flattens that data to a single dimensionality, e.g. all pairwise interactions or all three-way interactions. For a concatemer of 4 mappings A B C D, the pairwise contacts produced are AB, AC, AD, BC, BD, CD. The sort option sorts the fragments by their fragment ID in each entry. Hi-C platforms often require fragments to be in sorted order, however multiway contacts lose information potentially when the fragments are sorted, so this --sort option is included to be able to choose whether the output is sorted.

Example:
`poreCflatten.py --poreC 20180821_FAJ15224_DD.bwa.namesort.split.20.poreC --setSize 2 --sort > 20180821_FAJ15224_DD.bwa.namesort.split.20.hic.txt`

The output of the above command can be turned into a binary .hic file using the `juicer_tools.jar pre` operation.



### hicToHiCPro.py

This tool converts a text .hic file (produced by using poreCflatten.py on a .poreC file) into a format that can be loaded into HiCPro for analysis.



