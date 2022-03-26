# mark-nonconverted-reads
This program examines a methylation sequencing (EM-seq or bisulfite) BAM or SAM stream (or file) for reads that contain multiple nonconverted Cs in non-CpG context. Particularly in [FFPE samples](https://doi.org/10.1373/clinchem.2014.223040), such filtering can markedly reduce false-methylation detection. This treatment can be applied to BAMs produced using any alignment method, but it is not recommend for use in organisms that have significant methylation outside of CpG context (e.g. most plants).

mncr marks problematic reads with a tag (XX:Z:UC) and optionally sets the vendor failed bit to prevent these reads from being used in downstream methylation calling. It also prints a count of nonconverted reads per contig to stderr.

By default, the program will read from stdin, tag reads that have 3 or more nonconverted Cs, flag them, and output to STDOUT. These can all be changed with command line arguments.

*Note*: This does not take read pairing into account. Each read will be tagged or not independently of its mate.

Optional arguments:

| Argument | Details |
| ------- | ------- |
|--reference \<Filename> | Reference fasta file [default = searches the bam file for a bwameth command]|
|--bam \<Input file> | File must end in .bam or .sam [default = stdin]|
|--out \<Output name> | Name of the output sam file [default = stdout]|
|--c_count \<int> |Minimum number of nonconverted Cs on a read to consider it nonconverted [default = 3]|
|--flag_reads |Set the 'Failing platform / vendor quality check' flag [default = don't set it]|
