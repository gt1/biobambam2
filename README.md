biobambam2
==========

This package contains some tools for processing BAM files including

 - bamsormadup: parallel sorting and duplicate marking
 - bamcollate2: reads BAM and writes BAM reordered such that alignment
   or collated by query name
 - bammarkduplicates: reads BAM and writes BAM with duplicate alignments
   marked using the BAM flags field
 - bammaskflags: reads BAM and writes BAM while masking (removing)
   bits from the flags column
 - bamrecompress: reads BAM and writes BAM with a defined compression
   setting. This tool is capable of multi-threading.
 - bamsort: reads BAM and writes BAM resorted by coordinates or query
   name
 - bamtofastq: reads BAM and writes FastQ; output can be collated or
   uncollated by query name

A short list of options is available for each program by calling it
with the -h parameter, e.g.

	bamsort -h

Source
------

The biobambam2 source code is hosted on github:

	git@github.com:gt1/biobambam2.git

Compilation of biobambam2
-------------------------

biobambam2 needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then biobambam2 can be compiled and
installed in ${HOME}/biobambam2 using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/biobambam2
	- make install

Using bamsormapdup
------------------

bamsormadup is a new tool in biobambam2. In has two modes of operation. 
If SO=coordinate (as it is by default) then it expects a name collated (all reads for one name appear consecutively) input file.
It sorts this file by coordinate and marks duplicate reads and read pairs and outputs the sorted file in bam, cram or sam format.
If SO=queryname then it expects an input file in any order and sorts this file by query name.
In both cases the program can produce checksums at read set level (like bamseqchksum) and at file level (like md5sum). If the output file is
sorted by coordinate and in bam format, then program can also produce a bam index on the fly. Most stages of bamsormadup are parallelised. The number
of threads used can be set with the threads option (e.g. threads=4). If no such option is given then the number of logical CPUs (cores) of the machine is used
as the number of threads. Parallel sorting of alignment files is I/O heavy, so a fast I/O system is crucial. We recommend to store all data
(input, tempory and output) on solid state storage (SSD).
