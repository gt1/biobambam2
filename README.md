biobambam2
==========

This package contains some tools for processing BAM files including

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
