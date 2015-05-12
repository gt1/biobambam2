biobambam2

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
(input, temporary and output) on solid state storage (SSD).

### bamsormadup and cram output

CRAM output with bamsormadup requires libmaus to be built with support for io_lib version 1.3.11 (or newer) or a sufficiently recent svn revision of version 1.3.10.
The binaries provided on github have support for CRAM writing. The program will look for reference sequences in the following places while encoding cram:

- The directory stored in the environment variable REF_CACHE (if any). The md5 hash value is used as a file name within this directory. For instance if REF_CACHE is
  set to ${HOME}/ref and a reference sequence (as stated in a SQ header line) has an md5 hash 01234567890123456789012345678901 then the program would look up the file
  ${HOME}/ref/01234567890123456789012345678901 . If REF_CACHE contains finite length string references, then parts of the hash will be inserted before adding
  the rest of the hash in the end . If for instance REF_CACHE is set to ${HOME}/ref/%2s/%2s/%s then the program would look for the sequence above in the file
  ${HOME}/ref/01/23/4567890123456789012345678901 . REF_CACHE designates a read and write cache . The program will try to produce reference sequences not previously 
  stored in this directory if it is given a FastA or gzipped FastA file as a reference (either via the reference command line key or via the UR field of the corresponding
  sequence line).
- The list of directories and URL prefixes stored in the environment variable REF_PATH (if any). Multiple paths are separated by the colon symbol ':'.
  A colon sign in a path can be escaped by duplicating it, e.g. the URL http://www.ebi.ac.uk/ena/cram/md5/ would be escaped as http:://www.ebi.ac.uk/ena/cram/md5/ .
  The locations given in this list are considered as read only. URLs must be specified using the URL= prefix, e.g. URL=http://www.ebi.ac.uk/ena/cram/md5/ .
- A FastA or gzipped FastA file given as the UR parameter in the header of the input file. This file will be scanned to obtain the reference sequence. All newly found
  sequences will be stored in the REF_CACHE directory if the respective environment variable is set.
- A FastA file given via the reference key on the command line.

If the program cannot find a required reference sequence for encoding CRAM in any of these locations, then it fails. Note that if a reference sequence is only present
in a FastA file then the REF_CACHE environment variable must be set or CRAM encoding will fail in the current version.

