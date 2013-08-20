r240 | 2013-07-09 11:30:03 +0800 (Tue, 09 Jul 2013)

Fix a bug in reading files in 'map' step. This bug might lead to seg fault.

------------------------------------------------------------------------
r239 | 2013-06-26 09:41:39 +0800 (Wed, 26 Jun 2013)

1) Fix the bug of reading fasta file in map step. This bug was introduced when
   fixing a bug of reading fastq file in r238.

------------------------------------------------------------------------
r224 - r238 | 2013-06-13

1) Fix a serious bug in 'map' step of version r223. This bug can lead
   to incorrect pairing of PE reads in LIB of even order, e.g., the 
   2nd LIB, the 4th LIB and so on...And these affected LIBs may not 
   contribute to the construction of scaffold.
2) Merge 'standPregraph' and 'sparsePregraph'. Now, there are only two 
   executable programs: SOAPdenovo-63mer and SOAPdenovo-127mer. User 
   can choose to use 'pregraph' for standard Kmer graph or 
   'sparse_pregraph' for sparse Kmer graph.
3) Add an option for debug version compilation. User can use 
   'make debug=1' to obtain programs for debug.
4) Fix a bug in sorting edges in 'contig' step.
5) Fix a bug in reading files when using multi-kmer. Now the 
   'max_read_length' will change according to the LIB being red.

------------------------------------------------------------------------
r223 | 2012-12-28 10:11:43 +0800 (Fri, 28 Dec 2012) 

Fix the problem that parameter k doesn't work when k is larger than 63
for 127mer version.

------------------------------------------------------------------------
r222 | 2012-12-21 14:45:49 +0800 (Fri, 21 Dec 2012) 

1) Change some codes so that program can handle reads longer than 5000.
2) Add a new perl script which can seperate singletons from scaffolds in 
   *.scafSeq file.

------------------------------------------------------------------------
r221 | 2012-12-07 14:27:02 +0800 (Fri, 07 Dec 2012)

Fix a bug in reading files which might cause zombie process.

------------------------------------------------------------------------
r220 | 2012-11-26 10:09:45 +0800 (Mon, 26 Nov 2012)

Fix bug in aio that the buffer was not enough for fq for long reads.

------------------------------------------------------------------------
r219 | 2012-11-08 12:58:45 +0800 (Thu, 08 Nov 2012)

Fix a bug that using -r 1 will casuse the infomation loss of  MaxReadLen
and MinReadLen in *.preGraphBasic file in pregraph_sparse module.

------------------------------------------------------------------------
r218 | 2012-11-08 11:04:54 +0800 (Thu, 08 Nov 2012)

Output palindrome sequence only once now instead of twice before in
pregraph_sparse module.

------------------------------------------------------------------------
r217 | 2012-11-01 13:09:50 +0800 (Thu, 01 Nov 2012)  

Fix bug in scaffolding which may lead to scaffold consisting of none 
or only one contig.

------------------------------------------------------------------------
r216 | 2012-10-31 14:53:29 +0800 (Wed, 31 Oct 2012)

Fix a bug of 'pregraph-sparse' which may lead to segmentation fault in 
'contig' step if option -R is set and there are reads longer than 100bp.

------------------------------------------------------------------------
r215 | 2012-10-16 18:53:28 +0800 (Tue, 16 Oct 2012)

Fix a bug of aio which happens rarely in 'pregraph' step when there are 
reads shorter than Kmer.

------------------------------------------------------------------------
r214 | 2012-10-08 15:58:09 +0800 (Mon, 08 Oct 2012)

Modify usage description of '-V'.

------------------------------------------------------------------------
r213 | 2012-09-29 09:24:32 +0800 (Sat, 29 Sep 2012)

Fix a bug which might happen in 'contig' step if the 'pregraph-sparse' is
used to replace the regular 'pregraph'.

