# expHTS
========

Python application for "Experimental High Throughput Sequencing"

High Throughput Sequencing is often conducted within an "Experimental" framework, meaning multiple samples are sequenced within a single experiment and should be processed in an identicle manner. expHTS facilities this by processing multiple samples within a single frame work and generating output so that sample quality can be easily compared across experiment samples. The first sub-application associated with expHTS is sample preprocessing. The pipeline processes RawData in the following manner:

1. Screen for contaminants (minimally PhiX)
	Using a mapping based approach
2. Deduplicate Reads 
	Using Super-Deduper
3. Quality Trimming (polyA/T) trimming
	Using Scickle
4. Overlap paired-end reads
	Using Flash2
5. Modify read names, filter for too-short, QA
6. Produce multi-sample reports

