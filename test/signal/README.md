Test files for the signal track readers (tools/test-signal.js).

- From the IGV test set (github.com/igvteam/igv, test/data; MIT license): NA12878.SLX.egfr.sam.tdf,
  dm3_var_sample.wig.*.tdf, hg18_var_sample.wig.*.tdf, wholeGenomeTest.tdf with its bedGraph,
  recombRate.igv.txt.tdf (9 tracks), MNNBBJG-example4.entries.tdf.
- From the igv.js test set (github.com/igvteam/igv.js, test/data; MIT license): BRAF_C.bam.tdf,
  gstt1_sample.bam.tdf, fixedStep.bw with fixedStep.wig, bigEndian.bw, uncompressed.bw,
  chromTreeByteOrder.bw, nozooms.bw, bigtools.bigWig.
- Written by tools/make-signal-tests.py (synthetic values): af_sample.tdf, dense01_sample.tdf.
  With "big" it also writes large files into big/ (not kept in the repository).
