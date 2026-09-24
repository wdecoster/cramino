# Changelog

## 2.0.0

### Breaking changes to the output

If you parse cramino output in a workflow, please check the following:

* **`N50` and `N75` now use the full read length.** They are computed over the primary alignments, including soft and hard clipped bases. Before, they were computed over the aligned part (without clipped bases) of the primary and supplementary alignments. As a result, the N50 of a bam file is now the same as that of the ubam with the same reads, and is higher than before for the same file (by 13% for the test data in this repository).
* **`Median length` and `Mean length` now also use the full read length**, for the same reads as the N50 and N75.
* **New fields `N50 aligned`, `N75 aligned`, `Median length aligned` and `Mean length aligned`** (`n50_aligned`, `n75_aligned`, `median_length_aligned` and `mean_length_aligned` in tsv and json), with these statistics as they were computed before: over the aligned part of the primary and supplementary alignments. In the text and tsv output, these follow `Mean length`, so all columns after it (identity, phasing, ...) have moved 4 positions to the right, and anything parsing tsv columns by position has to be updated.
* **`Mean coverage` is `NA`** (`null` in json) when the header has no reference sequences, such as for a ubam, instead of `inf`. If no reads pass the filters, it is 0 (or `NA` without reference sequences).
* **`% from total alignments` is renamed to `% from total reads`** in the text output (the tsv/json key `percent_from_total` is unchanged), and now counts every read once. Before, supplementary alignments were counted in the total but not in the reads used, so a file in which all reads are mapped showed e.g. 83% instead of 100%.
* **Phaseblocks are now correct (`--phased`), which changes all phaseblock statistics.** A phaseblock now spans from the leftmost start to the rightmost end of the alignments (primary and supplementary) with the same phaseset on a chromosome. Before, a block ended at the end of the alignment with the highest start, rather than at the rightmost end, and a block was split each time reads of another, overlapping phaseset were encountered along the genome. For the test data in this repository, this gives 8 instead of 11 phaseblocks, with a median length of 503851 instead of 246860.
* **`Fraction reads phased` only counts reads** (primary alignments). Before, phased supplementary alignments were also counted, so the fraction could exceed 1.
* **The karyotype (`--karyotype`) output changed.**
  * It is consistent between formats: the json `normalized_count` is now the same median-normalized value as in the text output (before, it was the number of alignments per bp), and the tsv output now has a `karyotype_<chromosome>` column per chromosome.
  * It counts reads (primary alignments) instead of alignments, so supplementary alignments are no longer counted.
  * All chromosomes in the header are listed, in the order of the header, including those without reads (with 0). Before, those were left out, and the chromosomes were in alphabetical order in the text output and in random order in json. The median used for normalization is over the chromosomes with reads, and is now the average of the two middle values for an even number of chromosomes.
* **With `--ubam`, the identity is reported as estimated identity.** It is estimated from the base qualities, and now named accordingly in all formats: `Median estimated identity` (text, before `Median est. identity`), `median_estimated_identity` (tsv, before `median_est_identity`), and `estimated_identity_stats.median_estimated_identity` (json, before `identity_stats.median_identity` with `is_estimated: true`). The same goes for the mean and modal identity. Json `identity_stats` no longer has the `is_estimated` field.
* **Reads without an identity are left out of the identity statistics and histogram**, rather than counted as an identity of 0: with `--ubam` the reads without base qualities, and otherwise mapped reads without aligned bases. In the arrow output, their identity is null.
* **A missing NM tag is an error.** Without a de tag, the NM tag is required to calculate the identity of aligned reads. cramino now exits with an error naming the read, rather than panicking.
* **`--hist` and `--hist-count` require `=` for a file**, e.g. `--hist=histograms.txt`. Before, `cramino --hist sample.bam` took `sample.bam` as the histogram file and read the input from stdin. Also, with `--format json`, histograms without a file are only included in the json, which was invalid json before, as the histograms were printed after it. With `--format tsv`, a file is required.
* **The modal identity is deterministic.** If multiple identities are the most common, the highest is reported. Before, one was chosen at random.
* **`--min-read-len` now filters on the aligned length**, i.e. the same length that is reported, and keeps alignments of exactly that length (`>=` instead of `>`). Before, it filtered on the length of the stored sequence, so alignments shorter than the minimum could still be reported.

Other than described above, `Yield`, `Mean coverage` and the read length histogram are unchanged, and still computed over the aligned part of the primary and supplementary alignments.

### Performance

* Lengths are stored as 32-bit instead of 128-bit integers. Together with the 4 bytes per read for the new read length statistics, this reduces memory use by about 150 MB for a file with 17 million alignments of 14 million reads.

### Fixes

* Records without a stored sequence (SEQ `*`) but with soft clipped bases resulted in an overflow, reporting a nonsensical N50 of 340282366920938463463374607431768211436. Lengths are now taken from the CIGAR.
* Environment variables for remote (https/s3) access are now set before any threads are started.
* Remote input is recognized by its scheme (`s3://`, `https://`, `http://`, `ftp://`). Before, a local file with a name starting with `s3` was treated as a URL (and one starting with `http`, `ftp` or `s3` had no creation time), and `http://` and `ftp://` input could not be opened. Support for `s3://` is now also compiled in.
* A mapped read without aligned bases (CIGAR `*`), or with a NaN de tag, no longer causes a panic, and a de tag stored as a double is accepted.
* Piping the output into e.g. `head` no longer causes a panic.
* Errors opening or reading the input, e.g. a file that does not exist or a CRAM file of which the reference is not found, and a negative or non-integer PS tag are reported without a panic.
* `--reference` is also used when the input is not a file ending in `.cram`, e.g. a CRAM on stdin.
* The gap-compressed identity underflowed, reporting e.g. a median identity of -429067552.05, when the NM tag was lower than the number of inserted and deleted bases. Such an inconsistent NM is now taken as no mismatches.
* `--threads 0` gives an error message instead of a panic.
* The N50 (and N75, and phaseblock N50) is now the length for which *at least* half of the bases are in reads of that length or longer, rather than *more than* half. This only makes a difference when the cumulative length is exactly half of the total.

### Other

* Removed unused code.
