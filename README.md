# CRAMINO

A tool for quick quality assessment of cram and bam files, intended for long read sequencing.

## Installation

Preferably, for most users, download a ready-to-use binary for your system to add directory on your $PATH from the [releases](https://github.com/wdecoster/cramino/releases).  
You may have to change the file permissions to execute it with `chmod +x cramino`

Alternatively, use conda to install  
`conda install -c bioconda cramino`

Or for Rust developers, build cramino with cargo:  
`cargo install cramino`

## Usage

```text
cramino [OPTIONS] [INPUT]

Arguments:
  [INPUT]  cram or bam file to check [default: -]

Options:
  -t, --threads <THREADS>            Number of parallel decompression threads to use [default: 4]
      --reference <REFERENCE>        reference for decompressing cram
  -m, --min-read-len <MIN_READ_LEN>  Minimal aligned length of a read to be considered (full length for unmapped reads) [default: 0]
      --hist[=<FILE>]                If histograms have to be generated (optionally specify output file as --hist=FILE)
      --arrow <ARROW>                Write data to an arrow format file
      --karyotype                    Provide normalized number of reads per chromosome
      --phased                       Calculate metrics for phased reads
      --spliced                      Provide metrics for spliced data
      --ubam                         Also include unmapped reads, and estimate the identity from base qualities (disables --karyotype, --phased and --spliced)
      --format <FORMAT>              Output format (text, json, or tsv) [default: text]
      --scaled                       Scale histogram bins by total basepairs in each bin (not just read count)
      --hist-count[=<FILE>]          Output read length histogram bin counts in TSV format (optionally specify output file as --hist-count=FILE), cannot be combined with --hist
  -h, --help                         Print help
  -V, --version                      Print version
```

## Example output

```text
File name               small-test-phased.bam
Number of alignments    7416
% from total reads      100.00
Number of reads         6176
Yield [Gb]              0.08
Mean coverage           0.02
Yield [Gb] (>25kb)      0.03
N50                     21885
N75                     12356
Median length           8283.50
Mean length             12426.93
N50 aligned             19294
N75 aligned             10697
Median length aligned   6767.50
Mean length aligned     10367.90

Median identity         98.92
Mean identity           98.22
Modal identity          99.1

Path                    test-data/small-test-phased.bam
Creation time           18/11/2025 14:42:51
```

A 140Gbase bam file is processed in 12 minutes, using <1Gbyte of memory. Note that the identity score above is defined as the [gap-compressed identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity). The `--ubam` flag will provide metrics for all reads in the file, regardless of whether they are aligned or not. With `--ubam`, the identity is instead estimated from the base qualities, and reported as `estimated identity`; reads without base qualities are left out. `--ubam` disables `--karyotype`, `--phased` and `--spliced`.

The input can also be a remote file (`https://`, `http://`, `ftp://` or `s3://`), or `-` for stdin (the default).

### Which records are counted

Secondary alignments are never counted, and unmapped reads only with `--ubam`.

* `Number of reads` counts the reads, i.e. the primary alignments (and unmapped reads with `--ubam`). `Number of alignments` also includes the supplementary alignments.
* `% from total reads` is the percentage of all reads in the file that are used for this report. This depends on the `--min-read-len` and `--ubam` settings; without both of those, it is the % of reads that are mapped.
* `N50`, `N75`, `Median length` and `Mean length` are computed over the full length of each read, including soft and hard clipped bases. For the same reads, a bam file and its ubam therefore have the same values.
* The `aligned` statistics, `Yield`, `Mean coverage` and the read length histogram are computed over the aligned part of the primary and supplementary alignments, without clipped bases. A read that is split in several alignments therefore counts as several shorter alignments, and the `N50 aligned` is typically lower than the `N50`. For unmapped reads (with `--ubam`), the full read length is used, so for a ubam the `aligned` statistics are the same as the others.
* The identity statistics and histogram are computed over the primary and supplementary alignments.
* `Mean coverage` is `NA` (`null` in json) if the header has no reference sequences, e.g. for a ubam.
* The identity of aligned reads is calculated from the `de` tag, or otherwise from the `NM` tag and the CIGAR. cramino exits with an error if a read has neither, and `NM` tags can be added with `samtools calmd`.
* `--min-read-len` applies to the aligned length of each alignment (to the full length for unmapped reads), and keeps alignments of at least that length.

The N50 is the length for which half of all bases are in reads of that length or longer. Unlike the median or mean length, it is weighted by the number of bases, and therefore reflects the long reads that make up most of the data.

Other records, such as duplicates or QC-failed reads, can be filtered with samtools, and piped to cramino, which reads from stdin by default. Piping uncompressed BAM (`-u`) avoids spending time on compression in between:

```bash
samtools view -u -F 0x400 input.bam | cramino
```

### Optional output

* an arrow file for use within [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoComp](https://github.com/wdecoster/nanocomp) (`--arrow <filename>`)
* calculating a normalised number of reads per chromosome, e.g. to determine the sex or aneuploidies (`--karyotype`)
* information about the phase blocks. (`--phased`)
* the median and mean number of exons per read, and the fraction of unspliced reads. (`--spliced`)
* histograms of read lengths and read identities, as below. (`--hist`, or `--hist=FILE` to write them to a file). With `--phased`, also a histogram of phase block lengths, and with `--spliced` of the number of exons. With `--scaled`, read length and Phred accuracy histograms are basepair-weighted. Please let me know if the histograms look inappropriately scaled for your data.
* read length histogram bin counts in TSV format (`--hist-count`, or `--hist-count=FILE`), which cannot be combined with `--hist`. With `--scaled`, the TSV values are basepair totals instead of read counts.

With `--format tsv`, histograms can only be written to a file. With `--format json`, histograms are only written separately if a file is given. When `--hist` or `--hist-count` is set, JSON output includes the read length and identity histograms under `histograms.read_length` and `histograms.q_score`, each with its `step`, `max_value` and `bins`. Each bin has a `start`, `end` (absent for the overflow bin), `count` and `bases`, independent of `--scaled`. The phase block and exon histograms are not included in JSON.

The histograms below are for `test-data/small-test-phased.bam`:

```text
# Histogram for read lengths:
     0-2000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  2000-4000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  4000-6000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  6000-8000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 8000-10000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
10000-12000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
12000-14000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
14000-16000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
16000-18000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
18000-20000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
20000-22000 ∎∎∎∎∎∎∎∎∎∎∎∎∎
22000-24000 ∎∎∎∎∎∎∎∎∎∎∎∎∎
24000-26000 ∎∎∎∎∎∎∎∎∎∎
26000-28000 ∎∎∎∎∎∎∎∎
28000-30000 ∎∎∎∎∎∎∎∎
30000-32000 ∎∎∎∎∎∎
32000-34000 ∎∎∎∎∎∎
34000-36000 ∎∎∎
36000-38000 ∎∎∎∎
38000-40000 ∎∎∎
40000-42000 ∎∎
42000-44000 ∎
44000-46000 ∎
46000-48000 ∎
48000-50000 ∎
50000-52000 
52000-54000 
54000-56000 
56000-58000 
58000-60000 
     60000+ 


# Histogram for Phred-scaled accuracies:
  Q0-1 
  Q1-2 
  Q2-3 
  Q3-4 
  Q4-5 
  Q5-6 
  Q6-7 ∎
  Q7-8 ∎
  Q8-9 ∎
 Q9-10 ∎
Q10-11 ∎∎∎∎∎∎
Q11-12 ∎∎∎∎∎∎∎∎∎
Q12-13 ∎∎∎∎∎∎∎∎∎∎
Q13-14 ∎∎∎∎∎∎∎∎∎∎∎∎∎
Q14-15 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q15-16 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q16-17 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q17-18 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q18-19 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q19-20 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q20-21 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q21-22 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q22-23 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q23-24 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q24-25 ∎∎∎∎∎∎∎∎
Q25-26 ∎∎∎∎
Q26-27 ∎∎
Q27-28 ∎
Q28-29 
Q29-30 
Q30-31 
Q31-32 
Q32-33 
Q33-34 
Q34-35 
Q35-36 
Q36-37 
Q37-38 
Q38-39 
Q39-40 
  Q40+ ∎∎∎∎∎∎∎
```

Reproducible histogram output for `test-data/small-test-phased.bam` is available in `docs/histogram-example.txt` (unscaled) and `docs/histogram-example-scaled.txt` (scaled).

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
