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
cramino [OPTIONS] <INPUT>

Arguments:
  [INPUT]  cram or bam file to check [default: -]

Options:
  -t, --threads <THREADS>            Number of parallel decompression threads to use [default: 4]
      --reference <REFERENCE>        reference for decompressing cram
  -m, --min-read-len <MIN_READ_LEN>  Minimal length of read to be considered [default: 0]
      --hist [<FILE>]                If histograms have to be generated (optionally specify output file)
      --scaled                       Scale histogram bins by total basepairs in each bin (not just read count)
      --hist-count [<FILE>]          Output histogram bin counts in TSV format (optionally specify output file)
      --arrow <ARROW>                Write data to an arrow format file
      --karyotype                    Provide normalized number of reads per chromosome
      --phased                       Calculate metrics for phased reads
      --spliced                      Provide metrics for spliced data
      --ubam                         Provide metrics for unaligned reads
      --format <FORMAT>              Output format (text, json, or tsv) [default: text]
  -h, --help                         Print help
  -V, --version                      Print version
```

## Example output

```text
File name       example.cram
Number of reads 14108020
% from total reads  83.45
Yield [Gb]      139.91
N50     17447
Median length   6743.00
Mean length     9917
Median identity 94.27
Mean identity   92.53
Path    alignment/example.cram
Creation time   09/09/2022 10:53:36
```

A 140Gbase bam file is processed in 12 minutes, using <1Gbyte of memory. Note that the identity score above is defined as the [gap-compressed identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity). The `--ubam` flag will provide metrics for all reads in the file, regardless of whether they are aligned or not.
The `% from total reads` output field contains the percentage of reads used for this report, depending on the `--min-read-len` and `--ubam` settings. Without both of those, this indicates the % of reads that are mapped, primary or supplementary.

### Optional output

* a checksum to check if files were updated/changed or corrupted. (`--checksum`)
* an arrow file for use within [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoComp](https://github.com/wdecoster/nanocomp) (`--arrow <filename>`)
* calculating a normalised number of reads per chromosome, e.g. to determine the sex or aneuploidies (`--karyotype`)
* information about the phase blocks. (`--phased`)
* information about number of splice sites. (`--spliced`)
* histograms of read lengths and read identities, as below. (`--hist`). With `--phased`, also a histogram of phase block lengths. With `--scaled`, read length and Phred accuracy histograms are basepair-weighted. Please let me know if the histograms look inappropriately scaled for your data.
* histogram bin counts in TSV format (`--hist-count`). With `--scaled`, the TSV values are basepair totals instead of read counts.

When `--hist` or `--hist-count` is set, JSON output includes histogram bins under `histograms.read_length` and `histograms.q_score`. Each bin includes `start`, `end` (or `null` for overflow), `count`, and `bases`.

```text
# Histogram for read lengths:
     0-2000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  2000-4000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  4000-6000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  6000-8000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 8000-10000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
10000-12000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
12000-14000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
14000-16000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
16000-18000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
18000-20000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
20000-22000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
22000-24000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
24000-26000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
26000-28000 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
28000-30000 ∎∎∎∎∎∎∎∎∎∎∎∎
30000-32000 ∎∎∎∎∎∎∎∎∎
32000-34000 ∎∎∎∎∎∎
34000-36000 ∎∎∎∎
36000-38000 ∎∎
38000-40000 ∎
40000-42000 ∎
42000-44000 ∎
44000-46000 
46000-48000 
48000-50000 
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
  Q5-6 ∎∎∎
  Q6-7 ∎∎∎∎∎∎∎∎∎∎∎∎
  Q7-8 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  Q8-9 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 Q9-10 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q10-11 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q11-12 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q12-13 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q13-14 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q14-15 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q15-16 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q16-17 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q17-18 ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
Q18-19 ∎∎∎∎
Q19-20 ∎
Q20-21 
Q21-22 
Q22-23 
Q23-24 
Q24-25 
Q25-26 
Q26-27 
Q27-28 
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
  Q40+ 
```

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
