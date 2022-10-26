# CRAMINO

A tool for quick quality assessment of cram and bam files, intended for long read sequencing.

## Installation

### Download binaries from [releases](https://github.com/wdecoster/cramino/releases)

### With cargo

`cargo install cramino`

## Usage

```bash
cramino [OPTIONS] <INPUT>

ARGS:
    <INPUT>    cram or bam file to check

OPTIONS:
    -t, --threads <THREADS>    Number of parallel decompression threads to use [default: 4]
        --hist                 If histograms have to be generated
        --checksum             If a checksum has to be calculated
    -h, --help                 Print help information
    -V, --version              Print version information
```

## Example output

```text
File name       example.cram
Number of reads 14108020
Yield [Gb]      139.91
N50     17447
Median length   6743.00
Mean length     9917
Median identity 94.27
Mean identity   92.53
Path    alignment/example.cram
Creation time   09/09/2022 10:53:36
```

A 140Gbase bam file is processed in 12 minutes, using <1Gbyte of memory. Note that the identity score above is defined as the [gap-compressed identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

### Optional output

* a checksum to check if files were updated/changed or corrupted. (`--checksum`)
* an arrow file for use within [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoComp](https://github.com/wdecoster/nanocomp) (`--arrow <filename>`)
* calculating a normalised number of reads per chromosome, e.g. to determine the sex or aneuploidies (`--karyotype`)
* histograms of read lengths and read identities, as below. (`--hist`)

```text
 70.97195691947476 ..  71.97292392225151 [  122235 ]: ∎∎
 71.97292392225151 ..  72.97389092502823 [  136051 ]: ∎∎∎
 72.97389092502823 ..  73.97485792780498 [  145876 ]: ∎∎∎
 73.97485792780498 ..   74.9758249305817 [  157751 ]: ∎∎∎
  74.9758249305817 ..  75.97679193335844 [  179551 ]: ∎∎∎∎
 75.97679193335844 ..  76.97775893613516 [  171769 ]: ∎∎∎∎
 76.97775893613516 ..   77.9787259389119 [  159340 ]: ∎∎∎
  77.9787259389119 ..  78.97969294168863 [  151355 ]: ∎∎∎
 78.97969294168863 ..  79.98065994446536 [  146207 ]: ∎∎∎
 79.98065994446536 ..  80.98162694724209 [  142832 ]: ∎∎∎
 80.98162694724209 ..  81.98259395001882 [  140902 ]: ∎∎∎
 81.98259395001882 ..  82.98356095279556 [  143909 ]: ∎∎∎
 82.98356095279556 ..  83.98452795557229 [  149142 ]: ∎∎∎
 83.98452795557229 ..  84.98549495834902 [  158386 ]: ∎∎∎
 84.98549495834902 ..  85.98646196112576 [  176819 ]: ∎∎∎∎
 85.98646196112576 ..  86.98742896390249 [  199558 ]: ∎∎∎∎
 86.98742896390249 ..  87.98839596667922 [  234573 ]: ∎∎∎∎∎
 87.98839596667922 ..  88.98936296945595 [  280849 ]: ∎∎∎∎∎∎
 88.98936296945595 ..  89.99032997223267 [  348535 ]: ∎∎∎∎∎∎∎∎
 89.99032997223267 ..   90.9912969750094 [  445640 ]: ∎∎∎∎∎∎∎∎∎∎
  90.9912969750094 ..  91.99226397778614 [  583424 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎
 91.99226397778614 ..  92.99323098056287 [  776111 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 92.99323098056287 ..   93.9941979833396 [ 1051370 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  93.9941979833396 ..  94.99516498611634 [ 1414103 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 94.99516498611634 ..  95.99613198889307 [ 1833438 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 95.99613198889307 ..   96.9970989916698 [ 2084833 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
  96.9970989916698 ..  97.99806599444653 [ 1620179 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
 97.99806599444653 ..  98.99903299722327 [  416669 ]: ∎∎∎∎∎∎∎∎∎
 98.99903299722327 ..                100 [   39254 ]:
```
