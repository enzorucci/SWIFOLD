# SWIFOLD
Smith-Waterman Acceleration on Intel’s FPGA with OpenCL for Long DNA Sequences

## Description
SWIFOLD is a software to accelerate Smith-Waterman alignment of long DNA sequences. 

## Performance

We have executed SWIFOLD in different FPGA-based platforms. Some of the performance results are detailed below. 

**Experimental platform**: A server with two Intel Xeon CPU E5-2670 8-core 2.60GHz, 48 GB main memory and an Intel Arria 10 GX with 2 GByte DDR3. Intel’s ICC compiler (version 17.0.1.132). Intel FPGA OpenCL SDK 16.0

| Input sequence| Size | Target sequence  | Size | Matrix cells   | Time     | GCUPS  |
| :----------:  | :--: | :----------:     | :--: | :----:         | :-----:  | :----: |
| AF133821.1    | 10K  | AY352275.1       | 10K  | 100K           | 0.002s   | 49.81  |
| NC_000898     | 162K | NC 007605        | 172K | 28M            | 0.83s    | 122.94 |
| BA000035.2    | 3M   | BX927147.1       | 3M   | 9G             | 2m11.4s  | 131.45 |
| NC 005027.1   | 7M   | NC 003997.3      | 5M   | 35G            | 3m27.2s  | 131.98 |
| NT 033779.4   | 23M  | NT 037436.3      | 25M | 575G            | 1h11m7s  | 132.33 |
| NC 000020.11  | 65M  | NC 006487.4      | 67M  | 4.4T           | 8h59m37s | 132.43 |


**Experimental platform**: A server with two Intel Xeon CPU E5-2670 8-core 2.60GHz, 32 GB main memory and an Altera Stratix V GSD5 Half-Length PCIe Board with Dual DDR3 (two banks of 4 GByte DDR3). Intel’s ICC compiler (version 16.0.3). Quartus II DKE V12.0 2 with OpenCL SDK v14.0.


| Input sequence| Size | Target sequence  | Size | Matrix cells   | Time     | GCUPS  |
| :----------:  | :--: | :----------:     | :--: | :----:         | :-----:  | :----: |
| AF133821.1    | 10K  | AY352275.1       | 10K | 100K            | 0.007s | 14.47  |
| NC_000898     | 162K | NC 007605        | 172K | 28M            | 0.83s | 33.45 |
| BA000035.2    | 3M   | BX927147.1       | 3M | 9G               | 4m36.8s | 37.32 |
| NC 005027.1   | 7M   | NC 003997.3      | 5M | 35G              | 16m34.5s | 37.56 |
| NT 033779.4   | 23M  | NT 037436.3      | 25M | 575G            | 4h9m48s  | 37.68 |

## Usage

  `./swifold -i seq1.fasta -j -seq2.fasta  `
  
```
  -i,   --input=<string> Query sequence filename (must be in FASTA format) [REQUIRED]
  
  -j,   --target=<string> Target sequence filename (must be in FASTA format) [REQUIRED]
  
  -M,   --match=<integer> Match score (default: 1)

  -X,   --mismatch=<integer> Mismatch score (default: 3)
  
  -g,   --gap_open=<integer> Gap open penalty (default: 5)

  -e,   --gap_extend=<integer> Gap extend penalty (default: 2)

  -?,   --help Give this help list
        --usage Give a short usage message
```

## Reference

Rucci E., Garcia C., Botella G., De Giusti A., Naiouf M., Prieto-Matias M. (2017) SWIFOLD: Smith-Waterman Acceleration on Intel’s FPGA with OpenCL for Long DNA Sequences. *Under evaluation*

Rucci E., Garcia C., Botella G., De Giusti A., Naiouf M., Prieto-Matias M. (2017) Accelerating Smith-Waterman Alignment of Long DNA Sequences with OpenCL on FPGA. In: Rojas I., Ortuño F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2017. Lecture Notes in Computer Science, vol 10209. Springer, Cham. DOI https://doi.org/10.1007/978-3-319-56154-7_45

## Changelog
* September 01, 2017
Binary code released

## Contact
If you have any question or suggestion, please contact Enzo Rucci (erucci [at] lidi.info.unlp.edu.ar)

