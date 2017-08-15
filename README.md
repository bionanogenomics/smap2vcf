# smap2vcf
# Tool for Transforming SMAP to VCF Format

### Overview 
The SMAP to VCF converter tool is a standalone python script that converts insertion, deletion, duplication, and inversion and translocation breakpoint calls in an SMAP file to dbVar-compliant VCF v4.2 format. There are three required inputs to the script: 1) the SMAP file to convert, 2) the reference (_r.cmap file) which is generated with the SMAP, and 3) the map-to-reference alignment (.xmap file) which is generated with the SMAP. The QUAL score is calculated as -10 times the log base 10 of (1 minus the confidence)) where confidence is the SMAP confidence score for the given call. The QUAL ceiling is set at 20. The uncertainty of the start and end positions of the SV call are provided in the INFO column. This output VCF file can be used for further downstream analysis using any tools that take VCF files as input.

###Usage

usage: smap_to_vcf_v2.py [-h] [-s SMAPPATH] [-r REFCMAPPATH] [-x XMAPPATH]
                         [-n SAMPLE] [-o OUTPUT_PREFIX] [-a REF_ACCESSION]
                         [-b HUMAN_BOOL] [-i EXP_ID]

Stand-alone script to convert BioNano smap file format to vcf.

optional arguments:

    -h, --help        show this help message and exit

    -s SMAPPATH       Path to smap file to convert (required)

    -r REFCMAPPATH    Path to reference cmap which is refered to in the smap
                      supplied as -s argument (optional, default to _r.cmap in
                      smap dir)

    -x XMAPPATH       Path to xmap which is refered to in the smap supplied as
                      -s argument (optional, default to .xmap in smap dir)

    -n SAMPLE         Sample ID name for genotype data (optional, default
                      "Sample1")

    -o OUTPUT_PREFIX  Prefix for output vcf (optional, default to be same as
                      input smap)

    -a REF_ACCESSION  RefSeq assembly accession version (optional, default
                      "GCA_000001405.1")

    -b HUMAN_BOOL     Whether sample is human (optional, default "True")

    -i EXP_ID         dbVar-required experiemnt ID (optional, default "1")



Note:  `python smap_to_vcf.py -h` to see usage on command line

### Requirements
This tool was designed to run with Python 2.7.  

### License
We offer this tool for open source use under the [MIT Software License](https://opensource.org/licenses/MIT). 
