# Tool for Transforming SMAP to VCF Format

### Overview 
The SMAP to VCF converter tool is a stand alone python script that converts insertions and deletions in a SMAP file to VCF format. There are two required inputs to the scripts: 1) the SMAP file to convert and 2) the reference (_r.cmap file) which is generated with the SMAP. The QUAL score is calculated as -10 time the log base 10 of (1 minus the confidence)) where confidence is the SMAP confidence score for the given call. The uncertainty of the start and end positions of the SV call are provided in the INFO column. This output VCF file can be used for further downstream analysis using any tools that take VCF files as input.

##Usage

smap_to_vcf.py [-h] [-s SMAPPATH] [-r REFCMAPPATH] [-n SAMPLE]

Stand-alone script to convert BioNano smap file format to vcf.

optional arguments:
*  -h, --help      show this help message and exit
*  -s SMAPPATH     Path to smap file to convert (required)
*  -r REFCMAPPATH  Path to reference cmap which is reffered to in the smap supplied as -s argument (required)
*  -n SAMPLE       Sample ID name for genotype data (optional, default "Sample1")

Note:  `python smap_to_vcf.py -h` to see usage on command line

### License
We offer this tool for open source use under the [MIT Software License](https://opensource.org/licenses/MIT). 
