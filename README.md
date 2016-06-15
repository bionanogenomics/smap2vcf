# Tool for Transforming SMAP to VCF Format

### Overview 
cmap2fa.pl is a tool for converting a BioNano multi-color CMAP file of nicking enzyme label positions to an NGS FASTA file of nucleotides sequences. The tool translates each label position to an enzyme-specific sequence (as in its forward strand format), and fills every interval of two labels with “N”s.   

##Usage
pdo `python smap_to_vcf.py -h` to see usage:

wandrews; python smap_to_vcf.py -h
usage: smap_to_vcf.py [-h] [-s SMAPPATH] [-r REFCMAPPATH] [-n SAMPLE]

Stand-alone script to convert BioNano smap file format to vcf.

optional arguments:
*  -h, --help      show this help message and exit
*  -s SMAPPATH     Path to smap file to convert (required)
*  -r REFCMAPPATH  Path to reference cmap which is reffered to in the smap supplied as -s argument (required)
*  -n SAMPLE       Sample ID name for genotype data (optional, default "Sample1")

### License
We offer this tool for open source use under the [MIT Software License](https://opensource.org/licenses/MIT). 
