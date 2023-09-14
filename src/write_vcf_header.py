# https://github.com/getzlab/SNVmerger_TOOL/blob/9318b1d10680b0b6cb28684c5f948f55ee76d338/src/write_vcf_header.py from Qing Zhang
import pandas as pd
import argparse
import os

# parse inputs
def parse_args():
	parser = argparse.ArgumentParser(description = "Merging SNP on-cis from MuTect1 result to DNP/TNP/ONP")

	# input files
	parser.add_argument("-ref", required = True, help = "Path to ref.fa")
	parser.add_argument("-ref_fai", required = True, help = "Path to ref.fa.fai")

	# output files	
	parser.add_argument("-o_vcf", required = True, help = "Path to output VCF file that contains SNP/INDEL/XNP")
	parser.add_argument("-case_name", required=True, help = "case name matched with VCF header")
	parser.add_argument("-control_name", required=True, help = "control name matched with VCF header")
	
	args = parser.parse_args()

	# check args
	if not os.path.exists(args.ref):
		raise FileNotFoundError("fasta file not found!")
	if not os.path.exists(args.ref_fai):
		raise FileNotFoundError("fasta index file not found!")

	return args


args = parse_args()
header_config = {
    'fileformat' : 'VCFv4.1',
    'FILTER' : {
        'PASS' : {
            'Description' : "Accept as a confident somatic mutation"
        },
        'REJECT' : {
            'Description' : "Rejected as a confident somatic mutation"
        }
    },
    'FORMAT' : {
        'AD' : {
            'Number' : ".",
            'Type' : 'Integer',
            'Description' : 'Allelic depths for the ref and alt alleles in the order listed'
        }
    },
	'INFO' : {
		'SOMATIC' : {
			'Number' : 0,
			'Type' : 'Flag',
			'Description' : "Somatic event"
		}
	},
    'reference' : args.ref, 
	'normal_sample' : args.control_name,
	'tumor_sample' : args.case_name
}

with open(args.o_vcf, 'w') as vcf:
	for k,v in header_config.items():
		if k in ["FILTER", "FORMAT", "INFO"]:
			if v is None:
				vcf.write("##{}=<.>\n".format(k))
				continue
			for sk, sv in v.items():
				line = "##{}=<".format(k)
				line += "ID={},".format(sk)
				line += ",".join(["{}=\"{}\"".format(ssk, ssv) if ssk=="Description" else "{}={}".format(ssk, ssv) for ssk, ssv in sv.items()])
				vcf.write(line+'>\n')
		elif k == "reference":
			# we need contigs <name, length> from fasta index
			contiglen_df = pd.read_csv(args.ref_fai, sep = '\t', usecols=[0,1], names=['contig', 'len'])
			for i,r in contiglen_df.iterrows():
				vcf.write("##contig=<ID={},length={}>\n".format(r['contig'], r['len']))
		if k in ['fileformat', 'reference', 'normal_sample', 'tumor_sample']:
			vcf.write("##{}={}\n".format(k, v))
	vcf.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}	{}\n".format(args.control_name, args.case_name))