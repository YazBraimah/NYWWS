import os
import sys
import argparse

def aggregated_parse(input_file_lines):
	for lines in input_file_lines[1:]:
		items=lines.split("\t")
		sample_name=items[0].strip(".freyja.variants.tsv")
		strains_list=items[2].split(" ")
		abundances_list=items[3].split(" ")
		for strain_index,each_strain in enumerate(strains_list):

			print_line=[sample_name,each_strain,abundances_list[strain_index]]
			print(",".join(print_line))


def individual_file_parse(input_file_lines):
	sample_name=input_file_lines[0].split("/")[-1].strip(".freyja.variants.tsv\n")
	strains_list=input_file_lines[2].split("\t")[1].strip().split(" ")
	abundances_list=input_file_lines[3].split("\t")[1].strip().split(" ")
	for strain_index,each_strain in enumerate(strains_list):
		print_line=[sample_name,each_strain,abundances_list[strain_index]]
		print(",".join(print_line))

def main():

	#Setting input arguments. Expects an input file. Also optional flag to specify if the input file is an aggregated results file (default), or individual sample results file.
	parser=argparse.ArgumentParser(description="Script to parse Freyja Demix results from an individual or aggregated demix result")
	parser.add_argument("--input_format", help="the expected format for the input demix file. By default set to 'aggregated'. Set to 'individual' for parsing an individual file result. Any other input string and the code will fail", default="aggregated")
	parser.add_argument("filename")
	args=parser.parse_args()


	input_file_lines=open(args.filename).readlines()

	#Runs the aggregated file parser by default
	if args.input_format=="aggregated":
		aggregated_parse(input_file_lines)

	#Alternativ ely runs the individual file parser, if specified in input
	elif args.input_format=="individual":
		individual_file_parse(input_file_lines)

	#Doesn't run if unknown input file format in flag. However, if the wrong input format is given it will probably just output weird results. Should probably test/fix this.
	else:
		print("unknown input format flag. exiting.")
		sys.exit()

main()
