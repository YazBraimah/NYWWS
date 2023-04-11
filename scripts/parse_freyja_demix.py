
infile = snakemake.input[0]
outfile = snakemake.output[0]

with open(infile) as freyja, open(outfile, "w") as result:

    # write header line in output file
    header = ["sample_id", "variant", "frequency"]
    result.write(",".join(header) + "\n")

    # skip the header line of the input file
    freyja.readline() 

    # write data
    for line in freyja:
        items = line.split("\t")
        sample_name = items[0].strip(".freyja.variants.tsv")
        strains_list = items[2].split(" ")
        abundances_list = items[3].split(" ")
        for strain_index, each_strain in enumerate(strains_list):
            print_line = [sample_name, each_strain, abundances_list[strain_index]]
            result.write(",".join(print_line) + "\n")
