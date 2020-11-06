from sys import argv

if len(argv) != 3:
	print("Run with 2 arguments: input_GFF.gff3 output_FASTA.fasta")
	exit(0)
else:
	pass

script, input_gff, input_fasta = argv

gff_file = open(input_gff)
fasta_file = open(input_fasta)

gene_count = 0
gene_length_total = 0

# Go through lines in GFF.
for line in gff_file:
	if "#" in line:
		continue

	current_line = line.split('\t')
	if gene in current_line[2]:
		gene_count += 1
		start = current_line[3]
		end = current_line[4]
		gene_length = end - start
		gene_length_total += gene_length
	
	
		
line_length_total = 0	
	
for line in fasta_file:
	if ">" in line:
		continue
		
	line_length = len(line)
	line_length_total += line_length
	
print(f"The total length of the genome is {line_length_total}.")
print(f"The number of genes is {gene_count}.")
print(f"The total gene length is {gene_length_total}.")
print (f"The percentage of the genome which is coding is {gene_length_total}/{line_length_total}.")
		
