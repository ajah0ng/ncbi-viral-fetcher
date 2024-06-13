# Option 1: Python Script to annotate viral genomic sequences by Anna Hong (z5429018) on 27/10/23

import sys
from Bio import Entrez, SeqIO

# retrieve_sequence searches NCBI databases given accession_id and returns a genbank record.
# Arguments: accession_id
# Returns: record
def retrieve_sequence(accession_id):

    try:
        Entrez.email = "z5419028@ad.unsw.edu.au"
        handle = Entrez.efetch(db = "nuccore", id = accession_id, rettype= "gb")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record

    except Exception as e:
        print(f"Error: Unable to retrieve genomic sequence for accession ID {accession_id}.")
        sys.exit(1)

# save_fasta takes genbank record and outputs genbank in FASTA file format.
# Arguments: record
# Output: {accession_id}.fa
def save_fasta(record):

    try:
        fastafile_name = f"{accession_id}.fa"
        SeqIO.write(record, fastafile_name, "fasta")

    except Exception as e:
        print(f"Error: Unable to save {accession_id} as a FASTA FILE")
        print(e)
        sys.exit(1)

# create_gtf creates a GTF file given file.
# Arguments: record
# Output: {record.id}.gtf
def create_gtf(record):
    try: 
        gtf_file_name = f"{record.id}.gtf"

        with open(gtf_file_name, "w") as gtf_file:

            for feature in record.features:
                location_list = feature.location.parts # split start and end positions. 
                qualifiers = []  # Create string

                if feature.type == "gene":

                    start = int(location_list[0].start + 1) 

                    if int(location_list[0].end) == int(location_list[-1].end): # If only one start and end
                        end = location_list[0].end
                    else:  # else end codon is the first "end" value + the last "end" value
                        end = int(location_list[0].end) + int(location_list[-1].end)

                    score = "."
                    strand = "+" if feature.location.strand == 1 else "-"
                    frame = "."

                    gtf_file.write(f"{record.id}\t{feature.type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t")

                    # Check for each qualifier, if not, return empty string.
                    gene_id = feature.qualifiers.get("db_xref", [""])[0]
                    gene_name = feature.qualifiers.get("gene", [""])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    regulatory_class = feature.qualifiers.get("regulatory_class", [""])[0]
                    note = feature.qualifiers.get("note", [""])[0]

                     # If present add attribute to qualifier string, else leave empty.
                    if gene_id:
                        qualifiers.append(f"gene_id \"{gene_id}\"")
                    if gene_name:
                        qualifiers.append(f"gene_name \"{gene_name}\"")
                    if locus_tag:
                        qualifiers.append(f"locus_tag \"{locus_tag}\"")
                    if regulatory_class:
                        qualifiers.append(f"regulatory_class \"{regulatory_class}\"")
                    if note:
                        qualifiers.append(f"note \"{note}\"")

                    # Print the attributes
                    if qualifiers:
                        gtf_file.write(";\t ".join(qualifiers))
                        gtf_file.write("\n")

                if feature.type == "CDS":

                    for location in location_list:

                        start = location.start + 1
                        end = location.end
                        score = "."
                        strand = "+" if feature.location.strand == 1 else "-"
                        frame = 1

                        gtf_file.write(f"{record.id}\t{feature.type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t")

                        # Check for each qualifier, if not, return empty string.
                        gene_id = feature.qualifiers.get("db_xref", [""])[0]
                        gene_name = feature.qualifiers.get("gene", [""])[0]
                        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                        regulatory_class = feature.qualifiers.get("regulatory_class", [""])[0]
                        note = feature.qualifiers.get("note", [""])[0]

                        # Reset qualifiers for each CDS feature
                        qualifiers = []

                        if gene_id:
                            qualifiers.append(f"gene_id \"{gene_id}\"")
                        if gene_name:
                            qualifiers.append(f"gene_name \"{gene_name}\"")
                        if locus_tag:
                            qualifiers.append(f"locus_tag \"{locus_tag}\"")
                        if regulatory_class:
                            qualifiers.append(f"regulatory_class \"{regulatory_class}\"")
                        if note:
                            qualifiers.append(f"note \"{note}\"")

                        if qualifiers:
                            gtf_file.write(";\t ".join(qualifiers))
                            gtf_file.write("\n")

    except Exception as e:
        print(f"Error: Unable to create GTF annotation file for {record.id}.")
        print(e)
        sys.exit(1)


# createCodonUsage creates a codon usage table containing the count of each codon.
# Arguments: record
# Output: {accession_id}_codon_usage.txt
def createCodonUsage(record):

    try:
        coding_sequences = [] 

        # Extract CDS sequences into coding_sequences string
        for feature in record.features:
            if feature.type == "CDS":
                coding_sequences.append(feature)

        codon_count = {} # Stores the codon and the count. 

        # Count occurrences of each codon in the coding sequences
        for cds in coding_sequences:

            sequence = str(cds.extract(record.seq))

            # iterate in 3's. 
            for i in range(0, len(sequence), 3): 

                codon = sequence[i:i + 3] 

                if len(codon) == 3: 
                    # Update count in dictionary for codon
                    codon_count[codon] = codon_count.get(codon, 0) + 1

        codon_usage_file = f"{accession_id}_codon_usage.txt"

        with open(codon_usage_file, "w") as txt_file:
            
            # Write in codons in lowercase and count to text file.
            for codon, count in codon_count.items():
                txt_file.write(f"{codon.lower()}\t{count}\n")
                
    except Exception as e:
        print(f"Error: Unable to create codon usage file for {accession_id}.")
        print(e)
        sys.exit(1) 

if __name__ == "__main__":

    args = sys.argv 

    # Check if enough arguments, else exit. 
    if (len(args) == 1):
            print("Missing arguments, please specify accession/taxid")
            exit()

    if (len(args) > 2 ):
            print("Too many arguments, please specify")
            exit()

    # Extract accession_id argument
    accession_id = args[1]
    gb_record = retrieve_sequence(accession_id)
    save_fasta(gb_record)
    create_gtf(gb_record)
    createCodonUsage(gb_record)


