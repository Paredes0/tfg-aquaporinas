import os

def update_prune_file(fasta_path, prune_path):
    # Extract IDs from fasta
    new_ids = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # specific cleaning might be needed if there are descriptions
                    # taking everything up to the first whitespace
                    header = line.strip()[1:]
                    # The user said "los encabezados", which usually implies the whole ID line.
                    # looking at prune.txt, they are just IDs.
                    # looking at the fasta, some headers have spaces?
                    # No, looking at lines like ">mRNA_54367-Fxa2Ag00184" it seems valid as is.
                    # Standard fasta ID is until the first space.
                    # Let's double check if there are spaces in the headers I saw.
                    # >Fxa1Ag01329
                    # >mRNA_54367-Fxa2Ag00184
                    # >Fxa3Ag00839
                    # It seems safe to take the whole line after >
                    new_ids.append(header)
    except FileNotFoundError:
        print(f"Error: Could not find {fasta_path}")
        return

    # Read prune.txt content
    prune_lines = []
    try:
        with open(prune_path, 'r') as f:
            lines = f.readlines()
            # Find the DATA line
            data_index = -1
            for i, line in enumerate(lines):
                if line.strip() == 'DATA':
                    data_index = i
                    break
            
            if data_index != -1:
                # Keep everything up to DATA and include it
                prune_lines = lines[:data_index+1]
                # Add a newline if it was there or just for good formatting
                prune_lines.append('\n')
            else:
                print("Error: Could not find DATA line in prune.txt")
                return

    except FileNotFoundError:
        print(f"Error: Could not find {prune_path}")
        return

    # Write back to prune.txt
    with open(prune_path, 'w') as f:
        f.writelines(prune_lines)
        for new_id in new_ids:
            f.write(f"{new_id}\n")
    
    print(f"Successfully updated {prune_path} with {len(new_ids)} IDs from {fasta_path}")

# Define paths
fasta_file = r"c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\analisis proteinas aquaporina\consenso_aqp.fasta"
prune_file = r"c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\prune.txt"

update_prune_file(fasta_file, prune_file)
