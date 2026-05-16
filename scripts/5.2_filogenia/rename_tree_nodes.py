
import csv
import re
import os

def rename_newick_nodes(tree_path, tabular_path, output_path):
    # 1. Build Mapping
    # We need to map every possible ID that appears in the tree to its subfamily.
    # The tree might contain gene_ids or mRNA_ids or exonerate_ids.
    # We will build a map that covers all of them pointing to the subfamily.
    
    id_map = {}
    
    try:
        with open(tabular_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            for row in reader:
                if len(row) < 4:
                    continue
                
                # Columns: 0=gene_id, 1=mRNA_gff, 2=exonerate_id, 3=subfamily
                gene_id = row[0].strip()
                mrna_id = row[1].strip()
                exonerate_id = row[2].strip()
                subfamily = row[3].strip()
                
                if gene_id: id_map[gene_id] = subfamily
                if mrna_id: id_map[mrna_id] = subfamily
                if exonerate_id: id_map[exonerate_id] = subfamily
                
    except FileNotFoundError:
        print(f"Error: Tabular file not found: {tabular_path}")
        return

    # 2. Read Tree File
    try:
        with open(tree_path, 'r') as f:
            tree_content = f.read()
    except FileNotFoundError:
        print(f"Error: Tree file not found: {tree_path}")
        return

    # 3. Perform Replacements
    # To avoid substring replacement issues (e.g. replacing 'ID1' in 'ID10'),
    # we sort keys by length descending.
    # Also, since Newick format uses specific delimiters, likely we are safe, 
    # but let's be careful.
    
    # We only replace IDs that are actually in the tree content to save time/risk
    keys_in_tree = [k for k in id_map.keys() if k in tree_content]
    
    # Sort by length descending
    keys_in_tree.sort(key=len, reverse=True)
    
    print(f"Found {len(keys_in_tree)} IDs from the table in the tree file.")
    
    new_content = tree_content
    count = 0
    for key in keys_in_tree:
        # Simple replace is risky if IDs are substrings of each other AND not delimited?
        # Given the specific IDs (Fxa...), usually safe if sorted by length.
        # But let's check if the replacement (subfamily) contains the key? No, FaPIP2 doesn't look like Fxa...
        
        # We replace the ID with the Subfamily
        if key in new_content:
            new_content = new_content.replace(key, id_map[key])
            count += 1
            
    # 4. Write Output
    with open(output_path, 'w') as f:
        f.write(new_content)
        
    print(f"Replaced {count} IDs. Saved to {output_path}")

tree_file = r"c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\aqp_subfamily_newick.txt"
tab_file = r"c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\analisis proteinas aquaporina\tabla_Aquaporinas_traduccion.tabular"
output_file = r"c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\aqp_subfamily_renamed.newick"

rename_newick_nodes(tree_file, tab_file, output_file)
