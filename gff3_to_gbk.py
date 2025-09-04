####Small python file to add gff annotation into a genbank genome file
###  Used for visualising gff file into , e.g. SnapGene or geneious. 
###Mainly Padloc and defensefinder utility script 
### Usage : python3 gff3_to_gbk.py --> replace the corresponding file names in the last line of the python script. 

from BCBio import GFF
from Bio import SeqIO

def gff3_to_gbk(gff_file, gbk_in, gbk_out):
    """
    Insert features from a GFF3 file into a multi-record GenBank file,
    avoiding duplication of original features.
    """
    # Parse all GenBank records
    gbk_records = list(SeqIO.parse(gbk_in, "genbank"))
    gbk_dict = {r.id: r for r in gbk_records}

    # Parse GFF3 records
    with open(gff_file) as gff_handle:
        for gff_record in GFF.parse(gff_handle, base_dict=gbk_dict):
            if gff_record.id in gbk_dict:
                # Only add new features that were in the GFF
                new_features = [f for f in gff_record.features if f not in gbk_dict[gff_record.id].features]
                gbk_dict[gff_record.id].features.extend(new_features)
            else:
                print(f"⚠ Warning: {gff_record.id} not found in GBK file")

    # Write original GBK records (now with added features)
    SeqIO.write(gbk_records, gbk_out, "genbank")
    print(f"✅ Features from {gff_file} inserted into {gbk_out}")

# Example usage
if __name__ == "__main__":
    gff3_to_gbk("annot_file.gff3", "my_genome.gbff", "output.gbk")
