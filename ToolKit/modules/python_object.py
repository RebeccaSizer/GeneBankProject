#load the utils needed
from GeneBankProject.ToolKit.utils.entrez_efetch import fetch_transcript_record
from GeneBankProject.ToolKit.modules.dna_to_rna import convert_dna_to_rna
import logging

logging.basicConfig(filename="GeneBankProject/logs/python_object.log", level=logging.DEBUG, format="%(asctime)s:%(levelname)s:%(message)s")

class GenBankRecord():

    def __init__(self, record):

        #add a statement to show initiation. Returns Unknown if no GBSeq_accession-version found. 
        logging.info(f"Starting to process record: {record.get('GBSeq_accession-version', 'Unkown')}")


        # --- transcript version ---
        try:
            self.transcript_id = record["GBSeq_accession-version"]
            logging.debug(f"Found transcript ID: {self.transcript_id}")
        except KeyError:
            logging.warning(f"No transcript ID found")
        except AttributeError as e:
            logging.error(f"Transcript ID not a string: {e}")
    

        # --- protein ID for the translation (NP_...) ---
        try:
            for feature in record['GBSeq_feature-table']['GBFeature']:
                if feature['GBFeature_key'] == 'CDS':
                    for qual in feature.get('GBFeature_quals', {}).get('GBQualifier', []):
                        if qual.get('GBQualifier_name') == 'protein_id':
                            protein_id = qual.get('GBQualifier_value')
                            self.protein_id = protein_id
                            logging.debug(f"Found protein ID: {self.protein_id}")
        except KeyError:
            logging.warning(f"No Protein ID found")
        except AttributeError as e:
            logging.error(f"Protein ID not a string: {e}")


        # --- Gene symbol ---
        try:
            for line in record['GBSeq_feature-table']['GBFeature']:
                if line['GBFeature_key'] == 'gene':
                    for qual in line.get('GBFeature_quals', {}).get('GBQualifier', []):
                        if qual.get('GBQualifier_name') == 'gene':
                            self.gene_symbol = qual.get('GBQualifier_value')
                            logging.debug(f"Found gene symbol: {self.gene_symbol}")
                            break
        except KeyError:
            logging.warning(f"Gene symbol not found")
        except AttributeError as e:
            logging.error(f"Gene symbol not a string: {e}")


        # --- HGNC ID ---
        try:
            for line in record['GBSeq_feature-table']['GBFeature']:
                if line['GBFeature_key'] == 'gene':
                    for qual in line.get('GBFeature_quals', {}).get('GBQualifier', []):
                        if qual.get('GBQualifier_name') == 'db_xref' and 'HGNC' in qual.get('GBQualifier_value', ''):
                            value = qual.get('GBQualifier_value')
                            self.hgnc_id = value.replace("HGNC:HGNC:", "HGNC:")
                            logging.debug(f"Found HGNC ID: {self.hgnc_id}")
                            break
        except KeyError:
            logging.warning(f"HGNC ID not found")
        except AttributeError as e:
            logging.error(f"HGNC ID not a string: {e}")


        # --- Transcript Sequence ---
        try:
            dna_sequence_original = record['GBSeq_sequence']
            logging.debug(f"DNA Sequence found (Length: {len(dna_sequence_original)})")
        except KeyError:
            logging.warning(f"DNA sequence not found")
        except AttributeError as e:
            logging.error(f"DNA sequence not a string {e}")
        
        try:
            self.dna_sequence = dna_sequence_original.upper()
            logging.info(f" DNA sequence converted to upper case (length: {len(self.dna_sequence)})")
        except Exception as e:
            logging.error(f"could not convert DNA sequence to upper case: {e}")
            

        # --- turn the dna sequence to the rna sequence ---
        try:
            self.rna_sequence = convert_dna_to_rna(self.dna_sequence)
            logging.info(f"convert the dna sequence to an rna sequence")
        except Exception as e:
            logging.error(f"could not convert dna to rna: {e}")

        # --- protein sequence ---
        try:
            for line in record['GBSeq_feature-table']['GBFeature']:
                if line['GBFeature_key'] == 'CDS':
                    for qual in line.get('GBFeature_quals', {}).get('GBQualifier', []):
                        if qual.get('GBQualifier_name') == 'translation':
                            self.protein_sequence = qual.get('GBQualifier_value')
                            logging.debug(f"Protein sequence found (length: {len(self.protein_sequence)})")
                            break
        except KeyError:
            logging.warning(f"Protein sequence not found")
        except Exception as e:
            logging.error(f"Protein sequence is not a string: {e}")

    
        #set the fasta line 
        self.fasta_line = record['GBSeq_definition']      #this needs changing to make it more versatile

    #method to return transcript_id
    @property
    def get_transcript_id(self):
        return self.transcript_id

    #method to get the protein id
    @property
    def get_protein_id(self):
        return self.protein_id

    @property
    #method to return gene symbol
    def get_gene_symbol(self):
        return self.gene_symbol

    @property
    #method to returnhgnc id
    def get_hgnc_id(self):
        return self.hgnc_id
    
    @property
    def get_dna_sequence(self):
        return self.dna_sequence

    @property
    def get_rna_sequence(self):
        return self.rna_sequence

    @property
    def get_protein_sequence(self):
        return self.protein_sequence

    def as_fasta(self, transcript_seq):
        print(f"> {self.transcript_id} {self.fasta_line}")
        count = 0
        new_seq = ""

        for i in transcript_seq:
            count += 1
            new_seq += i

            if count%60 == 0:
                new_seq += f"\n"
        
        print(new_seq)

    def as_genbank(self, transcript_seq): 
        count = 0 
        gap = 10
        block = 6
        new_seq = ""
        sequence_lower = transcript_seq.lower()
        column_count = 0
        length = str(len(sequence_lower))
        column_width = int(len(length)) + 1

        for i in sequence_lower:
            count += 1

            if count == 1:
                new_seq += f"{count:< {column_width + 1}}"
            
            new_seq += i

            if count % gap == 0:
                new_seq += " "
            
            if count % (gap * block) == 0 and count < len(transcript_seq):
                new_seq += f"\n{count:< {column_width + 1}}"
        
        print(new_seq)


#test the script 
if __name__ == "__main__":

    #first you need to ask the user for what reference sequence they want 
    if False:
        transcript_id = str(input("Please type the transcript ID: "))
    else:
        transcript_id = "NM_000277.3" # Example usage: fetch a GenBank transcript record (e.g. COL5A1 mRNA RefSeq)

    #fetch the GenBank record using the reference_seq specified above 
    record = fetch_transcript_record(transcript_id)
    output = GenBankRecord(record)
    
    print(output.get_transcript_id)
    print(output.get_gene_symbol)
    print(output.get_hgnc_id)
    print(output.get_protein_id)
    print(output.get_protein_sequence)
    print(output.get_dna_sequence)
    print(output.get_rna_sequence)
    output.as_genbank(output.get_rna_sequence)

"""
    print(record)
    # Top-level metadata fields from the GenBank record
    print(record['GBSeq_accession-version'])   # Accession with version, e.g. "NM_000093.5"
    print(record['GBSeq_definition'])          # Definition line, e.g. "collagen alpha-1(V) chain (COL5A1), mRNA"
    print(record['GBSeq_keywords'])            # Keywords list, e.g. ["RefSeq", "mRNA", "collagen"]

    # Print the raw nucleotide sequence in uppercase
    print(record['GBSeq_sequence'].upper())

    # Iterate over all annotated features in the feature table
    for line in record['GBSeq_feature-table']['GBFeature']:

        # --- Gene feature ---
        if line['GBFeature_key'] == 'gene':
            # Gene symbol, e.g. "COL5A1"
            print(line['GBFeature_quals']['GBQualifier'][0]['GBQualifier_value'])
            # HGNC ID (sometimes appears as "HGNC:HGNC:2197", so fix formatting)
            print(line['GBFeature_quals']['GBQualifier'][4]['GBQualifier_value']
                  .replace("HGNC:HGNC:", "HGNC:"))

        # --- Coding sequence (CDS) feature ---
        elif line['GBFeature_key'] == 'CDS':
            # Start coordinate of coding region on the transcript
            print(line['GBFeature_intervals']['GBInterval']['GBInterval_from'])
            # End coordinate of coding region on the transcript
            print(line['GBFeature_intervals']['GBInterval']['GBInterval_to'])

    # Keywords again (redundant, but ensures access outside the loop)
    print(record['GBSeq_keywords'])

    # Print all top-level keys in the record dict for exploration/debugging
    print(record.keys())
"""