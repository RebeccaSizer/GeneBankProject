#load the utils needed
from GeneBankProject.ToolKit.utils.entrez_efetch import fetch_transcript_record
from GeneBankProject.ToolKit.modules.dna_to_rna import convert_dna_to_rna
class GenBankRecord():

    def __init__(self, record):

        # --- transcript version ---
        self.transcript_id = record["GBSeq_accession-version"]

        # --- Gene symbol ---
        for line in record['GBSeq_feature-table']['GBFeature']:
            if line['GBFeature_key'] == 'gene':
                self.gene_symbol = line['GBFeature_quals']['GBQualifier'][0]['GBQualifier_value']
            
        # --- HGNC ID ---
                self.hgnc_id = line['GBFeature_quals']['GBQualifier'][4]['GBQualifier_value'].replace("HGNC:HGNC:", "HGNC:")

        # --- protein sequence ---
        
        for line in record['GBSeq_feature-table']['GBFeature']:
            if line['GBFeature_key'] == 'CDS':
                self.protein_sequence = line['GBFeature_quals']['GBQualifier'][11]['GBQualifier_value']    

        # --- Transcript Sequence ---
        self.dna_sequence = record['GBSeq_sequence'].upper()

        # --- turn the dna sequence to the rna sequence ---
        self.rna_sequence = convert_dna_to_rna(self.dna_sequence)
    
        #set the fasta line 
        self.fasta_line = record['GBSeq_definition']

    #method to return transcript_id
    @property
    def get_transcript_id(self):
        return self.transcript_id

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
        new_seq = ""
        sequence_lower = transcript_seq.lower()
        column_count = 0
        number_width = str(transcript_seq/60)
        column_width = len(number_width)

        for i in sequence_lower:
            count += 1

            if count == 1:
                new_seq += f"{count:<column_width}"
            
            new_seq += i

            if count % gap == 0:
                new_seq += " "
            
            if count % (gap * block) == 0 and count < len(transcript_seq):
                new_seq += f"\n{count:<column_width - 1}"
        
        print(new_seq)


#test the script 
if __name__ == "__main__":

    #first you need to ask the user for what reference sequence they want 
    if False:
        transcript_id = str(input("Please type the transcript ID: "))
    else:
        transcript_id = "NM_000093.5" # Example usage: fetch a GenBank transcript record (e.g. COL5A1 mRNA RefSeq)

    #fetch the GenBank record using the reference_seq specified above 
    record = fetch_transcript_record(transcript_id)
    output = GenBankRecord(record)
    
    print(output.get_transcript_id)
    output.as_genbank(output.rna_sequence)
    



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