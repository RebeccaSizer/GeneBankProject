#load the utils needed
from GeneBankProject.ToolKit.utils.entrez_efetch import fetch_transcript_record

class GenBankRecord():

    def __init__(self, record):

        # --- transcript version ---
        self.transcript_id = record["GBSeq_accession-version"]

        # --- Gene symbol ---
        for line in record['GBSeq_feature-table']['GBFeature']:
            if line['GBFeature_key'] == 'gene':
                self.gene_symbol = line['GBFeature_quals']['GBQualifier'][0]['GBQualifier_value']
            
        # --- HGNC ID ---
                self. hgnc_id = line['GBFeature_quals']['GBQualifier'][4]['GBQualifier_value'].replace("HGNC:HGNC:", "HGNC:")

        # --- Transcript Sequence ---
        self.transcript_sequence = record['GBSeq_sequence'].upper()

    def __str__(self):
        return f"{self.transcript_id}\n{self.gene_symbol}\n{self.hgnc_id}\n{self.transcript_sequence}"


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
    print(output)





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