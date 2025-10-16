import pytest
from GeneBankProject.ToolKit.modules.genbank_record import GenBankRecord

# Example minimal GenBank record for testing
minimal_record = {
    "GBSeq_accession-version": "NM_000001.1",
    "GBSeq_feature-table": {
        "GBFeature": [
            {
                "GBFeature_key": "gene",
                "GBFeature_quals": {
                    "GBQualifier": [
                        {"GBQualifier_name": "gene", "GBQualifier_value": "GENE1"},
                        {"GBQualifier_name": "db_xref", "GBQualifier_value": "HGNC:1234"}
                    ]
                }
            },
            {
                "GBFeature_key": "CDS",
                "GBFeature_quals": {
                    "GBQualifier": [
                        {"GBQualifier_name": "protein_id", "GBQualifier_value": "NP_000001"},
                        {"GBQualifier_name": "translation", "GBQualifier_value": "MKT"}
                    ]
                }
            }
        ]
    },
    "GBSeq_sequence": "ATGGAAAAC",
    "GBSeq_definition": "Test gene definition"
}


#function to test the assignment of values to attributes 
def test_initialization():
    record = GenBankRecord(minimal_record)

    #check the attributes are correct and present
    assert record.get_transcript_id == "NM_000001.1"
    assert record.get_protein_id == "NP_000001"
    assert record.get_gene_symbol == "GENE1"
    assert record.get_hgnc_id == "HGNC:1234"
    assert record.get_dna_sequence == "ATGGAAAAC"
    assert record.get_rna_sequence == "AUGGAAAAC"  
    assert record.get_protein_sequence == "MKT

#function to test the fasta output
def test_genbank_output(capsys):
    record = GenBankRecord(minimal_record)
    record.as_genbank("ATGGAAAAC")
    captured = capsys.readouterr()
    assert "1 ATGGAAAAC" in captured.out or "atggaaaac" in captured.out.lower()