import logging

def convert_dna_to_rna(sequence:str):
	logging.info("Starting the conversion of dna to rna")

	if not isinstance(sequence, str):
		raise TypeError("Sequence must be a string")
		logging.error("DNA sequence is not a string, cannot convert to RNA sequence")

	rna_sequence =""
	sequence_lower = sequence.lower()
	bases = ["a", "g", "c", "t"]

	try:
		for i in sequence_lower:
			if i not in bases:
				raise ValueError(f"Invalid base '{i}' found in DNA sequence")
			if i != "t":
				rna_sequence += i
			
			if i == "t":
				rna_sequence += "u"
	except ValueError as e:
		logging.error(f"Sequence conversion failed")
	
	return rna_sequence