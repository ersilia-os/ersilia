# Martino suggested:
# https://pypi.org/project/crc64iso/
from bioservices.uniprot import UniProt


class ProteinIdentifier(object):
    def __init__(self):
        pass
        self.seguid = seguid
        self.uniprot = UniProt(verbose=False)

    def sequence_from_uniprot(self, uniprot_ac):
        """Returns protein sequence from uniprot identifier"""
        try:
            return self.uniprot.get_fasta_sequence(uniprot_ac)
        except ValueError:
            return None

    @staticmethod
    def protein_identifier_resolver():
        """Returns protein sequence of a given identifier, using"""
        pass  # TODO

    @staticmethod
    def encode(sequence):
        """Protein seguid checksum based on amino-acid sequence"""
        return str(self.seguid(sequence))
