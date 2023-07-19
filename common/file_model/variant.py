from typing import Any, Mapping
import re
import operator
from functools import reduce

def reduce_allele(acc,alt):
    return (len(acc.value) < 2) & (len(alt.value) < 2)    

def reduce_allele_length(acc,alt):
    return len(acc.value) == len(alt.value)

class Variant ():
    def __init__(self, record: Any) -> None:
        self.name = record.ID[0]
        self.record = record 
        self.chromosome = record.CHROM         ###TODO: convert to numerical if not
        self.position = record.POS
        self.alts = record.ALT
        self.ref = record.REF
        self.info = record.INFO
        self.primary_source = self.get_primary_source()
        self.type = "Variant"
        self.allele_type = self.get_allele_type()
        self.slice = self.get_slice()
        self.alleles = self.get_alleles()
    
    def get_primary_source(self) -> Mapping:
        if re.search("^rs", self.name):
            source_id = "dbSNP"
            source_name = "dnSNP"
            source_description = "NCBI db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/snp/"
            source_release =154
            
        else:
            source_id = "ClinVar"
            source_name = "ClinVar"
            source_description = "ClinVar db of human variants"
            source_url = "https://clinvar-url"
            source_release = ""

        ### TODO: This is not elegant, need to serialise?
        return {
            "accession_id": self.name,
            "name": self.name,
            "description": "",
            "assignment_method": {
                                "type": "DIRECT",
                                "description": "A reference made by an external resource of annotation to an Ensembl feature that Ensembl imports without modification"
                            },
            "url": f"{source_url}{self.name}",
            "source": {
                        "id" : f"{source_id}",
                        "name": f"{source_name}",
                        "description": f"{source_description}",
                        "url":  f"{source_url}",
                        "release": f"{source_release}"
                        }
        }

    def set_allele_type(self, condition1, condition2, condition3):         
        match [condition1,condition2,condition3 ]:
            case [True, True, True]: 
                allele_type = "SNV" 
                SO_term = "SO:0001483"

            case [True, False, False]: 
                allele_type = "deletion" 
                SO_term = "SO:0000159"

            case [False, True, False]: 
                allele_type = "insertion"
                SO_term = "SO:0000667"

            case [False, False, False]: 
                allele_type = "indel"
                SO_term = "SO:1000032"

            case [False, False, True]: 
                allele_type = "substitution"
                SO_term = "SO:1000002"
                
        return allele_type, SO_term
    
    def get_allele_type(self, allele=None) -> Mapping :
        if allele:
            allele_type, SO_term = self.set_allele_type(len(allele)<2, len(self.ref)<2, len(allele) == len(self.ref))
        else:
            allele_type, SO_term = self.set_allele_type(reduce(reduce_allele,self.alts), len(self.ref)<2, reduce(reduce_allele_length,self.alts) == len(self.ref))
        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"https://sequenceontology.org/browser/current_release/term/{SO_term}",
            "source": {
                    "id": "",
                    "name": "Sequence Ontology",
                    "url": "www.sequenceontology.org",
                    "description": "The Sequence Ontology..."
                    }

        }
    
    
    
    def get_slice(self) -> Mapping :
        start = self.position
        end = self.position
        return {
            "location": {
                "start": start,
                "end": end,
                "length": end - start
            },
            "region": {
                "name": self.chromosome,
                "code": "chromosome",
                "topology": "linear",
                "so_term": "SO:0001217"
            },
            "strand": {
                "code": "forward",
                "value": 1
            }
        }
    
    def get_alleles(self):
        variant_allele_list = []
        for alt in self.alts:
            variant_allele = self.create_variant_allele(alt)
            variant_allele_list.append(variant_allele)
        return variant_allele_list

    def create_variant_allele(self, alt):
        name = f"{self.chromosome}:{self.position}:{self.ref}:{alt.value}"
        return {
            "name": name,
            "allele_sequence": alt.value,
            "reference_sequence": self.ref,
            "type": "VariantAllele",
            "allele_type": self.get_allele_type(alt.value),
            "slice": self.get_slice()
        }
        


