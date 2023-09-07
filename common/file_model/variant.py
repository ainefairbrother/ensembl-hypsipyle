from typing import Any, Mapping, List, Union
import re
import os
import json
import operator
from functools import reduce

def reduce_allele_length(allele_list: List):
    allele_length = -1
    for allele in allele_list:
        if len(allele.value) > allele_length:
            allele_length = len(allele.value)
    return allele_length


class Variant ():
    def __init__(self, record: Any, header: Any) -> None:
        self.name = record.ID[0]
        self.record = record 
        self.header = header
        self.chromosome = record.CHROM         ###TODO: convert the contig name in the file to match the chromosome id given in the payload 
        self.position = record.POS
        self.alts = record.ALT
        self.ref = record.REF
        self.info = record.INFO
        self.type = "Variant"
    
    def get_alternative_names(self) -> List:
        return []
    
    def get_primary_source(self) -> Mapping:
        source = self.header.get_lines("source")[0].value
        if re.search("^dbSNP", source):
            source_id = "dbSNP"
            source_name = "dbSNP"
            source_description = "NCBI db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/snp/"
            source_release =154
            
        elif re.search("^ClinVar", source):
            source_id = "ClinVar"
            source_name = "ClinVar"
            source_description = "ClinVar db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
            source_release = ""

        else:
            source_id = "test"
            source_name = "test"
            source_description = "test db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/test/variation/"
            source_release = ""

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

    def set_allele_type(self, alt_one_bp: bool, ref_one_bp: bool, ref_alt_equal_bp: bool):         
        match [alt_one_bp, ref_one_bp, ref_alt_equal_bp]:
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
    
    def get_allele_type(self, allele: Union[str, List]) -> Mapping :
        if isinstance(allele, str):
            if allele == self.ref:
                allele_type, SO_term = "biological_region","SO:0001411"
            else:
                allele_type, SO_term = self.set_allele_type(len(allele)<2, len(self.ref)<2, len(allele) == len(self.ref))
        elif isinstance(allele, list):
            alt_length = reduce_allele_length(allele)
            allele_type, SO_term = self.set_allele_type(alt_length < 2 , len(self.ref)<2, alt_length == len(self.ref))

        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"http://sequenceontology.org/browser/current_release/term/{SO_term}",
            "source": {
                    "id": "",
                    "name": "Sequence Ontology",
                    "url": "www.sequenceontology.org",
                    "description": "The Sequence Ontology..."
                    }

        }  
    
    def get_slice(self, allele: Union[str, List] ) -> Mapping :

        start = self.position
        length = len(self.ref)
        end = start + length -1
        if allele != self.ref:
            allele_type = self.get_allele_type(allele)
            if allele_type["accession_id"] == "insertion":
                length = 0
                end = start + 1
        
        return {
            "location": {
                "start": start,
                "end": end,
                "length": length
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
    
    def get_alleles(self) -> List:
        variant_allele_list = []
        info_map = self.traverse_csq_info()

        frequency_map = self.format_frequency(",".join(map(str,self.info["FREQ"])).split("|")) if self.info["FREQ"] else {}
        for index,alt in enumerate(self.alts):
            if index+1 <= len(self.alts):
                variant_allele = self.create_variant_allele(info_map, frequency_map, index+1,  alt.value)
                variant_allele_list.append(variant_allele)
        reference_allele = self.create_variant_allele(info_map, frequency_map, 0, self.ref)
        variant_allele_list.append(reference_allele)
        self.set_frequency_flags(variant_allele_list)
        return variant_allele_list

    def create_variant_allele(self, info_map: Mapping, frequency_map: List, allele_index: str, alt: str) -> Mapping:
        """
        The function currently does not include parent-child resolvers.
        Ideally, each field should be a function called upon demand and does lazy loading
        This needs to be rewritten similar to Variant
        """
        name = f"{self.chromosome}:{self.position}:{self.ref}:{alt}"
        min_alt = self.minimise_allele(alt)
        return {
            "name": name,
            "allele_sequence": alt,
            "reference_sequence": self.ref,
            "type": "VariantAllele",
            "allele_type": self.get_allele_type(alt),
            "slice": self.get_slice(alt),
            "phenotype_assertions": info_map[min_alt]["phenotype_assertions"] if min_alt in info_map else [],
            "predicted_molecular_consequences": info_map[min_alt]["predicted_molecular_consequences"] if min_alt in info_map else [],
            "population_frequencies": self.get_population_allele_frequencies(frequency_map, allele_index)
        }
    

    def format_frequency(self, raw_frequency_list: List) -> Mapping:
        freq_map = {}
        for freq in raw_frequency_list:
            key = freq.split(":")[0]
            freq_list = freq.split(":")[1].split(",")
            freq_map[key] = freq_list
        return freq_map
    
    def get_population_allele_frequencies(self, population_map: Mapping, allele_index: int) -> List:
        population_allele_frequencies = []
        for key, pop_list in population_map.items():
            ## Adding only GnomAD population
            if key == "GnomAD":
                if pop_list[allele_index] not in ["None", "."]:
                    population_allele_frequencies.append({
                        "population": key,
                        "allele_frequency": pop_list[allele_index],
                        "is_minor_allele": False,
                        "is_hpmaf": False
                    })
            })
        return population_allele_frequencies
            

    def set_frequency_flags(self, allele_list: List):
        """
        Calculates minor allele frequency by iterating through each allele 
        Assumption: Considers that population is only gnomAD (genomes) for now
        Sets the maf as hpmaf as gnomAD is the only population at the moment
        """
        maf_frequency = -1
        maf_index = -1
        highest_frequency = -1
        highest_frequency_index = -1
        highest_maf_frequency = -1
        highest_maf_frequency_index = -1
        maf_map = {}
        for allele_index, allele in enumerate(allele_list):
            if(len(allele["population_frequencies"]) > 0 ):
                pop = allele["population_frequencies"][0]   
                pop_allele_frequency = float(pop["allele_frequency"])
                if ( pop_allele_frequency > maf_frequency and pop_allele_frequency < highest_frequency ):
                    maf_frequency = pop_allele_frequency
                    maf_index = allele_index
                elif ( pop_allele_frequency > highest_frequency ):
                    maf_frequency = highest_frequency
                    maf_index = highest_frequency_index
                    highest_frequency = pop_allele_frequency
                    highest_frequency_index = allele_index

        if maf_frequency>=0:
            allele_list[maf_index]["population_frequencies"][0]["is_minor_allele"]  = True
            allele_list[maf_index]["population_frequencies"][0]["is_hpmaf"]  = True    
                    
                         


    def get_info_key_index(self, key: str, info_id: str ="CSQ") -> int:
        info_field = self.header.get_info_field_info(info_id).description
        csq_list = info_field.split("Format: ")[1].split("|")
        for index, value in enumerate(csq_list):
            if value == key:
                return index
            
    def get_most_severe_consequence(self) -> Mapping:
        consequence_index = self.get_info_key_index("Consequence")
        consequence_map = {}
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory,'variation_consequence_rank.json')) as rank_file:
            consequence_rank = json.load(rank_file)
        for csq_record in self.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            for cons in csq_record_list[consequence_index].split("&"):
                rank = consequence_rank[cons]
                consequence_map[rank] = cons
        return{
                    "result": consequence_map[min(consequence_map.keys())]  ,
                    "analysis_method": {
                        "tool": "Ensembl VEP",
                        "qualifier": "most severe consequence"
                    }
        } 
        
    def minimise_allele(self, alt: str):
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele = alt
        if len(alt) > len(self.ref):
            minimised_allele = alt[1:] 
        elif len(alt) < len(self.ref):
            minimised_allele = "-"
        return minimised_allele
        

    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        allele_index = self.get_info_key_index("Allele")
        phenotypes_index = self.get_info_key_index("VAR_SYNONYMS")
        feature_type_index = self.get_info_key_index("Feature_type")
        feature_index = self.get_info_key_index("Feature")
        sift_index = self.get_info_key_index("SIFT")
        polyphen_index = self.get_info_key_index("PolyPhen")
        consequence_index = self.get_info_key_index("Consequence")
        spdi_index = self.get_info_key_index("SPDI")
        cadd_index = self.get_info_key_index("CADD_PHRED")
        info_map = {}
        for csq_record in self.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            phenotype = None
            if re.search("--OMIM",csq_record_list[phenotypes_index]):
                phenotype = csq_record_list[phenotypes_index].split("--")[1]
            if csq_record_list[allele_index] in info_map.keys():
                if phenotype:
                    info_map[csq_record_list[allele_index]]["phenotype_assertions"].append(self.create_allele_phenotype_assertion(csq_record_list[feature_index], csq_record_list[feature_type_index], phenotype ))
                info_map[csq_record_list[allele_index]]["predicted_molecular_consequences"].append(self.create_allele_predicted_molecular_consequence(csq_record_list[spdi_index], 
                                                                                                                                               csq_record_list[feature_index], 
                                                                                                                                               csq_record_list[feature_type_index], 
                                                                                                                                               csq_record_list[consequence_index],
                                                                                                                                               csq_record_list[sift_index], 
                                                                                                                                               csq_record_list[polyphen_index],
                                                                                                                                               csq_record_list[cadd_index]
                                                                                                                                   ))
            else:
                info_map[csq_record_list[allele_index]] = {"phenotype_assertions": [], "predicted_molecular_consequences": []} 
                if phenotype:
                    info_map[csq_record_list[allele_index]]["phenotype_assertions"].append(self.create_allele_phenotype_assertion(csq_record_list[feature_index], csq_record_list[feature_type_index], phenotype))
                info_map[csq_record_list[allele_index]]["predicted_molecular_consequences"].append(self.create_allele_predicted_molecular_consequence(csq_record_list[spdi_index], 
                                                                                                                                               csq_record_list[feature_index], 
                                                                                                                                               csq_record_list[feature_type_index], 
                                                                                                                                               csq_record_list[consequence_index],
                                                                                                                                               csq_record_list[sift_index], 
                                                                                                                                               csq_record_list[polyphen_index],
                                                                                                                                               csq_record_list[cadd_index]
                                                                                                                                               ))
        return info_map
    
    def create_allele_phenotype_assertion(self, feature: str, feature_type: str , phenotype: str) -> Mapping:
        return {
            "feature": feature,
            "feature_type": {
                "accession_id": feature_type                
            } ,
            "phenotype": {
                "name": phenotype
            },
            "evidence": []

        }
    
    def create_allele_predicted_molecular_consequence(self, allele: str, feature: str, feature_type: str, consequences: str, sift_score: str, polyphen_score: str, cadd_score: str ) -> Mapping:
        """
        This needs to be designed better, currently all the scores come as args
        Steve suggested that we add prediction results per Gene instead of transcript
        Currently, CADD returns empty
        """
        consequences_list = []
        for cons in consequences.split("&"):
            consequences_list.append(
                {
                    "accession_id": cons
                }
            )
        prediction_results = []

        if cadd_score: 
            cadd_prediction_result = {
                    "result": cadd_score ,
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": "CADD"
                        }

            }
            prediction_results.append(cadd_prediction_result)
        if sift_score:
            sift_prediction_result = {
                    "result": sift_score ,
                    "analysis_method": {
                        "tool": "SIFT",
                        "qualifier": "SIFT"
                    }
                }
            prediction_results.append(sift_prediction_result)
        if polyphen_score:
            polyphen_prediction_result = {
                    "result": polyphen_score ,
                    "analysis_method": {
                        "tool": "PolyPhen",
                        "qualifier": "PolyPhen"
                    }
                }
            prediction_results.append(polyphen_prediction_result)

        
        return {
            "allele_name": allele,
            "feature_stable_id": feature,
            "feature_type": {
                "accession_id": feature_type                
            } ,
            "consequences": consequences_list,
            "prediction_results": prediction_results
        }
    



            

                

            
        


        


