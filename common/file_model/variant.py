from typing import Any, Mapping, List, Union
import re
import os
import json
import operator
from functools import reduce
from common.file_model.variant_allele import VariantAllele

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
        self.vep_version = re.search("v\d+", self.header.get_lines("VEP")[0].value).group()
        self.population_map = {}
    
    def get_alternative_names(self) -> List:
        return []
    
    def get_primary_source(self) -> Mapping:
        """
        Fetches source from variant INFO columns
        Fallsback to fetching from the header
        """
        
        try:
            if "SOURCE" in self.info:
                source = self.info["SOURCE"]
            else:
                source = self.header.get_lines("source")[0].value

            if re.search("^dbSNP", source):
                source_id = "dbSNP"
                source_name = "dbSNP"
                source_description = "NCBI db of human variants"
                source_url = "https://www.ncbi.nlm.nih.gov/snp/"
                source_url_id = source_url
                source_release = 156
                variant_id = self.name
            
            elif re.search("^EVA", source):
                source_id = "EVA"
                source_name = "EVA"
                source_description = "European Variation Archive"
                source_url = "https://www.ebi.ac.uk/eva"
                source_url_id = "https://www.ebi.ac.uk/eva/?variant&accessionID="
                source_release = "release_5"
                variant_id = self.name

            elif re.search("^Ensembl", source):
                source_id = "Ensembl"
                source_name = "Ensembl"
                source_description = "Ensembl"
                source_url = "https://beta.ensembl.org"
                source_url_id = "https://beta.ensembl.org/"
                source_release = "110" # to be fetched from the file
                variant_id = f"{self.chromosome}:{self.position}:{self.name}"
            
            

        except Exception as e:
            return None 

        return {
            "accession_id": self.name,
            "name": self.name,
            "description": f"{source_description}",
            "assignment_method": {
                                "type": "DIRECT",
                                "description": "A reference made by an external resource of annotation to an Ensembl feature that Ensembl imports without modification"
                            },
            "url": f"{source_url_id}{variant_id}",
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

        
        for index,alt in enumerate(self.alts):
            if index+1 <= len(self.alts):
                variant_allele = VariantAllele(index+1,  alt.value, self)
                variant_allele_list.append(variant_allele)
        reference_allele = VariantAllele(0, self.ref, self)
        variant_allele_list.append(reference_allele)
        return variant_allele_list
    
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
                consequence_map[int(rank)] = cons
        return{
                    "result": consequence_map[min(consequence_map.keys())]  ,
                    "analysis_method": {
                        "tool": "Ensembl VEP",
                        "qualifier": "most severe consequence"
                    }
        } 

    def get_gerp_score(self) -> Mapping:
        csq_record = self.info["CSQ"]
        csq_record_list = csq_record[0].split("|")
        if self.get_info_key_index("Conservation") is not None:
            gerp_index = self.get_info_key_index("Conservation") 
            gerp_prediction_result = {
                    "score": csq_record_list[gerp_index] ,
                    "analysis_method": {
                        "tool": "GERP",
                        "qualifier": "GERP"
                    }
                } if csq_record_list[gerp_index] else {}
            return gerp_prediction_result
    
    def get_ancestral_allele(self) -> Mapping:
        csq_record = self.info["CSQ"]
        csq_record_list = csq_record[0].split("|")
        if self.get_info_key_index("AA") is not None:
            aa_index = self.get_info_key_index("AA")
            aa_prediction_result = {
                    "result": csq_record_list[aa_index] ,
                    "analysis_method": {
                        "tool": "AncestralAllele",
                        "qualifier": "",
                        "version": "110" #self.vep_version
                    }
                } if csq_record_list[aa_index] and csq_record_list[aa_index]!="."  else {}
            return aa_prediction_result
    
    def get_info_key_index(self, key: str, info_id: str ="CSQ") -> int:
            info_field = self.header.get_info_field_info(info_id).description
            csq_list = info_field.split("Format: ")[1].split("|")
            for index, value in enumerate(csq_list):
                if value == key:
                    return index   
                
    def traverse_population_info(self) -> Mapping:
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory,'populations.json')) as pop_file:
            pop_mapping = json.load(pop_file)
        population_frequency_map = {}
        for csq_record in self.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            allele_index = self.get_info_key_index("Allele") 
            if csq_record_list[allele_index] is not None:
                if csq_record_list[allele_index] not in population_frequency_map.keys():
                    population_frequency_map[csq_record_list[allele_index]] = {}
                for pop_key, pop in pop_mapping.items():
                    for sub_pop in pop:
                        if sub_pop["name"] in population_frequency_map[csq_record_list[allele_index]]:
                            continue
                        allele_count = allele_number = allele_frequency = None
                        for freq_key, freq_val in sub_pop["frequencies"].items():
                            col_index = self.get_info_key_index(freq_val)
                        if col_index is not None and csq_record_list[col_index] is not None:
                            if freq_key == "af":
                                allele_frequency = csq_record_list[col_index] 
                            elif freq_key == "an":
                                allele_number = csq_record_list[col_index] or None
                            elif freq_key == "ac":
                                allele_count = csq_record_list[col_index] or None
                            else:
                                raise Exception('Frequency metric is not recognised')
                            population_frequency = {
                                                "population": sub_pop["name"],
                                                "allele_frequency": allele_frequency or -1,
                                                "allele_count": allele_count,
                                                "allele_number": allele_number,
                                                "is_minor_allele": False,
                                                "is_hpmaf": False
                                            }

                            population_frequency_map[csq_record_list[allele_index]][sub_pop["name"]] = population_frequency
        return population_frequency_map
    
    def set_frequency_flags(self):
        """
        Calculates MAF (minor allele frequency) and  HPMAF by iterating through each allele 
        """
        allele_list = self.get_alleles()
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory,'populations.json')) as pop_file:
            pop_mapping = json.load(pop_file)
        pop_names = []
        for pop in pop_mapping.values():
            pop_names.extend([sub_pop["name"] for sub_pop in pop])
        hpmaf = []
        pop_frequency_map = self.traverse_population_info()

        for pop_name in pop_names:
            by_population = []
            for allele in allele_list:
                if allele.minimise_allele(allele.alt) in pop_frequency_map:
                    pop_freqs = pop_frequency_map[allele.minimise_allele(allele.alt)].values()
                # Calculate reference allele from parsed frequency info of all alleles
                if allele.get_allele_type()["accession_id"] == "biological_region" and len(by_population)>0 :
                     allele_frequency_ref = 1 - sum(list(zip(*by_population))[0])
                     if allele_frequency_ref <= 1 and allele_frequency_ref >= 0:
                        population_frequency_ref = {
                                                "population": pop_name,
                                                "allele_frequency": allele_frequency_ref ,
                                                "allele_count": None,
                                                "allele_number": None,
                                                "is_minor_allele": False,
                                                "is_hpmaf": False
                                            }
                        minimised_allele = allele.minimise_allele(allele.alt)
                        if minimised_allele not in pop_frequency_map:
                            pop_frequency_map[minimised_allele] = {}
                        pop_frequency_map[minimised_allele][pop_name] = population_frequency_ref
                        by_population.append([allele_frequency_ref,minimised_allele,pop_name])
                
                else:       
                    for pop_freq in pop_freqs:
                        if pop_name == pop_freq["population"] and float(pop_freq["allele_frequency"]) > -1:
                            minimised_allele = allele.minimise_allele(allele.alt)
                            by_population.append([float(pop_freq["allele_frequency"]),minimised_allele, pop_name])
            
                

            by_population_sorted = sorted(by_population, key=lambda item: item[0])
            if len(by_population_sorted) >= 2:
                maf_frequency, maf_allele, maf_population = by_population_sorted[-2]
                pop_frequency_map[maf_allele][maf_population]["is_minor_allele"] = True
                hpmaf.append([maf_frequency,maf_allele,maf_population])
        if len(hpmaf) > 0:
            hpmaf_sorted = sorted(hpmaf, key=lambda item: item[0])
            _, hpmaf_allele, hpmaf_population = hpmaf_sorted[-1]
            pop_frequency_map[hpmaf_allele][hpmaf_population]["is_hpmaf"] = True
        return pop_frequency_map
    


