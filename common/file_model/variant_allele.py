from typing import Any, Mapping, List, Union
import re
import os
import json
import operator
from functools import reduce

class VariantAllele():
    def __init__(self, allele_index: str, alt: str, variant:dict) -> None:
        
        self.name = f"{variant.chromosome}:{variant.position}:{variant.ref}:{alt}"
        self.variant = variant
        self.allele_index = allele_index
        self.alt = alt
        self.type = "VariantAllele"
        self.allele_sequence = alt
        self.reference_sequence = variant.ref
        self.population_map = []
        self.info_map = self.traverse_csq_info()

    def get_allele_type(self):
        #TODO: change this to VariantAllele level
        return self.variant.get_allele_type(self.alt)

    def get_alternative_names(self):
        return self.variant.get_alternative_names()

    def get_slice(self):
        #TODO: review this to change to VariantAllele level
        return self.variant.get_slice(self.alt)
    
    def get_phenotype_assertions(self):
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["phenotype_assertions"] if min_alt in self.info_map else []

    def get_predicted_molecular_consequences(self):
        ## compute info_map only for these methods, first check if non-empty and compute only then
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["predicted_molecular_consequences"] if min_alt in self.info_map else []
    
    def get_prediction_results(self):
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["prediction_results"] if min_alt in self.info_map else []
    
    def get_population_allele_frequencies(self):
        min_alt = self.minimise_allele(self.alt)
        population_map = self.variant.set_frequency_flags()
        return population_map[min_alt].values() if min_alt in population_map else []


    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        column_list = ["Allele", "PHENOTYPES", "Feature_type", "Feature", "Consequence",
                        "SIFT", "PolyPhen", "SPDI", "CADD_PHRED","Conservation", "Gene", "SYMBOL",
                        "BIOTYPE","cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons"]
        prediction_index_map = {}
        for col in column_list:
                if self.variant.get_info_key_index(col) is not None:
                    prediction_index_map[col.lower()] = self.variant.get_info_key_index(col) 

        info_map = {}
        for csq_record in self.variant.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            allele = csq_record_list[prediction_index_map["allele"]]

            if allele not in info_map.keys():
                info_map[allele] = {"phenotype_assertions": [], "predicted_molecular_consequences": [], "prediction_results": []} 

            # parse and form phenotypes - adding phenotype from any of the csq record would be enough for adding only variant-linked phenotypes
            if not info_map[allele]["phenotype_assertions"]:
                phenotypes = csq_record_list[prediction_index_map["phenotypes"]].split("&") if "phenotypes" in prediction_index_map.keys() else []   
                for phenotype in phenotypes:
                    phenotype_assertions = self.create_allele_phenotype_assertion(phenotype) if phenotype else []
                    if (phenotype_assertions):
                        info_map[allele]["phenotype_assertions"].append(phenotype_assertions)
            
            # parse and form molecular consequences
            predicted_molecular_consequences = self.create_allele_predicted_molecular_consequence(csq_record_list, prediction_index_map)
            if (predicted_molecular_consequences):
                info_map[allele]["predicted_molecular_consequences"].append(predicted_molecular_consequences)
            
            # parse and form prediction results
            current_prediction_results = info_map[allele]["prediction_results"] 
            info_map[allele]["prediction_results"] += self.create_allele_prediction_results(current_prediction_results, csq_record_list, prediction_index_map)
        return info_map
     
    def create_allele_prediction_results(self, current_prediction_results: Mapping, csq_record: List, prediction_index_map: dict) -> list:
        prediction_results = []
        if "cadd_phred" in prediction_index_map.keys():
            if not self.prediction_result_already_exists(current_prediction_results, "CADD"):
                cadd_prediction_result = {
                        "score": csq_record[prediction_index_map["cadd_phred"]] ,
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": "CADD"
                        }

                } if csq_record[prediction_index_map["cadd_phred"]] else None
                if cadd_prediction_result:
                    prediction_results.append(cadd_prediction_result)

        
        return prediction_results
    
    def prediction_result_already_exists(self, current_prediction_results: Mapping, tool: str) -> bool:
        for prediction_result in current_prediction_results:
            if prediction_result["analysis_method"]["tool"] == tool:
                return True

        return False
    
    def create_allele_predicted_molecular_consequence(self, csq_record: List, prediction_index_map: dict) -> Mapping:
        consequences_list = []
        if "consequence" in prediction_index_map.keys():
            for cons in csq_record[prediction_index_map["consequence"]].split("&"):
                if cons in ["downstream_gene_variant", "upstream_gene_variant", "intergenic_variant", "regulatory_region_variant"]:
                    consequences_list = []
                    break
                consequences_list.append(
                    {
                        "value": cons
                    }
                )

        prediction_results = []
        if "sift" in prediction_index_map.keys():
            sift_score = csq_record[prediction_index_map["sift"]]
            if sift_score:
                (label, score) = self.format_sift_polyphen_output(sift_score)
                if label is not None and score is not None:
                    classification = {
                        "label": label,
                        "definition": ""
                    }
                    sift_prediction_result = {
                            "classification": classification,
                            "score": score,
                            "analysis_method": {
                                "tool": "SIFT",
                                "qualifier": "SIFT"
                            }
                        }
                    prediction_results.append(sift_prediction_result)

        if "polyphen" in prediction_index_map.keys():
            polyphen_score = csq_record[prediction_index_map["polyphen"]]
            if polyphen_score:
                (label, score) = self.format_sift_polyphen_output(polyphen_score)
                if label is not None and score is not None:
                    classification = {
                        "label": label,
                        "definition": ""
                    }
                    polyphen_prediction_result = {
                            "classification": classification,
                            "score": score,
                            "analysis_method": {
                                "tool": "PolyPhen",
                                "qualifier": "PolyPhen"
                            }
                        }
                    prediction_results.append(polyphen_prediction_result)
        
        cdna_location = cds_location = protein_location = None

        ###parse cdna location
        cdna_position = csq_record[prediction_index_map["cdna_position"]]
        if cdna_position:
            cdna_position_list  = cdna_position.split("-")
            cdna_start = cdna_position_list[0]
            cdna_end = cdna_start if len(cdna_position_list) < 2 else cdna_position_list[1]
            cdna_length = int(cdna_end) - int(cdna_start) + 1
            ref_sequence = self.minimise_allele(self.reference_sequence)
            alt_sequence = csq_record[prediction_index_map["allele"]]
            cdna_location = {
                "start": cdna_start,
                "end": cdna_end, 
                "length": cdna_length, 
                "ref_sequence": ref_sequence,
                "alt_sequence": alt_sequence
            }

        ###parse cds location
        cds_position = csq_record[prediction_index_map["cds_position"]]
        codons = csq_record[prediction_index_map["codons"]]
        if cds_position:
            cds_position_list  = cds_position.split("-")
            cds_start = cds_position_list[0]
            cds_end = cds_start if len(cds_position_list) < 2 else cds_position_list[1]
            cds_length = int(cds_end) - int(cds_start) + 1
            ref_sequence = codons.split("/")[0]
            alt_sequence = codons.split("/")[1]
            cds_location = {
                "start": cds_start,
                "end": cds_end, 
                "length": cds_length, 
                "ref_sequence": ref_sequence,
                "alt_sequence": alt_sequence
            }

        ###parse protein location
        protein_position = csq_record[prediction_index_map["protein_position"]]
        amino_acids = csq_record[prediction_index_map["amino_acids"]]
        if protein_position:
            protein_position_list  = protein_position.split("-")
            protein_start = protein_position_list[0]
            protein_end = protein_start if len(protein_position_list) < 2 else protein_position_list[1]
            protein_length = int(protein_end) - int(protein_start) + 1
            ref_sequence = amino_acids.split("/")[0]
            alt_sequence = amino_acids.split("/")[1]
            protein_location = {
                "start": protein_start,
                "end": protein_end, 
                "length": protein_length, 
                "ref_sequence": ref_sequence,
                "alt_sequence": alt_sequence
            }
        
        if consequences_list:
            return {
                "allele_name": csq_record[prediction_index_map["allele"]],
                "stable_id": csq_record[prediction_index_map["feature"]],
                "feature_type": {  
                    "value": csq_record[prediction_index_map["feature_type"]]     
                } ,
                "consequences": consequences_list,
                "gene_stable_id" : csq_record[prediction_index_map["gene"]], 
                "gene_symbol": csq_record[prediction_index_map["symbol"]], 
                "transcript_biotype": csq_record[prediction_index_map["biotype"]], 
                "protein_stable_id": "protein_id_placeholder",
                "prediction_results": prediction_results,
                "cdna_location": cdna_location,
                "cds_location": cds_location,
                "protein_location": protein_location
            }
            

    
    def format_sift_polyphen_output(self, output: str) -> tuple:
        try:
            (result, score) = re.split(r"[()]", output)[:2]
        except:
            return (None, None)

        if result not in [
            'probably damaging',
            'possibly damaging',
            'benign',
            'unknown',
            'tolerated',
            'deleterious',
            'tolerated - low confidence',
            'deleterious - low confidence',
        ]:
            result = None

        try:
            score = float(score)
        except:
            # need to log something here
            score = None

        return (result, score)
    
    def create_allele_phenotype_assertion(self, phenotype: str) -> Mapping:
        splits = phenotype.split("+")
        if len(splits) != 3 or splits[2].startswith("ENS"):
            return None

        phenotype_name,source,feature = splits
        evidence_list = []
        if re.search("^ENS.*G\d+", feature):
            feature_type = "Gene"
        elif re.search("^ENS.*T\d+", feature):
            feature_type = "Transcript"
        elif re.search("^rs", feature):
            feature_type = "Variant"  
        else:
            feature_type = None
        
        if phenotype:
            return {
                "feature": feature,
                "feature_type": {
                    "value": feature_type   

                } ,
                "phenotype": {
                    "name": phenotype_name
                },
                "evidence": evidence_list
            }

    
    def minimise_allele(self, alt: str) -> str:
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele_string = alt
        if len(alt) > len(self.reference_sequence):
            minimised_allele_string = alt[1:] 
        elif len(alt) < len(self.reference_sequence):
            minimised_allele_string = "-"
        elif alt == self.reference_sequence and len(self.alt) != len(self.reference_sequence):
            minimised_allele_string = alt[1:] if len(alt) > 1 else "-"
        return minimised_allele_string
    