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
        print("Hi")
        population_map = self.traverse_population_info()
        return population_map[min_alt] if min_alt in population_map else []

    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        column_list = ["Allele", "PHENOTYPES", "Feature_type", "Feature", "Consequence", "SIFT", "PolyPhen", "SPDI", "CADD_PHRED","Conservation"]
        prediction_index_map = {}
        for col in column_list:
                if self.variant.get_info_key_index(col) is not None:
                    prediction_index_map[col.lower()] = self.variant.get_info_key_index(col) 

        info_map = {}
        for csq_record in self.variant.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            if csq_record_list[prediction_index_map["allele"]] not in info_map.keys():
                info_map[csq_record_list[prediction_index_map["allele"]]] = {"phenotype_assertions": [], "predicted_molecular_consequences": [], "prediction_results": []} 

            phenotypes = csq_record_list[prediction_index_map["phenotypes"]].split("&") if "phenotypes" in prediction_index_map.keys() else []   
            for phenotype in phenotypes:
                phenotype_assertions = self.create_allele_phenotype_assertion(phenotype) if phenotype else []
                if (phenotype_assertions):
                    info_map[csq_record_list[prediction_index_map["allele"]]]["phenotype_assertions"].append(phenotype_assertions)
            predicted_molecular_consequences = self.create_allele_predicted_molecular_consequence(csq_record_list, prediction_index_map)
            if (predicted_molecular_consequences):
                info_map[csq_record_list[prediction_index_map["allele"]]]["predicted_molecular_consequences"].append(predicted_molecular_consequences)
            current_prediction_results = info_map[csq_record_list[prediction_index_map["allele"]]]["prediction_results"] 
            info_map[csq_record_list[prediction_index_map["allele"]]]["prediction_results"] += self.create_allele_prediction_results(current_prediction_results, csq_record_list, prediction_index_map)
        return info_map
    
    def traverse_population_info(self) -> Mapping:
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory,'populations.json')) as pop_file:
            pop_names = json.load(pop_file)
    
        kg_column_list = ["AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF"]

        gnomADe_column_map = { 
                       "ALL": ["gnomAD_exomes_AF","gnomAD_exomes_AC","gnomAD_exomes_AN"],
                       "afr": ["gnomAD_exomes_AF_afr","gnomAD_exomes_AC_afr","gnomAD_exomes_AN_afr"],
                       "amr": ["gnomAD_exomes_AF_amr","gnomAD_exomes_AC_amr","gnomAD_exomes_AN_amr"],
                       "asj": ["gnomAD_exomes_AF_asj","gnomAD_exomes_AC_asj","gnomAD_exomes_AN_asj"],
                       "eas": ["gnomAD_exomes_AF_eas","gnomAD_exomes_AC_eas","gnomAD_exomes_AN_eas"],
                       "fin": ["gnomAD_exomes_AF_fin","gnomAD_exomes_AC_fin","gnomAD_exomes_AN_fin"],
                       "nfe": ["gnomAD_exomes_AF_nfe","gnomAD_exomes_AC_nfe","gnomAD_exomes_AN_nfe"],
                       "oth": ["gnomAD_exomes_AF_oth","gnomAD_exomes_AC_oth","gnomAD_exomes_AN_oth"],
                       "sas": ["gnomAD_exomes_AF_sas","gnomAD_exomes_AC_sas","gnomAD_exomes_AN_sas"] 
                       }
        
        gnomADg_column_map = { 
                       "ALL": ["gnomAD_genomes_AF","gnomAD_genomes_AC","gnomAD_genomes_AN"],
                       "afr": ["gnomAD_genomes_AF_afr","gnomAD_genomes_AC_afr","gnomAD_genomes_AN_afr"],
                       "ami": ["gnomAD_genomes_AF_ami","gnomAD_genomes_AC_ami","gnomAD_genomes_AN_ami"],
                       "amr": ["gnomAD_genomes_AF_amr","gnomAD_genomes_AC_amr","gnomAD_genomes_AN_amr"],
                       "asj": ["gnomAD_genomes_AF_asj","gnomAD_genomes_AC_asj","gnomAD_genomes_AN_asj"],
                       "eas": ["gnomAD_genomes_AF_eas","gnomAD_genomes_AC_eas,gnomAD_genomes_AN_eas"],
                       "fin": ["gnomAD_genomes_AF_fin","gnomAD_genomes_AC_fin","gnomAD_genomes_AN_fin"],
                       "mid": ["gnomAD_genomes_AF_mid","gnomAD_genomes_AC_mid","gnomAD_genomes_AN_mid"],
                       "nfe": ["gnomAD_genomes_AF_nfe","gnomAD_genomes_AC_nfe","gnomAD_genomes_AN_nfe"],
                       "oth": ["gnomAD_genomes_AF_oth","gnomAD_genomes_AC_oth","gnomAD_genomes_AN_oth"],
                       "sas": ["gnomAD_genomes_AF_sas","gnomAD_genomes_AC_sas","gnomAD_genomes_AN_sas"]
        }
        
        population_frequency_map = {}
        frequency_stats = {}
        
        for csq_record in self.variant.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            allele_index = self.variant.get_info_key_index("Allele") 
            if allele_index is not None and csq_record_list[allele_index] is not None:
                if csq_record_list[allele_index] not in population_frequency_map.keys():
                    population_frequency_map[csq_record_list[allele_index]] = []
                for col in kg_column_list:
                    col_index = self.variant.get_info_key_index(col)
                    if col_index is not None and csq_record_list[col_index] is not None:
                            if pop_names[col] not in frequency_stats.keys():
                                frequency_stats[pop_names[col]] = {}
                            frequency_stats[pop_names[col]][csq_record_list[allele_index]] = csq_record_list[col_index]  
                            print(csq_record_list[col_index])  
                            if csq_record_list[col_index]: 
                                print(csq_record_list[col_index])  
                                population_frequency = {
                                                "population": pop_names[col],
                                                "allele_frequency": csq_record_list[col_index],
                                                "is_minor_allele": False,
                                                "is_hpmaf": False
                                            }
                                if population_frequency_map[csq_record_list[allele_index]]
                                population_frequency_map[csq_record_list[allele_index]].append(population_frequency)

                for key,gnomADe_column_list in gnomADe_column_map.items():
                    for col in gnomADe_column_list:
                        population_frequency_key = col.split("_")[2]
                        allele_count = allele_number = allele_frequency = None
                        col_index = self.variant.get_info_key_index(col)
                        if col_index and csq_record_list[col_index]:
                            if population_frequency_key == "AF":
                                allele_frequency = csq_record_list[col_index] 
                            elif population_frequency_key == "AN":
                                allele_number = csq_record_list[col_index]
                            elif population_frequency_key == "AC":
                                allele_count = csq_record_list[col_index]
                            else:
                                raise Exception('gnomAD exomes column is not recognised')

                    if allele_frequency is not None:
                        if pop_names["gnomAD_exomes"][key] not in frequency_stats.keys():
                            frequency_stats[pop_names["gnomAD_exomes"][key]] = {}
                        frequency_stats[pop_names["gnomAD_exomes"][key]][csq_record_list[allele_index]] = allele_frequency     
                        population_frequency = {
                                        "population": pop_names["gnomAD_exomes"][key],
                                        "allele_frequency": allele_frequency,
                                        "allele_count": allele_count,
                                        "allele_number": allele_number,
                                        "is_minor_allele": False,
                                        "is_hpmaf": False
                                    }
                        population_frequency_map[csq_record_list[allele_index]].append(population_frequency)
            
                for key,gnomADg_column_list in gnomADg_column_map.items():
                    for col in gnomADg_column_list:
                        allele_count = allele_number = allele_frequency = None
                        population_frequency_key = col.split("_")[2]
                        col_index = self.variant.get_info_key_index(col)
                        if col_index is not None:
                            if population_frequency_key == "AF":
                                allele_frequency = csq_record_list[col_index] or None
                            elif population_frequency_key == "AN":
                                allele_number = csq_record_list[col_index] or None
                            elif population_frequency_key == "AC":
                                allele_count = csq_record_list[col_index] or None
                            else:
                                raise Exception('gnomAD genomes column is not recognised')

                    if allele_frequency is not None:
                        if pop_names["gnomAD_genomes"][key] not in frequency_stats.keys():
                            frequency_stats[pop_names["gnomAD_genomes"][key]] = {}
                        frequency_stats[pop_names["gnomAD_genomes"][key]][csq_record_list[allele_index]] = allele_frequency    
                        population_frequency = {
                                        "population": pop_names["gnomAD_genomes"][key],
                                        "allele_frequency": allele_frequency,
                                        "allele_count": allele_count,
                                        "allele_number": allele_number,
                                        "is_minor_allele": False,
                                        "is_hpmaf": False
                                    }
                        population_frequency_map[csq_record_list[allele_index]].append(population_frequency)
        print(frequency_stats)
        return population_frequency_map
     
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
                consequences_list.append(
                    {
                        "accession_id": cons
                    }
                )

        prediction_results = []
        if "sift" in prediction_index_map.keys():
            sift_score = csq_record[prediction_index_map["sift"]]
            if sift_score:
                (result, score) = self.format_sift_polyphen_output(sift_score)
                if result is not None and score is not None:
                    sift_prediction_result = {
                            "result": result,
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
                (result, score) = self.format_sift_polyphen_output(polyphen_score)
                if result is not None and score is not None:
                    polyphen_prediction_result = {
                            "result": result,
                            "score": score,
                            "analysis_method": {
                                "tool": "PolyPhen",
                                "qualifier": "PolyPhen"
                            }
                        }
                    prediction_results.append(polyphen_prediction_result)
        if consequences_list:
            return {
                "allele_name": csq_record[prediction_index_map["allele"]],
                "feature_stable_id": csq_record[prediction_index_map["feature"]],
                "feature_type": {
                    "accession_id": csq_record[prediction_index_map["feature_type"]]               
                } ,
                "consequences": consequences_list,
                "prediction_results": prediction_results
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
        phenotype_name,source,feature = phenotype.split("+")
        evidence_list = []
        if re.search("^ENS.*G\d+", feature):
            feature_type = "Gene"
        elif re.search("^ENS.*T\d+", feature):
            feature_type = "Transcript"
        elif re.search("^rs", feature):
            feature_type = "Variant"  
        else:
            feature_type = None 
        # if source:
        #     evidence =  {
        #         "source": {
        #             "id": source,
        #             "name": source.replace("_"," ") 
        #         }
        #     }
        #     evidence_list.append(evidence)
        
        if phenotype:
            return {
                "feature": feature,
                "feature_type": {
                    "accession_id": feature_type                
                } ,
                "phenotype": {
                    "name": phenotype_name
                },
                "evidence": evidence_list
            }

    
    def minimise_allele(self, alt: str):
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele = alt
        if len(alt) > len(self.variant.ref):
            minimised_allele = alt[1:] 
        elif len(alt) < len(self.variant.ref):
            minimised_allele = "-"
        return minimised_allele
    
    # def format_frequency(self, raw_frequency_list: List) -> Mapping:
    #     freq_map = {}
    #     for freq in raw_frequency_list:
    #         key = freq.split(":")[0]
    #         freq_list = freq.split(":")[1].split(",")
    #         freq_map[key] = freq_list
    #     return freq_map
    
    # def get_population_allele_frequencies(self) -> List:
    #     frequency_map = {}
    #     if "FREQ" in self.variant.info:
    #         frequency_map = self.format_frequency(",".join(map(str,self.variant.info["FREQ"])).split("|"))
    #     population_allele_frequencies = []
    #     for key, pop_list in frequency_map.items():
    #         ## Adding only GnomAD population
    #         if key == "GnomAD":
    #             if pop_list[self.allele_index] not in ["None", "."]:
    #                 population_allele_frequencies.append({
    #                     "population": key,
    #                     "allele_frequency": pop_list[self.allele_index],
    #                     "is_minor_allele": False,
    #                     "is_hpmaf": False
    #                 })
    #     return population_allele_frequencies

  
