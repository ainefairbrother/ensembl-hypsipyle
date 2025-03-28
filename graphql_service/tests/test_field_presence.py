"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from string import Template
from ariadne import graphql
import pytest
from .test_utils import setup_test

# -----------------------------------------------------------------------
# General setup: define test cases and set-up schema/context

GENOME_ID="a7335667-93e7-11ec-a39d-005056b38ce3"

VARIANT_TEST_CASES = [
    ("1:10007:rs1639538116", GENOME_ID),
    ("20:35174034:rs2069945", GENOME_ID),
    # ("1:2193890:rs1688548931", GENOME_ID), # revisit this one
    ("1:2193889:rs1688548931", GENOME_ID),
    ("17:39531219:rs760014508", GENOME_ID),
    ("1:924510:rs1405511870", GENOME_ID), # wasn't found, so looked up on dbsnp and changed 924511>924510
    ("19:375858:rs1973467575", GENOME_ID), # wasn't found, so looked up on dbsnp and changed 375859>375858
    ("19:474690:rs569586139", GENOME_ID),
    ("19:375786:rs747148860", GENOME_ID),
    ("19:685717:rs1452804066", GENOME_ID),
    ("19:695430:rs750573136", GENOME_ID),
    ("19:111026:rs1268711754", GENOME_ID),
    ("14:91234215:rs1555409827", GENOME_ID),
    ("1:999609:rs1463383567", GENOME_ID),
    ("1:964529:rs1642816219", GENOME_ID)
]

# Set up in-memory schema and context
executable_schema, context = setup_test()

# -----------------------------------------------------------------------
# Define helper functions

# Define query template with basic structure that queries the variant table
def query_template(additional_fields=""):
    """Query template"""
    
    query = """{
        variant(
            by_id: { genome_id: "$genome_id", variant_id: "$variant_id" }
        ) {
            %s
        }
    }""" % additional_fields
    return Template(query)

async def execute_query(genome_id, variant_id, additional_fields):
    """Execute the query with given parameters and return (query, success, result)."""
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id, variant_id=variant_id)
    query_data = {"query": query}
    success, result = await graphql(
        executable_schema, query_data, context_value=context(request={})
    )
    return query, success, result

# -----------------------------------------------------------------------
# Test:: Root Variant Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_root_fields_present(variant_id, genome_id):
    """Test that all top-level Variant fields are returned."""
    
    additional_fields = """
        name
        alternative_names { __typename }
        primary_source { __typename }
        type
        allele_type { __typename }
        slice { __typename }
        alleles { __typename }
        prediction_results { __typename }
        ensembl_website_display_data { __typename }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Variant Root] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    variant = result.get("data", {}).get("variant")
    assert variant is not None, f"[Variant Root] Variant is None for variant {variant_id}. Errors: {result.get('errors', 'No error info')}. Query: {query}"
    
    expected_fields = [
        "name", "alternative_names", "primary_source", "type", "allele_type", 
        "slice", "alleles", "prediction_results", "ensembl_website_display_data"
    ]
    missing_fields = [field for field in expected_fields if field not in variant]
    assert not missing_fields, f"[Variant Root] Missing fields in variant {variant_id}: {missing_fields}. Query: {query}. Result: {result}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_fields_present(variant_id, genome_id):
    """Test that all allele-level fields are returned for each variant."""
    
    additional_fields = """
        alleles {
            name
            allele_sequence
            reference_sequence
            alternative_names { __typename }
            type
            allele_type { __typename }
            slice { __typename }
            phenotype_assertions { __typename }
            prediction_results { __typename }
            population_frequencies { __typename }
            predicted_molecular_consequences { __typename }
            ensembl_website_display_data { __typename }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    assert alleles is not None, f"[Allele] alleles field is missing for variant {variant_id}. Query: {query}. Result: {result}"
    
    expected_fields = [
        "name", "allele_sequence", "reference_sequence", "alternative_names", "type",
        "allele_type", "slice", "phenotype_assertions", "prediction_results",
        "population_frequencies", "predicted_molecular_consequences", "ensembl_website_display_data"
    ]
    for index, allele in enumerate(alleles):
        missing = [field for field in expected_fields if field not in allele]
        assert not missing, f"[Allele] Missing fields in allele {index} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Phenotype Assertions Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_phenotype_assertions_present(variant_id, genome_id):
    """Test that phenotype_assertions in each allele are either empty or contain the expected subfields."""
    
    additional_fields = """
        alleles {
            phenotype_assertions {
                feature
                feature_type { __typename }
                phenotype { __typename }
                evidence { __typename }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Phenotype Assertions] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles", [])    
    expected_fields = [
        "feature", "feature_type", "phenotype", "evidence"
    ]
    for idx, allele in enumerate(alleles):
        phenotype_assertions = allele.get("phenotype_assertions", [])
        # If phenotype_assertions isn't empty, it should have expected fields
        if phenotype_assertions:
            for phenotype_assertion in phenotype_assertions:
                missing = [field for field in expected_fields if field not in phenotype_assertion]
                assert not missing, f"[Allele: Phenotype Assertions] Missing fields in allele {idx} phenotype_assertions for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
        else:
            # Accept empty phenotype_assertions
            assert phenotype_assertions == [], f"[Allele: Phenotype Assertions] Expected empty phenotype_assertions for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant Allele Phenotype Assertions Phenotype Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_phenotype_assertions_phenotype_present(variant_id, genome_id):
    """Test that each allele's phenotype_assertions, if present, include a phenotype field with relevant subfields."""
    
    additional_fields = """
        alleles {
            phenotype_assertions {
                phenotype { 
                    name
                    source { __typename }
                    ontology_terms { __typename }
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Phenotype Assertions: Phenotype] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    expected_fields = [
        "name", "source", "ontology_terms"
    ]    
    
    # Iterate over each allele
    for idx, allele in enumerate(alleles):
        phenotype_assertions = allele.get("phenotype_assertions", [])
        
        # If phenotype_assertions isn't empty, it should have expected fields
        if phenotype_assertions:
            
            # Iterate over each phenotype assertion
            for assertion in phenotype_assertions:
                
                # Get phenotype
                phenotype = assertion.get("phenotype")
                
                # Phenotype should not be none and should have expected fields
                assert phenotype is not None, (f"[Allele: Phenotype Assertions: Phenotype] Allele {idx} in variant {variant_id} has phenotype_assertions but is missing the 'phenotype' field. Query: {query}")
                
                # Check for missing fields in phenotype
                missing = [field for field in expected_fields if field not in phenotype]
                assert not missing, f"[Allele: Phenotype Assertions: Phenotype] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
        else:
            # Accept empty phenotype_assertions
            assert phenotype_assertions == [], f"[Allele: Phenotype Assertions: Phenotype] Expected empty phenotype_assertions for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Phenotype Assertions Evidence Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_phenotype_assertions_evidence_present(variant_id, genome_id):
    """Test that phenotype_assertions in each allele is either empty or contains an evidence subfield with expected subfields."""
    
    additional_fields = """
        alleles {
            phenotype_assertions {
                evidence{
                    source{ __typename }
                    assertion{ __typename }
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Phenotype Assertions: Evidence] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    expected_fields = [
        "source", "assertion"
    ]
    
    # Iterate over each allele
    for idx, allele in enumerate(alleles):
        phenotype_assertions = allele.get("phenotype_assertions", [])
        
        # If phenotype_assertions isn't empty, it should have expected fields
        if phenotype_assertions:
            
            # Iterate over each phenotype assertion
            for assertion in phenotype_assertions:
                
                # Get evidence
                evidence = assertion.get("evidence")
                
                # If evidence isn't empty, it should have expected fields
                if evidence:
                
                    # Check for missing fields in evidence
                    missing = [field for field in expected_fields if field not in evidence]
                    assert not missing, f"[Allele: Phenotype Assertions: Evidence] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
                
                # Accept empty evidence
                else:
                    assert evidence == []
        else:
            # Accept empty phenotype_assertions
            assert phenotype_assertions == [], f"[Allele: Phenotype Assertions: Evidence] Expected empty phenotype_assertions for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Prediction Results Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_prediction_results_present(variant_id, genome_id):
    """Test that prediction_results in each allele is either empty or contains the expected subfields."""
    
    additional_fields = """
        alleles {
            prediction_results{
                score
                result
                classification{ __typename }
                analysis_method{ __typename }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Prediction Results] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles", [])
    expected_fields = [
        "score", "result", "classification", "analysis_method"
    ]
    for idx, allele in enumerate(alleles):
        prediction_results = allele.get("prediction_results", [])
        if prediction_results:
            for prediction in prediction_results:
                missing = [field for field in expected_fields if field not in prediction]
                assert not missing, f"[Allele: Prediction Results] Missing fields in allele {idx} prediction_results for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
        else:
            assert prediction_results == [], f"[Allele: Prediction Results] Expected empty prediction_results for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Population Frequencies Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_population_frequencies_present(variant_id, genome_id):
    """Test that population_frequencies in each allele is either empty or contains the expected subfields."""
    
    additional_fields = """
        alleles {
            population_frequencies{
                population_name
                allele_frequency
                allele_count
                allele_number
                is_minor_allele
                is_hpmaf
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Population Frequencies] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles", [])
    expected_fields = [
        "population_name", "allele_frequency", "allele_count", "allele_number", "is_minor_allele", "is_hpmaf"
    ]
    for idx, allele in enumerate(alleles):
        population_frequencies = allele.get("population_frequencies", [])
        # If population_frequencies isn't empty, it should have expected fields
        if population_frequencies:
            for population_frequency in population_frequencies:
                missing = [field for field in expected_fields if field not in population_frequency]
                assert not missing, f"[Allele: Population Frequencies] Missing fields in allele {idx} population_frequencies for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
        else:
            # Accept empty population_frequencies
            assert population_frequencies == [], f"[Allele: Population Frequencies] Expected empty phenotype_aspopulation_frequenciessertions for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Predicted Molecular Consequences Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_predicted_molecular_consequences_present(variant_id, genome_id):
    """Test that predicted_molecular_consequences in each allele is either empty or contains the expected subfields."""
    
    additional_fields = """
        alleles {
            predicted_molecular_consequences{
                allele_name
                stable_id
                feature_type{ __typename }
                consequences{ __typename }
                prediction_results{ __typename }
                gene_stable_id
                gene_symbol
                protein_stable_id
                transcript_biotype
                cdna_location{ __typename }
                cds_location{ __typename }
                protein_location{ __typename }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Template] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [ 
        "allele_name", "stable_id", "feature_type", "consequences", 
        "prediction_results", "gene_stable_id", "gene_symbol", 
        "protein_stable_id", "transcript_biotype", "cdna_location", 
        "cds_location", "protein_location"
    ]
    
    for idx, allele in enumerate(alleles):
        predicted_molecular_consequences = allele.get("predicted_molecular_consequences", [])
        # If predicted_molecular_consequences isn't empty, it should have expected fields
        if predicted_molecular_consequences:
            for consequence in predicted_molecular_consequences:
                missing = [field for field in expected_fields if field not in consequence]
                assert not missing, f"[Allele: Predicted Molecular Consequences] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
        else:
            # Accept empty predicted_molecular_consequences
            assert predicted_molecular_consequences == [], f"[Allele: Predicted Molecular Consequences] Expected empty predicted_molecular_consequences field for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Predicted Molecular Consequences: Prediction Results Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_predicted_molecular_consequences_predicted_results_present(variant_id, genome_id):
    """Test that prediction_results in predicted_molecular_consequences in each allele contains the expected subfields."""
    
    additional_fields = """
        alleles {
            predicted_molecular_consequences{
                prediction_results{ 
                    score
                    result
                    classification{ __typename }
                    analysis_method{ __typename }
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Predicted Molecular Consequences: Predicted Results] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [
        "score", "result", "classification", "analysis_method"
    ]
    
    for idx, allele in enumerate(alleles):
        predicted_molecular_consequences = allele.get("predicted_molecular_consequences", [])
        # If predicted_molecular_consequences isn't empty, it should have expected fields
        if predicted_molecular_consequences:
            for predicted_molecular_consequence in predicted_molecular_consequences:
                prediction_results = predicted_molecular_consequence.get("prediction_results", [])
                # If prediction_results isn't empty, it should have expected fields
                if prediction_results:
                    for prediction_result in prediction_results:
                        missing = [field for field in expected_fields if field not in prediction_result]
                        assert not missing, f"[Allele: Predicted Molecular Consequences: Predicted Results] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
                # Accept empty prediction_result
                else:
                    assert prediction_results == []
        else:
            # Accept empty predicted_molecular_consequences
            assert predicted_molecular_consequences == [], f"[Allele: Predicted Molecular Consequences: Predicted Results] Expected empty predicted_molecular_consequences field for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Predicted Molecular Consequences: cDNA Location Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_predicted_molecular_consequences_cdna_location_present(variant_id, genome_id):
    """Test that cdna_location in predicted_molecular_consequences in each allele contains the expected subfields."""
    
    additional_fields = """
        alleles {
            predicted_molecular_consequences{
                cdna_location{
                    relation { __typename }
                    start
                    end
                    length
                    percentage_overlap
                    ref_sequence
                    alt_sequence
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Predicted Molecular Consequences: cDNA Location] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [
        "relation", "start", "end", "length", "percentage_overlap", "ref_sequence", "alt_sequence"
    ]
    
    for idx, allele in enumerate(alleles):
        predicted_molecular_consequences = allele.get("predicted_molecular_consequences", [])
        # If predicted_molecular_consequences isn't empty, it should have expected fields
        if predicted_molecular_consequences:
            for predicted_molecular_consequence in predicted_molecular_consequences:
                cdna_location = predicted_molecular_consequence.get("cdna_location")
                # If cdna_location isn't empty, it should have expected fields
                if cdna_location:
                    missing = [field for field in expected_fields if field not in cdna_location]
                    assert not missing, f"[Allele: Predicted Molecular Consequences: cDNA Location] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
                # Accept cdna_location as None
                else:
                    assert cdna_location is None
        else:
            # Accept empty predicted_molecular_consequences
            assert predicted_molecular_consequences == [], f"[Allele: Predicted Molecular Consequences: cDNA Location] Expected empty predicted_molecular_consequences field for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Predicted Molecular Consequences: CDS Location Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_predicted_molecular_consequences_cds_location_present(variant_id, genome_id):
    """Test that cds_location in predicted_molecular_consequences in each allele contains the expected subfields."""
    
    additional_fields = """
        alleles {
            predicted_molecular_consequences{
                cds_location{
                    relation { __typename }
                    start
                    end
                    length
                    percentage_overlap
                    ref_sequence
                    alt_sequence
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Predicted Molecular Consequences: CDS Location] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [
        "relation", "start", "end", "length", "percentage_overlap", "ref_sequence", "alt_sequence"
    ]
    
    for idx, allele in enumerate(alleles):
        predicted_molecular_consequences = allele.get("predicted_molecular_consequences", [])
        # If predicted_molecular_consequences isn't empty, it should have expected fields
        if predicted_molecular_consequences:
            for predicted_molecular_consequence in predicted_molecular_consequences:
                cds_location = predicted_molecular_consequence.get("cds_location")
                # If cds_location isn't empty, it should have expected fields
                if cds_location:
                    missing = [field for field in expected_fields if field not in cds_location]
                    assert not missing, f"[Allele: Predicted Molecular Consequences: CDS Location] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
                # Accept cds_location as None
                else:
                    assert cds_location is None
        else:
            # Accept empty predicted_molecular_consequences
            assert predicted_molecular_consequences == [], f"[Allele: Predicted Molecular Consequences: CDS Location] Expected empty predicted_molecular_consequences field for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Predicted Molecular Consequences: Protein Location Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_predicted_molecular_consequences_protein_location_present(variant_id, genome_id):
    """Test that protein_location in predicted_molecular_consequences in each allele contains the expected subfields."""
    
    additional_fields = """
        alleles {
            predicted_molecular_consequences{
                protein_location{
                    relation { __typename }
                    start
                    end
                    length
                    percentage_overlap
                    ref_sequence
                    alt_sequence
                }
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Predicted Molecular Consequences: Protein Location] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [
        "relation", "start", "end", "length", "percentage_overlap", "ref_sequence", "alt_sequence"
    ]
    
    for idx, allele in enumerate(alleles):
        predicted_molecular_consequences = allele.get("predicted_molecular_consequences", [])
        # If predicted_molecular_consequences isn't empty, it should have expected fields
        if predicted_molecular_consequences:
            for predicted_molecular_consequence in predicted_molecular_consequences:
                protein_location = predicted_molecular_consequence.get("protein_location")
                # If protein_location isn't empty, it should have expected fields
                if protein_location:
                    missing = [field for field in expected_fields if field not in protein_location]
                    assert not missing, f"[Allele: Predicted Molecular Consequences: Protein Location] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"
                # Accept protein_location as None
                else:
                    assert protein_location is None
        else:
            # Accept empty predicted_molecular_consequences
            assert predicted_molecular_consequences == [], f"[Allele: Predicted Molecular Consequences: Protein Location] Expected empty predicted_molecular_consequences field for allele {idx} in variant {variant_id}. Query: {query}"

# -----------------------------------------------------------------------
# Test:: Variant: Alleles: Ensembl Website Display Data Fields
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_variant_allele_ensembl_website_display_data_present(variant_id, genome_id):
    """Test that ensembl_website_display_data in each allele contains the expected subfields."""
    
    additional_fields = """
        alleles {
            ensembl_website_display_data{
                count_transcript_consequences
                count_overlapped_genes
                count_regulatory_consequences
                count_variant_phenotypes
                count_gene_phenotypes
                representative_population_allele_frequency
            }
        }
    """
    query, success, result = await execute_query(genome_id, variant_id, additional_fields)
    assert success, f"[Allele: Ensembl Website Display Data] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    alleles = result.get("data", {}).get("variant").get("alleles")
    
    expected_fields = [
        "count_transcript_consequences", "count_overlapped_genes", "count_regulatory_consequences",
        "count_variant_phenotypes", "count_gene_phenotypes", "representative_population_allele_frequency"
    ]
    
    for idx, allele in enumerate(alleles):
        ensembl_website_display_data = allele.get("ensembl_website_display_data")
        missing = [field for field in expected_fields if field not in ensembl_website_display_data]
        assert not missing, f"[Allele: Ensembl Website Display Data] Missing fields in allele {idx} for variant {variant_id}: {missing}. Query: {query}. Result: {result}"











# # -----------------------------------------------------------------------
# # Test: Template 
# @pytest.mark.asyncio
# @pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
# async def template(variant_id, genome_id):
#     """Test that all allele-level fields are returned for each variant."""
    
#     additional_fields = """
#         template {

#         }
#     """
#     query, success, result = await execute_query(genome_id, variant_id, additional_fields)
#     assert success, f"[Template] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
#     template = result.get("data", {}).get("variant")#.get...
#     assert template is not None, f"[Template] 'template' field is missing for variant {variant_id}. Query: {query}. Result: {result}"
    
#     expected_fields = [

#     ]