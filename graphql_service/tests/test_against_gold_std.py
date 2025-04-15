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

# -----------------------------------------------------------------------
# Setup

from string import Template
from ariadne import graphql
import pytest
from .test_utils import setup_test
import os
import json
import difflib

# -----------------------------------------------------------------------
# Define helper functions

# Generate test cases from the gold standard JSON file names
def get_test_case_ids():
    """Dynamically generate test cases from gold standard JSON files.
    
    Expected file naming: <variant_id>_broken.json
    Uses the same genome_id for all test cases.
    """
    gold_std_dir = "/app/graphql_service/tests/gold_std_query_results"
    genome_id = "a7335667-93e7-11ec-a39d-005056b38ce3"
    cases = []
    for filename in os.listdir(gold_std_dir):
        if filename.endswith(".json"):
            variant_id = filename.replace(".json", "")
            cases.append((variant_id, genome_id))
    return cases

TEST_CASES = get_test_case_ids()

# Set up in-memory schema and context
executable_schema, context = setup_test()

# Define master query template
def master_query_template():
    """Query template"""
    
    query = """{
        variant(
            by_id: { genome_id: "$genome_id", variant_id: "$variant_id" }
        ) {
            name
            alternative_names { __typename }
            primary_source { __typename }
            type
            allele_type { __typename }
            slice  { __typename }
            alleles { 
            name
            allele_sequence
            reference_sequence
            alternative_names { __typename }
            type
            allele_type { __typename }
            slice { __typename }
            phenotype_assertions {
                feature
                feature_type { __typename }
                phenotype {
                name
                source { __typename }
                ontology_terms { __typename }
                }
                evidence {
                source { __typename }
                assertion { __typename }
                }
            }
            prediction_results {
                score
                result
                classification { __typename }
                analysis_method { __typename }
            }
            population_frequencies {
                population_name
                allele_frequency
                allele_count
                allele_number
                is_minor_allele
                is_hpmaf
            }
            predicted_molecular_consequences {
                allele_name
                stable_id
                feature_type { __typename }
                consequences { __typename }
                prediction_results {
                score
                result
                classification { __typename }
                analysis_method { __typename }
                }
                gene_stable_id
                gene_symbol
                protein_stable_id
                transcript_biotype
                cdna_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
                cds_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
                protein_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
            }
            ensembl_website_display_data {
                count_transcript_consequences
                count_overlapped_genes
                count_regulatory_consequences
                count_variant_phenotypes
                count_gene_phenotypes
                representative_population_allele_frequency
            }
            }
            prediction_results {
            score
            result
            classification { __typename }
            analysis_method { __typename }
            }
            ensembl_website_display_data {
            count_citations
            }
        }
    }"""
    return Template(query)

# Define helper function to submit query and receive response
async def execute_query(genome_id, variant_id):
    """Execute the query with given parameters and return (query, success, result)."""
    template = master_query_template()
    query = template.substitute(genome_id=genome_id, variant_id=variant_id)
    query_data = {"query": query}
    success, result = await graphql(
        executable_schema, query_data, context_value=context(request={})
    )
    return query, success, result

# -----------------------------------------------------------------------
# Run tests

@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", TEST_CASES)
async def test_query_results_against_gold_std(variant_id, genome_id):
    """Test present query result for variant X against gold 
    standard query result for variant X."""

    query, success, result = await execute_query(genome_id, variant_id)
    assert success, f"[Variant Root] Query execution failed for variant {variant_id}. Query: {query}. Result: {result}"
    
    # Define the path to the gold standard JSON file
    gold_std_path = os.path.join("/app/graphql_service/tests/gold_std_query_results", f"{variant_id}.json")
    
    # Load the gold standard JSON file
    with open(gold_std_path, "r") as f:
        gold_standard = json.load(f)
    
    # Convert both JSON objects to strings for diffing, ensure sorting of keys
    result_str = json.dumps(result.get("data", {}), indent=4, sort_keys=True)
    gold_str = json.dumps(gold_standard.get("data", {}), indent=4, sort_keys=True)
    
    # Generate diff of the gold standard and actual variant
    # The difflib module will give a git-diff-like output, showing precisely where the 
    # differences are between the actual query result and the gold standard query result
    diff = ''.join(
        difflib.unified_diff(
            gold_str.splitlines(keepends=True),
            result_str.splitlines(keepends=True),
            fromfile="Gold Standard",
            tofile="Result"
        )
    )
    
    # Log the diff if there are any differences
    if result_str != gold_str: # if diff?
        print(f"[Variant Root] Diff for variant {variant_id}:\n{diff}")
    
    # Compare the result with the gold standard - assert that they're identical
    assert result_str == gold_str