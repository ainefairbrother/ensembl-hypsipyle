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

import logging
from string import Template
from ariadne import graphql
import pytest
from .test_utils import setup_test

# define test parameters
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

# run setup fn.
executable_schema, context = setup_test()

# define query template with basic structure that queries the variant table with
# a genome_id and variant_id
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

# -----------------------------------------------------------------------

@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_basic(variant_id, genome_id):
    """Test that the very simplest query is working - are name and type
    being returned in response to a variant-based query."""
    
    additional_fields = """
        name
        type
    """
    
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id,variant_id=variant_id)
    query_data = {"query": query}
    (success, result) = await graphql(executable_schema, query_data, context_value=context(request={}))
    assert success, f"Query execution failed for variant {variant_id}. Result: {result}"
    
    # get fields and check that they exist
    variant = result.get("data", {}).get("variant")
    
    # variant should not be none, a name and type should be returned
    assert variant is not None, f"Variant is None. Errors: {result.get('errors', 'No error info')}"
    
    # given that variant isn't none, check that name and type exist
    assert all(field in variant for field in ["name", "type"]), "Missing 'name' or 'type' field in query result."
    
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_alternative_names(variant_id, genome_id):
    """Test for the presence and validity of Variant fields.
    Specifically tests the alternative_names field which can
    either be an empty array, or else it should contain a number 
    of fields. This test aditionally checks for the presence of 
    all sub-fields in the assignment_method and source fields."""
    
    additional_fields = """
        alternative_names{
            accession_id
            name
            description
            assignment_method{
                type
                description
            }
            url
            source{
                id
                name
                description
                url
                release
            }
        }
    """
    
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id,variant_id=variant_id)
    query_data = {"query": query}
    (success, result) = await graphql(executable_schema, query_data, context_value=context(request={}))
    assert success, f"Query execution failed for variant {variant_id}. Result: {result}"
    
    # test for existance of the alternative_names field
    variant = result.get("data", {}).get("variant")
    assert "alternative_names" in variant, "Missing 'alternative_names' field in query result."
    
    # get the alternative_names field
    alternative_names = variant.get("alternative_names")

    # if alternative names isn't empty, check for fields and subfields, otherwise accept empty as valid return
    if alternative_names:
        assert all(field in alternative_names for field in ['accession_id', 'name', 'description', 'assignment_method', 'url', 'source']), "One or more expected fields are missing in variant.alternative_names"
        
        # then check that the assignment_method and source subfields have their fields
        # assignment_method
        assignment_method = alternative_names.get("assignment_method")
        assert all(field in assignment_method for field in ["type", "description"]), "One or more expected fields are missing in variant.alternative_names.assignment_method."
        # source
        source = alternative_names.get("source")
        assert all(field in source for field in ["id", "name", "description", "url", "release"]), "One or more expected fields are missing in variant.alternative_names.source."
        
    else:
        assert alternative_names == []
        
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_primary_source(variant_id, genome_id):
    """DOCSTRING"""
    
    additional_fields = """
        primary_source{
            accession_id
            name
            assignment_method{
                type
            }
            source{
                id
                name
            }
        }
    """
    
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id,variant_id=variant_id)
    query_data = {"query": query}
    (success, result) = await graphql(executable_schema, query_data, context_value=context(request={}))
    assert success, f"Query execution failed for variant {variant_id}. Result: {result}"
    
    # Make sure primary_source field is returned
    variant = result.get("data", {}).get("variant")
    assert variant is not None, f"Variant is None. Errors: {result.get('errors', 'No error info')}"
    assert "primary_source" in variant, "Missing 'primary_source' field in query result."
    
    # Make sure mandatory fields are present in primary_source
    primary_source = variant.get("primary_source")
    assert all(field in primary_source for field in ['accession_id', 'name', 'assignment_method', 'source']), "One or more expected fields are missing in variant.primary_source."
    
    # Make sure mandatory fields are present in primary_source subfields
    # assignment_method
    assignment_method = primary_source.get("assignment_method")
    assert "type" in assignment_method, "Missing 'type' field in variant.primary_source.assignment_method query result."
    # source
    source = primary_source.get("source")
    assert all(field in source for field in ["id", "name"]), "One or more expected fields are missing in variant.primary_source.source."

@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_allele_type(variant_id, genome_id):
    """DOCSTRING"""
    
    additional_fields = """
        allele_type{
            accession_id
            value
            url
            source{
                id
                name
                description
                url
                release
            }
            
        }
    """
    
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id,variant_id=variant_id)
    query_data = {"query": query}
    (success, result) = await graphql(executable_schema, query_data, context_value=context(request={}))
    assert success, f"Query execution failed for variant {variant_id}. Result: {result}"
    
    variant = result.get("data", {}).get("variant")
    assert "allele_type" in variant, "Missing 'allele_type' field in query result."
    
    # Make sure fields are present in allele_type subfields
    # assignment_method
    allele_type = variant.get("allele_type")
    assert all(field in allele_type for field in ["accession_id", "value", "url", "source"]), "One or more expected fields are missing in variant.allele_type."
    # source
    source = allele_type.get("source")
    assert all(field in source for field in ["id", "name", "description", "url", "release"]), "One or more expected fields are missing in variant.allele_type.source."
    
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", VARIANT_TEST_CASES)
async def test_basic_slice(variant_id, genome_id):
    """Test that the very simplest query is working - are name and type
    being returned in response to a variant-based query."""
    
    additional_fields = """
        slice{
            region{
                name
            }
            location{
                start
            }
            strand{
                code
            }
            default
        }
        
    """
    
    template = query_template(additional_fields)
    query = template.substitute(genome_id=genome_id,variant_id=variant_id)
    query_data = {"query": query}
    (success, result) = await graphql(executable_schema, query_data, context_value=context(request={}))
    assert success, f"Query execution failed for variant {variant_id}. Result: {result}"
    
    # get fields and check that they exist
    variant = result.get("data", {}).get("variant")

    # variant should not be none
    assert variant is not None, f"Variant is None. Errors: {result.get('errors', 'No error info')}"
    
    # given that variant isn't none, check that name and type exist
    slice_field = variant.get("slice")
    assert all(field in slice_field for field in ["region", "location", "strand", "default"]), "One or more expected fields are missing in variant.allele_type.source."
    