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

from typing import Dict, Optional, List, Any
import json
from ariadne import QueryType, ObjectType
from graphql import GraphQLResolveInfo
import subprocess

from graphql_service.resolver.exceptions import (
    VariantNotFoundError
)

# Define Query types for GraphQL
# Don't forget to import these into ariadne_app.py if you add a new type

QUERY_TYPE = QueryType()
VARIANT_TYPE = ObjectType("Variant")
VARIANT_ALLELE_TYPE = ObjectType("VariantAllele")

@QUERY_TYPE.field("variant")
async def resolve_variant(
        _,
        info: GraphQLResolveInfo,
        by_id: Dict[str, str] = None,
) -> Dict:
    "Load variants via variant id"
    
    query = {
        "type": "Variant",
        "variant_id": by_id["variant_id"],
        "genome_id": by_id["genome_id"],
    }
    file_client = info.context["file_client"]
    result = file_client.get_variant_record(by_id["variant_id"])
    if not result:
        raise VariantNotFoundError(by_id["variant_id"])
    return result

@VARIANT_TYPE.field("primary_source")
def primary_source(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load primary source for variant
    """
    return variant.get_primary_source()

@VARIANT_TYPE.field("allele_type")
def allele_type(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load allele_type for variant
    """
    return variant.get_allele_type(variant.alts)

@VARIANT_TYPE.field("slice")
def slice(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load slice for variant
    """
    return variant.get_slice(variant.alts)

@VARIANT_TYPE.field("alleles")
def alleles(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load alleles for variant
    """
    return variant.get_alleles()

@QUERY_TYPE.field("version")
def resolve_api(
    _: None, info: GraphQLResolveInfo
):  # the second argument must be named `info` to avoid a NameError
    return {"api": {"major": "0", "minor": "1", "patch": "0-beta"}}

