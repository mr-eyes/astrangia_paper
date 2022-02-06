from posixpath import split
import subprocess
import sys
import os
from typing import List
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi
from ncbi.datasets import GenomeApi
from ncbi.datasets.openapi.api import taxonomy_api
from ncbi.datasets.openapi.api import genome_api
import json


def get_summary_json(tax_id):
    with DatasetsApiClient() as api_client:
        api = GenomeApi(api_client)
        response = api.assembly_descriptors_by_taxon(
            taxon=tax_id, 
            async_req=True,
            filters_reference_only = True,
            )
        result = response.get().to_dict()
        if "total_count" in result:
            return int(result["total_count"])
        else:
            return 0
        

def search_up(lineages_ids, tax_id):
    for tx in lineages_ids[::-1]:
        if tx == tax_id:
            continue
        else:
            found_matches = get_summary_json(tx)
            if found_matches:
                return(tx)          

def name_to_taxid(species_name):
    ps_species = subprocess.Popen(("echo", species_name), stdout=subprocess.PIPE)
    ps_name2tax = subprocess.Popen(
            ("taxonkit", "name2taxid"), stdin=ps_species.stdout, stdout=subprocess.PIPE
        )
    ps_species.stdout.close()
    output = str(ps_name2tax.communicate()[0], encoding="utf-8").strip()
    return(output.split()[-1])

def tax_to_lineage(tax_id):
    """get lineage list from species name or tax ID

    Args:
        tax_id (str): tax_id or species name

    Returns:
        tuple([lineage_names], [lineage_ids]): Tuple of two lists containing lineage names and IDs
    """
    ps_name2tax = subprocess.Popen(("echo", tax_id), stdout=subprocess.PIPE)
    ps_lineage = subprocess.Popen(
        ("taxonkit", "lineage", "-i", "100", "-t"),
        stdin=ps_name2tax.stdout,
        stdout=subprocess.PIPE,
    )
    ps_name2tax.stdout.close()
    output = str(ps_lineage.communicate()[0], encoding="utf-8").strip()
    splitted = output.split("\t")
    lineage_names = splitted[1].split(";")
    lineage_ids = splitted[2].split(";")
    return (
        lineage_names,
        lineage_ids,
    )


species_name = "Scrippsiella"
tax_id = name_to_taxid(species_name)
lineage_names, lineage_ids, = tax_to_lineage(tax_id)
nearest_tax = search_up(lineage_ids, tax_id)
print(f"original_tax: {tax_id}, nearest_available_tax: {nearest_tax}")
