from os import access
import sys
import json
import subprocess
from tqdm import tqdm
import sys
import ncbi.datasets.openapi
from ncbi.datasets.openapi.api import taxonomy_api
from ncbi.datasets.openapi.model.v1_taxonomy_metadata_request import (
    V1TaxonomyMetadataRequest,
)
from ncbi.datasets.openapi.model.rpc_status import RpcStatus
from ncbi.datasets.openapi.model.v1_taxonomy_metadata_request_content_type import (
    V1TaxonomyMetadataRequestContentType,
)
from ncbi.datasets.openapi.model.v1_taxonomy_metadata_response import (
    V1TaxonomyMetadataResponse,
)
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets import GenomeApi


configuration = ncbi.datasets.openapi.Configuration(
    host="https://api.ncbi.nlm.nih.gov/datasets/v1"
)


def nearest_available_taxon(lineage: list):
    with DatasetsApiClient() as api_client:
        api = GenomeApi(api_client)
        taxon = str(lineage.pop())
        response = api.assembly_descriptors_by_taxon(
            taxon=taxon,
            async_req=True,
            filters_reference_only=True,
        )
        result = response.get().to_dict()
        if "total_count" in result:
            accessions = [
                assembly["assembly"]["assembly_accession"]
                for assembly in result["assemblies"]
            ]
            return (taxon, accessions)  # , result)
        else:
            return nearest_available_taxon(lineage)


def tax_to_lineage(tax_id):
    """get lineage list from species name or tax ID

    Args:
        tax_id (str): tax_id or species name

    Returns:
        tuple([lineage_names], [lineage_ids]): Tuple of two lists containing lineage names and IDs
    """
    tax_id = str(tax_id)
    ps_lineage = subprocess.Popen(
        [
            "sh",
            "-c",
            f" echo {tax_id} |"
            + ' taxonkit reformat --taxid-field 1 --format "{k},{p},{c},{o},{f},{g},{s},{t}"',
        ],
        stdout=subprocess.PIPE,
    ).communicate()
    lineage = str(ps_lineage[0], encoding="utf8").split("\t")[-1]
    return lineage


def fetch_tax_info(tax_or_name):
    with ncbi.datasets.openapi.ApiClient(configuration) as api_client:
        api_instance = taxonomy_api.TaxonomyApi(api_client)
        v1_taxonomy_metadata_request = V1TaxonomyMetadataRequest(
            taxons=[str(tax_or_name)],
            returned_content=V1TaxonomyMetadataRequestContentType("COMPLETE"),
        )
        try:
            # Use taxonomic identifiers to get taxonomic metadata by post
            api_response = api_instance.taxonomy_metadata_post(
                v1_taxonomy_metadata_request
            )
            return api_response.to_dict()
        except ncbi.datasets.openapi.ApiException as e:
            print(f"Exception when calling TaxonomyApi->taxonomy_metadata_post: {e}")


blast_results = "blast_results.tsv"

# species scan
species = set()
with open(
    blast_results,
) as blast:
    next(blast)
    for line in blast:
        _species = []
        for _ in line.strip().split()[8].split("_"):
            if len(_) > 3:
                _species.append(_)
        species.add(" ".join(_species))

print(f"Total number of organisms: {len(species)}")

taxon_to_lineage = dict()

with open("nearest_taxon.tsv", "w") as NEAREST_TSV, open(
    "sourmash_lineages.csv", "w"
) as SMASH_LIN:
    NEAREST_TSV.write(
        "organism\ttaxon\tnearest_organism\tnearest_taxon\toriginal_found?\taccessions\n"
    )
    SMASH_LIN.write(
        "accession,taxid,superkingdom,phylum,class,order,family,genus,species,strain\n"
    )
    for sp in tqdm(species):
        # species to tax
        tax_info = fetch_tax_info(sp)["taxonomy_nodes"][0]["taxonomy"]
        # tax to lineag
        organism_name = tax_info["organism_name"]
        organism_taxon = str(tax_info["tax_id"])
        full_lineage = tax_info["lineage"] + [organism_taxon]
        nearest_taxon, accessions = nearest_available_taxon(full_lineage)
        nearest_organism_name = fetch_tax_info(nearest_taxon)["taxonomy_nodes"][0][
            "taxonomy"
        ]["organism_name"]
        NEAREST_TSV.write(
            f"{organism_name}\t{organism_taxon}\t{nearest_organism_name}\t{nearest_taxon}\t{organism_taxon == nearest_taxon}\t{','.join(accessions)}\n"
        )
        if nearest_taxon not in taxon_to_lineage:
            taxon_to_lineage[nearest_taxon] = tax_to_lineage(nearest_taxon)

        for accession in accessions:
            accession = accession.split(".")[0]
            lineage = taxon_to_lineage[nearest_taxon]
            sourmash_lin = f"{accession}, {nearest_taxon}, {lineage}"
            SMASH_LIN.write(sourmash_lin)
