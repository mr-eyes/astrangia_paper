import subprocess
from tqdm import tqdm
import ncbi.datasets.openapi
from ncbi.datasets.openapi.api import taxonomy_api
from ncbi.datasets.openapi.model.v1_taxonomy_metadata_request import (
    V1TaxonomyMetadataRequest,
)
from ncbi.datasets.openapi.model.v1_taxonomy_metadata_request_content_type import (
    V1TaxonomyMetadataRequestContentType,
)
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets import GenomeApi


configuration = ncbi.datasets.openapi.Configuration(
    host="https://api.ncbi.nlm.nih.gov/datasets/v1"
)


def nearest_available_taxon(lineage):
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
            accessions = []
            for assembly in result["assemblies"]:
                if assembly["assembly"]["assembly_category"] == "representative genome":
                    tax_id = assembly["assembly"]["org"]["tax_id"]
                    sci_name = (
                        assembly["assembly"]["org"]["sci_name"]
                        .lower()
                        .replace(" ", "_")
                    )
                    accession = assembly["assembly"]["assembly_accession"]
                    accessions.append((accession, sci_name, tax_id))
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
) as SMASH_LIN, open("accession_to_info.tsv", "w") as NAMES, open(
    "prepare_genomes.sh", "w"
) as SH:
    NEAREST_TSV.write(
        "organism\ttaxon\tnearest_organism\tnearest_taxon\toriginal_found?\taccessions\n"
    )
    SMASH_LIN.write(
        "accession,taxid,superkingdom,phylum,class,order,family,genus,species,strain\n"
    )
    NAMES.write("accession\tparent_taxon\ttaxon\tsci_name\n")
    SH.write("set -e\n")
    for sp in tqdm(species):
        # species to tax
        tax_info = fetch_tax_info(sp)["taxonomy_nodes"][0]["taxonomy"]
        # tax to lineag
        organism_name = tax_info["organism_name"]
        organism_taxon = str(tax_info["tax_id"])
        full_lineage = tax_info["lineage"] + [organism_taxon]
        nearest_taxon, accessions_info = nearest_available_taxon(full_lineage)
        accessions = [acc for acc, _, _ in accessions_info]
        nearest_organism_name = fetch_tax_info(nearest_taxon)["taxonomy_nodes"][0][
            "taxonomy"
        ]["organism_name"]
        NEAREST_TSV.write(
            f"{organism_name}\t{organism_taxon}\t{nearest_organism_name}\t{nearest_taxon}\t{organism_taxon == nearest_taxon}\t{','.join(accessions)}\n"
        )

        for acc in accessions_info:
            assembly_accession, sci_name, tax_id = acc
            NAMES.write(f"{assembly_accession}\t{nearest_taxon}\t{tax_id}\t{sci_name}\n")
            url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{assembly_accession}/download"
            SH.write(f"# {assembly_accession}\t{nearest_taxon}\t{tax_id}\t{sci_name}\n")
            SH.write(f"wget -N {url} -O {assembly_accession}.zip\n")
            SH.write(f"unzip {assembly_accession}.zip -d {assembly_accession}\n")
            SH.write(
                f"cat {assembly_accession}/ncbi_dataset/data/{assembly_accession}/*fna > {assembly_accession}.fna\n"
            )
            SH.write(f"rm -rf {assembly_accession}/\n\n")

        for acc in accessions_info:
            assembly_accession, sci_name, tax_id = acc
            accession = assembly_accession.split(".")[0]
            lineage = tax_to_lineage(tax_id)
            sourmash_lin = f"{accession},{tax_id},{lineage}"
            SMASH_LIN.write(sourmash_lin)
