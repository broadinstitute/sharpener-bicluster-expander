
from swagger_server.models.transformer_info import TransformerInfo
from swagger_server.models.parameter import Parameter
from swagger_server.models.transformer_query import TransformerQuery
from swagger_server.models.gene_info import GeneInfo
from swagger_server.models.gene_info_identifiers import GeneInfoIdentifiers
from swagger_server.models.attribute import Attribute

from ncats.translator.modules.gene.gene.gene_to_gene_bicluster_RNAseqDB import GeneToGeneBiclusters

import requests

def expander_info():
    """
        Return information for this expander
    """
    return TransformerInfo(
        name='Gene To Gene Bicluster',
        function='expander',
        description='Gene To Gene Bicluster',
        parameters=[],
        required_attributes = []
    )


def expand(query: TransformerQuery):
    """
        Execute this expander, find all genes correlated to query genes.
    """
    controls = {control.name: control.value for control in query.controls}

    gene_list = []
    genes = {}
    for gene in query.genes:
        gene_list.append(gene)
        genes[gene.gene_id] = gene

    gtgb = GeneToGeneBiclusters(query.genes)
    results = gtgb.results.to_dict('records') if len(gtgb.results)>0 else []

    for result in results:
        gene_id = result['hit_id']
        symbol = result['hit_symbol']
        hgnc = gene_id if gene_id.startswith('HGNC:') else None
        input_gene_id = result['input_id']
        input_symbol = result['input_symbol']
        interaction = result['score']

        if gene_id not in genes:
            gene = GeneInfo(
                gene_id=gene_id,
                identifiers=GeneInfoIdentifiers(hgnc=hgnc),
                attributes=[]
            )
            gene_list.append(gene)
            genes[gene_id] = gene

        gene = genes[gene_id]
        gene.attributes.append(
            Attribute(
                name='interaction score',
                value=str(interaction),
                source=NAME
            )
        )
        gene.attributes.append(
            Attribute(
                name='Interaction Input Gene Symbol',
                value=str(input_symbol),
                source=NAME
            )
        )
        gene.attributes.append(
            Attribute(
                name='Interaction Input Gene ID',
                value=str(input_gene_id),
                source=NAME
            )
        )

    return gene_list
