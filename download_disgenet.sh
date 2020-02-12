#!/bin/bash
wget -c https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_to_disease_ALL.tsv.gz -O data/disease_to_disease_ALL.tsv.gz
wget -c https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz -O data/all_gene_disease_associations.tsv.gz
wget -c https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_variant_disease_associations.tsv.gz -O data/all_variant_disease_associations.tsv.gz
wget -c https://www.disgenet.org/static/disgenet_ap1/files/downloads/variant_to_gene_mappings.tsv.gz -O data/variant_to_gene_mappings.tsv.gz
