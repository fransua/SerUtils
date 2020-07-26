#! /usr/bin/env/python
"""
"""
from subprocess import Popen

species = 'hsapiens'
species = 'mmusculus'

xml = f'''
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "GFF" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "{species}_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "strand" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "cds_start" />
		<Attribute name = "cds_end" />
		<Attribute name = "exon_chrom_start" />
		<Attribute name = "exon_chrom_end" />
		<Attribute name = "external_gene_name" />
		<Attribute name = "gene_biotype" />
	</Dataset>
</Query>
'''

xml = ''.join(l.strip() for l in xml.split('\n'))

Popen(f"axel -n 10 -a --output {species}_result.txt 'http://www.ensembl.org/biomart/martservice?query={xml}'", shell=True).communicate()

