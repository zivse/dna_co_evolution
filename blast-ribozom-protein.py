import ssl  # monkey patch for BioPython 1.68 & 1.69
import csv
import requests as r
from Bio import SeqIO
from io import StringIO
from bs4 import BeautifulSoup
from Bio.Blast import NCBIWWW
ssl._create_default_https_context = ssl._create_unverified_context


with open('protein-ribozom.txt') as f:
    proteins = f.readlines()
with open('names.csv', 'w', newline='') as csvfile:
    fieldnames = ['query_id', 'protein_name', 'query_len']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for protein in proteins:
        print(protein)
        baseUrl = "http://www.uniprot.org/uniprot/"
        currentUrl = baseUrl + protein.strip() + ".fasta"
        print(currentUrl)
        response = r.post(currentUrl)
        cData = ''.join(response.text)
        Seq = StringIO(cData)
        pSeq = list(SeqIO.parse(Seq, 'fasta'))[0].seq
        result_handle = NCBIWWW.qblast("blastp", "nt", pSeq)
        blast_results = result_handle.read()
        print(blast_results)
        Bs_data = BeautifulSoup(blast_results, "xml")
        query_id = Bs_data.find_all('BlastOutput_query-ID')[0].contents[0]
        protein_name = Bs_data.find_all('BlastOutput_query-def')[0].contents[0]
        query_len = Bs_data.find_all('BlastOutput_query-len')[0].contents[0]
        writer.writerow({'query_id': query_id, 'protein_name':protein_name, 'query_len': query_len})








