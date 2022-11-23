description:
------

Bioinformatics project in Prof. Dan Mishmar and Prof. Chen Ceasar labs in Ben Gurion university.
The project deals with co-evolution between the nuclear genome to the mitochondria genome.
In the project we find the 500 closest homology proteins to each of the mitochondria proteins.
We built phylogenetic tree for each protein and check the changes through the evolution.
in addition, we compare all the genes with the 16-s and 12-s RNA from the nuclear genome using muscle.

project goal:
----
- research the co-evolution between the gene of the two subunits of the ribosome.

project structures:
------------
- blast-ribosome-protein:
     - finding homologies using BLASTp to all the mitochondria proteins using bio-python.
     - convert all the data into csv using pandas.

- check_results:
     - go over all tha csv files and check that in each file all the proteins have the same name.
     - finding proteins with wrong name and delete them using regex.
     - making a list of all the organisms that share the results of all the proteins.
     - create the data to the 12-s and 16-s RNA.