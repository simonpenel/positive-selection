# positive-selection
Serveur web dédié à l'affichage de la sélection positive dans un arbre phylogénétique et un alignement.

# genere l'arbre a afficher
python3 genere_xml.py  data/trees.demo.CLU_000022_0_9  data/HOGENOMID_SPECIES.demo  data/syntLeft.dna.demo  data/syntRight.dna.demo
cp data/trees.demo.CLU_000022_0_9.xml input_tree.xml  

# importe les packges nodes necessaire au serveur node
nmp install

# lance le serveur
npm start
aller voir http://0.0.0.0:8080/
