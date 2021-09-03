# Sélection positive
Serveur web dédié à l'affichage de la sélection positive dans un arbre phylogénétique et un alignement.

# Cloner le dépôt


# Générer l'arbre à afficher:
python3 genere_xml.py  data/trees.demo.CLU_000022_0_9  data/HOGENOMID_SPECIES.demo  data/syntLeft.dna.demo  data/syntRight.dna.demo
cp data/trees.demo.CLU_000022_0_9.xml input_tree.xml  

# Importer les packages nodes nécéssaires au serveur node
npm install

# Lancer le serveur
npm start
aller voir http://0.0.0.0:8080/
