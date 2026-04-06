file1='~/project/Au/regen/Xenopus/Xenopus_pep.fa'
type1='prot'
id1='Xe'
file2='~/project/Au/Au_genome/Au_pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60
Running blastx from 1 to 2 and tblastn from 2 to 1


file1='~/project/Au/regen/Schmidtea_mediterranea/GSE72389/GSE72389_smed_20140614.fa'
type1='nucl'
id1='Sc'
file2='~/project/Au/Au_genome/Au_pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60

file1='~/project/Au/regen/axolotl/Am_2.2_protein.fa'
type1='prot' #or 'prot' if file1 is a proteome
id1='Ax' #2-character ID (e.g. 'hu' for human)
file2='~/project/Au/Au_genome/combine.pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60

file1='~/project/Au/regen/Danio_rerio/Danio_rerio.GRCz10.pep.all.fa'
type1='prot' #or 'prot' if file1 is a proteome
id1='Da' #2-character ID (e.g. 'hu' for human)
file2='~/project/Au/Au_genome/combine.pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60


file1='~/project/Au/regen/lizard/Anolis_carolinensis.AnoCar2.0v2.pep.all.fa'
type1='prot' #or 'prot' if file1 is a proteome
id1='La' #2-character ID (e.g. 'hu' for human)
file2='~/project/Au/Au_genome/combine.pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60

file1='~/project/Au/regen/Xenopus/Xenopus_pep.fa'
type1='prot' #or 'prot' if file1 is a proteome
id1='XP' #2-character ID (e.g. 'hu' for human)
file2='~/project/Au/Au_genome/combine.pep.fa'
type2='prot'
id2='Au' #2-character ID (e.g. 'mo' for mouse)
bash ~/project/cross_species/SAMap/map_genes.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2 --threads 60

