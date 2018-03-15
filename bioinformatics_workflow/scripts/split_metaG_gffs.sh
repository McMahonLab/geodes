#!/bin/bash

mkdir metaG_gffs
cd metaG_gffs
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES005.assembled.gff GEODES005.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES006.assembled.gff GEODES006.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES057.assembled.gff GEODES057.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES058.assembled.gff GEODES058.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES117.assembled.gff GEODES117.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES118.assembled.gff GEODES118.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES166.assembled.gff GEODES165.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES166.assembled.gff GEODES166.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES167.assembled.gff GEODES167.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES168.assembled.gff GEODES168.

# make list of files to run
ls GEODES005* > ../005metaG_gffs.txt
ls GEODES006* > ../006metaG_gffs.txt
ls GEODES057* > ../057metaG_gffs.txt
ls GEODES058* > ../058metaG_gffs.txt
ls GEODES117* > ../117metaG_gffs.txt
ls GEODES118* > ../118metaG_gffs.txt
ls GEODES165* > ../165metaG_gffs.txt
ls GEODES166* > ../166metaG_gffs.txt
ls GEODES167* > ../167metaG_gffs.txt
ls GEODES168* > ../168metaG_gffs.txt
cd ..

