#!/bin/bash
cat fastaheaders*-nrdb.gff > /mnt/gluster/amlinz/unprocessed-nrdb.gff
#sed -e '9311290d' /mnt/gluster/amlinz/unprocessed-nrdb.gff > /mnt/gluster/amlinz/semiprocessed-nrdb.gff
sed -i -e 's/ID=>/ID=/' /mnt/gluster/amlinz/unprocessed-nrdb.gff
gzip /mnt/gluster/amlinz/unprocessed-nrdb.gff
