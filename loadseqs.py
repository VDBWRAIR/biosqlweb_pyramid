#!/usr/bin/env python

from Bio import GenBank
from Bio import Entrez
from BioSQL import BioSeqDatabase
import sys

# Should read these from settings at some point
dbpath = 'biosql.sqlite3'
dbname = 'local_db'
Entrez.email = 'your.email@example.com'

server = BioSeqDatabase.open_database(driver='sqlite3', db=dbpath)
db = server[dbname]

parser = GenBank.FeatureParser()
loadgb = lambda _id: db.load(GenBank.Iterator(Entrez.efetch(db='nucleotide', id=_id, rettype='gb', retmode='text'), parser))

ACCESSIONS_FILE = 'accessions.lst' if len(sys.argv) < 2 else sys.argv[1]
for id in open(ACCESSIONS_FILE):
    print "Loading %s" % id
    loadgb(id)
server.adaptor.commit()
