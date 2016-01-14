from Bio import GenBank
from Bio import Entrez
from BioSQL import BioSeqDatabase

# Should read these from settings at some point
dbpath = 'biosql.sqlite3'
dbname = 'local_db'
Entrez.email = 'your.email@example.com'

server = BioSeqDatabase.open_database(driver='sqlite3', db=dbpath)
db = server[dbname]

parser = GenBank.FeatureParser()
loadgb = lambda _id: db.load(GenBank.Iterator(Entrez.efetch(db='nucleotide', id=_id, rettype='gb', retmode='text'), parser))

for id in open('accession.lst').read().splitlines():
    print "Loading %s" % id
    loadgb(id)
server.adaptor.commit()
