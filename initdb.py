#!/usr/bin/env python

import os
import os.path
import sys

from BioSQL import BioSeqDatabase
import sh

# Should read these from settings at some point
dbpath = 'biosql.sqlite3'
dbname = 'local_db'

if os.path.exists(dbpath):
    print "%s already exists" % dbpath
    print "You can remove it and rerun this if you want to reinitialize"
    sys.exit(1)

print "Downloading BioSQL Schema"
print sh.wget(
    'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-sqlite.sql',
    '-O',
    'biosqldb-sqlite.sql'
)
print "Initializing database with schema"
print sh.sqlite3('-init', 'biosqldb-sqlite.sql', dbpath, _in='.exit')

print "Creating BioSQL database"
server = BioSeqDatabase.open_database(driver='sqlite3', db='biosql.sqlite3')
db = server.new_database(dbname)
server.adaptor.commit()
