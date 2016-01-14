###Imports
```python
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio import Entrez
from toolz import compose
```

### Getting databse handle
```python
server = BioSeqDatabase.open_database(driver="sqlite3",db="biodatabase")
dengue = server.open_database("Dengue")
db = server.new_database("Influenza", description="Just for testing")
In[]: db.adaptor.list_biodatabase_names()
out[] ['Dengue', 'Influenza']
```
###Queries
queries can be done manually, from a number of different tables.
```python
sql = compose(list, db.adaptor.cursor.execute)
sql('SELECT * FROM bioentry') # basic info accession, desc.
sql(r"select term_id from term where name='source'")
```
or via the lookup method:
```python
db.lookup(accession='EU569704')
```

Valid lookup kwargs:
```
display_id
name
accession
primary_id
version
gi
```
also:
```python
def get_all_records(db):
    return map(db.__getitem__, db.get_all_primary_ids())
```
Get All countries/locations:
```python
In[206]: sql(r"""SELECT * FROM seqfeature_qualifier_value WHERE term_id =
    (select term_id from term where name='country')""")

Out[206]: [(1, 15, 1, u'USA: Puerto Rico')]
```
from the SeqRecord object, usually the first feature, and usually of type source:
```python
In[]: source = filter(lambda x: x.type == 'source', rec.features)
In[]: source[0].qualifiers['country']
Out[]: 'USA: Puerto Rico'
```

###RegEx search
sqlite does not have regex search by default so have to register the function

```python
import re
def re_fn(expr, item):
    reg = re.compile(expr, re.I)
    return reg.search(item) is not None
db.adaptor.conn.create_function("REGEXP", 2, re_fn)
```
```python
# works now
sql(r"select * from seqfeature_qualifier_value where  value REGEXP '.*' ")

sql(r"select value from seqfeature_qualifier_value where value REGEXP '.*USA.*'")
# USA: Puerto Rico
```



###Inserting Records
```python
Entrez.email = "whatever@gmail.com"
handle = Entrez.efetch(db="nuccore", id="EU569704", rettype="gb", retmode="text") #id=accession
count = db.load(SeqIO.parse(handle, "genbank"))
server.commit()
```


### Other Stuff
map schema for where those features are and the like
Tables: seqfeature, location, comment, etc.

