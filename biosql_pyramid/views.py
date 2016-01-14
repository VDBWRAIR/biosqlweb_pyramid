import collections

from pyramid.view import view_config
from mako.template import Template

from BioSQL import BioSeqDatabase

dbpath = 'biosql.sqlite3'
dbname = 'local_db'

def _get_db(dbpath=dbpath, db=dbname):
    server = BioSeqDatabase.open_database(
        driver='sqlite3', db=dbpath
    )
    return server[db]


#@view_config(route_name='home', renderer='templates/mytemplate.pt')
#def my_view(request):
    #return {'project': 'biosql_pyramid'}

@view_config(route_name='summary', renderer='templates/summary.mako')
class SummaryView(object):
    def __init__(self, request):
        self.request = request

    def __call__(self):
        return {'records': self._get_records('')}

    def _get_records(self, biodb_name):
        """Retrieve BioSQL records in the database.
        """
        bioentries = self._get_bioentries(biodb_name)
        records = []
        for bioentry in bioentries:
            key = bioentry[0]
            entry = bioentry[1]
            b_tmpl = Template(bioentry_template)
            retrieve_url = "bioentry_details?bioentry_key=%s" % key
            records.append(
                b_tmpl.render(
                    id=entry.id,
                    description=entry.description, 
                    retrieve_url=retrieve_url
                )
            )
        return records

    def _get_bioentries(self, biodb_name, start=0, limit=10):
        """Retreive bioentries associated with the database.
        """
        start = int(self.request.params.get('start', start))
        limit = int(self.request.params.get('limit', limit))
        biodb = _get_db()
        return biodb.items()[start:start+limit]

@view_config(route_name='bioentry_details', renderer='templates/bioentry_details.mako')
def bioentry_details(request):
    """Retrieve full details for a bioentry based on the internal key.

    This could also use accession numbers or unique identifiers here.
    """
    bioentry_key = request.params.get('bioentry_key', '')
    biodb = _get_db()
    bioentry = biodb.get_Seq_by_primary_id(bioentry_key)
    annotations = _build_annotations(bioentry)
    annotations.sort()
    seqstr = str(bioentry.seq)
    seqstr = [
        seqstr[pos:pos+80] for pos in range(0, len(seqstr), 80)
    ]
    return {'annotations': annotations, 'sequence': seqstr}

def _build_annotations(bioentry):
    '''
    build a simple (key, value) list from dictionary but for all values that are
    iterable, " ,".join(value) them
    '''
    ann = {}
    for key, value in bioentry.annotations.items():
        if key.lower() == 'references':
            strval = [(ref.title, ref.authors) for ref in value]
        elif isinstance(value, list):
            strval = ", ".join(map(str, value))
        else:
            strval = str(value)
        ann[key] = strval
    return ann.items()

bioentry_template = """
<h3><a href="${retrieve_url}">${id} ${description}</a></h3>
<div>
</div>
"""
