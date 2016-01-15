import collections
import json
from StringIO import StringIO

from pyramid.view import view_config
from pyramid.response import Response, FileIter
import pyramid.httpexceptions as httpexceptions
from mako.template import Template

from BioSQL import BioSeqDatabase
from Bio import SeqIO

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
        #return {'records': self._get_records('')}
        return {'records': self._get_bioentries('')}

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
    info = _build_info(bioentry.annotations) + _build_info(bioentry.features[0].qualifiers)
    info.sort()
    seqstr = str(bioentry.seq)
    seqstr = [
        seqstr[pos:pos+80] for pos in range(0, len(seqstr), 80)
    ]
    return {'id': bioentry.id, 'info': info, 'sequence': seqstr}

@view_config(route_name='download_file')
def download_file(request):
    biodb = _get_db()
    bioentry_id = request.params.get('accession')
    outformat = request.params.get('format', 'fasta')
    if not bioentry_id:
        raise httpexceptions.HTTPBadRequest('No Accession given')
    try:
        bioentry = biodb.get_Seq_by_ver(bioentry_id)
    except Exception as e:
        raise httpexceptions.HTTPNotFound('No entry found for %s' % bioentry_id)
    try:
        out = bioentry.format(outformat)
    except ValueError as e:
        raise httpexceptions.HTTPBadRequest('%s is not a supported format' % outformat)
    output = StringIO(bioentry.format(outformat))
    response = request.response
    response.content_type = 'text/plain'
    response.content_disposition = 'attachment; filename="%s.%s"' % (bioentry_id, outformat)
    response.app_iter = FileIter(output)
    return response
    
def _build_info(info_dict):
    '''
    build a simple (key, value) list from dictionary but for all values that are
    iterable, " ,".join(value) them
    '''
    info = {}
    for key, value in info_dict.items():
        if key.lower() == 'references':
            strval = [(ref.title, ref.authors) for ref in value]
        elif isinstance(value, list):
            strval = ", ".join(map(str, value))
        else:
            strval = str(value)
        info[key] = strval
    return info.items()

@view_config(route_name='genbank_upload')
def genbank_upload(request):
    # XXX hack for os.linesep not being present; where did it go?
    # os.linesep = "\n"
    biodb = _get_db()
    handle = request.POST['upload_file'].file
    try:
        biodb.load(SeqIO.parse(handle, "genbank"))
        biodb.adaptor.commit()
    except Exception as e: # I would catch the IntegrityError, however, would have to somehow dynamically detect db driver and import that exception
        # columns identifier, biodatabase_id are not unique
        raise httpexceptions.HTTPClientError("Failed Saving record to database. Likely duplicate record")

    handle.close()
    return Response('Success', content_type='text/javascript')
