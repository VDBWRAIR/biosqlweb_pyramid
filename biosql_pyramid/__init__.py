from pyramid.config import Configurator


def main(global_config, **settings):
    """ This function returns a Pyramid WSGI application.
    """
    config = Configurator(settings=settings)
    config.include('pyramid_chameleon')
    config.add_static_view('static', 'static', cache_max_age=3600)
    #config.add_route('home', '/')
    config.add_route('summary', '/')
    config.add_route('bioentry_details', 'bioentry_details')
    config.add_route('genbank_upload', 'genbank_upload')
    config.add_route('download_file', 'download_file')
    config.scan()
    return config.make_wsgi_app()
