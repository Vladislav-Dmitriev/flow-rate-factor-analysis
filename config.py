from os import getenv

KPROD_API_HOST = getenv('KPROD_API_HOST', default='kprod_services_api')
KPROD_API_PORT = getenv('KPROD_API_PORT', default='8080')

CLIENT_NAME = getenv('CLIENT_NAME', default='kprod')

LOG_LEVEL = getenv('LOG_LEVEL', default='DEBUG')

