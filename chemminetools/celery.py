from __future__ import absolute_import, unicode_literals
import os
from django.conf import settings
from celery import Celery

import tools

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'chemminetools.settings')

brokerHost="172.31.37.64" # production IP
if settings.DEBUG:
    brokerHost = "172.31.35.118" #development IP



#app = Celery('chemminetools', backend='cache+memcached://127.0.0.1:11211/', broker='amqp://guest@localhost//')
#app = Celery('chemminetools', backend='file:///srv/shared_jobs/celery/results', broker='amqp://guest@localhost//')
#app = Celery('chemminetools', backend='file:///srv/shared_jobs/celery/results', broker='amqp://remote_worker:askfj4l3nbb43@'+brokerHost+'//')
app = Celery('chemminetools')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()


#@app.task(bind=True)
#def debug_task(self):
#    print('Request: {0!r}'.format(self.request))


