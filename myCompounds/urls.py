#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from views import *

urlpatterns = patterns('', url(r'^addCompounds/?(?P<resource>\w*)/?(?P<job_id>\d*)/?$', uploadCompound),
                       url(r'^(?P<resource>\S*)$', showCompounds))
