#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.conf import settings
import logging
from views import *

urlpatterns = patterns('', 
    url(r'^query/$', search),
    url(r'^getStructures/(?P<job_id>\d+)/(?P<format>workbench|smiles|sdf)', getStructures)
    )
