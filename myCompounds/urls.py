#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *

urlpatterns = [url(r'^addCompounds/?(?P<resource>\w*)/?(?P<job_id>\d*)/?$', uploadCompound),
                       url(r'^(?P<resource>\S*)$', showCompounds)]
