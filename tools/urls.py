#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *

urlpatterns = [
    url(r'^launch_job/(?P<category>\w*)/$', launch_job),
    url(r'^list_jobs/$', list_jobs),
    url(r'^view_job/(?P<job_id>\d+)/(?P<resource>\w*)/(?P<filename>\S*)$'
        , view_job),
    url(r'^view_job/(?P<job_id>\d+)/(?P<resource>\w*)$', view_job),
    url(r'^', list_jobs),
    ]
