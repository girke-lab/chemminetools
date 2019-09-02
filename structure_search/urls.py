#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from django.conf import settings
import logging
from .views import *

urlpatterns =[
    url(r'^query/$', search),
    url(r'^getStructures/(?P<job_id>\d+)/(?P<format>workbench|smiles|sdf)', getStructures)
    ]
