#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
#from .views import compound_detail
from . import views

urlpatterns = [url(r'^cid_lookup/?$', views.cid_lookup),
               url(r'^(?P<id>\d+)/png/?(?P<filename>\S*)$', views.render_image),
               url(r'^(?P<id>\d+)/(?P<resource>\w*)/?(?P<filename>\S*)$' , views.compound_detail, name='compound_detail')]
