#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from views import compound_detail

urlpatterns = patterns('', url(r'^cid_lookup/?$',
                       'compounddb.views.cid_lookup'),
                       url(r'^(?P<id>\d+)/png/?(?P<filename>\S*)$',
                       'compounddb.views.render_image'),
                       url(r'^(?P<id>\d+)/(?P<resource>\w*)/?(?P<filename>\S*)$'
                       , 'compounddb.views.compound_detail',
                       name='compound_detail'))
