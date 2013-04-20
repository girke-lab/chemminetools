#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from views import *
urlpatterns = patterns('', url(r'^$', 'similarityworkbench.views.ui',
                       name='similarity-workbench-ui'),
                       url(r'^renderer/(?P<smiles>\S*)$', renderer))
