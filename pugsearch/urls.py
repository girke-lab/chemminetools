#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.conf import settings
import logging

urlpatterns = patterns('', url(r'^(?P<sid>[0-9a-z]+)/$',
                       'pugsearch.views.search', name='pug_search'))
