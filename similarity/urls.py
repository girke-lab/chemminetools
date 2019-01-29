#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *

urlpatterns = [
        url(r'^$', showCompoundGrid),
        url(r'^mcs/$',mcs),
        url(r'^ap/$',ap),
        url(r'^renderer/(?P<smiles>\S*)$', renderer)
        ]
