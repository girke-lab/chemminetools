#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *

urlpatterns = [
        url(r'^$', home, name='home'),
        url(r'^getTargets/$', getTargets, name='getTargets'),
        url(r'^getTargetsTSV/$', getTargetsTSV, name='getTargetsTSV'),
        url(r'^getChembl/$', getChembl, name='getChembl'),
        url(r'^getChemblTSV/$', getChemblTSV, name='getChemblTSV'),
        ]
