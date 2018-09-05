#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *
urlpatterns = [ url(r'^$', ui, name='similarity-workbench-ui'),
               url(r'^renderer/(?P<smiles>\S*)$', renderer)]
