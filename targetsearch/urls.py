#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from .views import *

urlpatterns = [
        url(r'^$', search, name='search'),
        url(r'^new', newTS, name='newTS'),
        url(r'^bs4test', bs4test, name='bs4test'),
        ]
