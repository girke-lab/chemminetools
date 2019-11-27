#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django.conf.urls import *
from django.urls import path
from .views import *

urlpatterns = [
        url(r'^$', newTS, name='newTS'),
        path('compoundNames/<query>',compoundNames),
        ]
