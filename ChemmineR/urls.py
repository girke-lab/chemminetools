#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.conf import settings
import logging

urlpatterns = patterns('', url(r'^runapp(.*)$', 'ChemmineR.views.runapp'
                       , name='runapp'))

