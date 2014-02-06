#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.conf import settings
import logging

urlpatterns = patterns('', url(r'^runapp(.*)$', 'ChemmineR.views.runapp'
                       , name='runapp'),
    url(r'^listCMTools(.*)$', 'ChemmineR.views.listCMTools'
                       , name='listCMTools'),
    url(r'^jobStatus(.*)$', 'ChemmineR.views.jobStatus'
                       , name='jobStatus'),
    url(r'^jobResult(.*)$', 'ChemmineR.views.jobResult'
                       , name='jobResult'),
    url(r'^launchCMTool(.*)$', 'ChemmineR.views.launchCMTool'
                       , name='launchCMTool'))
