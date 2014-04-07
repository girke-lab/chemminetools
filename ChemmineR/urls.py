#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.conf import settings
import logging

urlpatterns = patterns('', url(r'^runapp(.*)$', 'ChemmineR.views.runapp'
                       , name='runapp'),
    url(r'^listCMTools(.*)$', 'ChemmineR.views.listCMTools'
                       , name='listCMTools'),
    url(r'^showJob/(?P<task_id>.*)/$', 'ChemmineR.views.showJob', name='showJob'),
    url(r'^jobStatus(.*)$', 'ChemmineR.views.jobStatus'
                       , name='jobStatus'),
    url(r'^getConverter(.*)$', 'ChemmineR.views.getConverter'
                       , name='getConverter'),
    url(r'^jobResult(.*)$', 'ChemmineR.views.jobResult'
                       , name='jobResult'),
    url(r'^toolDetails(.*)$', 'ChemmineR.views.toolDetails'
                       , name='toolDetails'),
    url(r'^launchCMTool(.*)$', 'ChemmineR.views.launchCMTool'
                       , name='launchCMTool'))
