#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls import *
from django.conf import settings
import logging
from .views import *

urlpatterns = [url(r'^runapp(.*)$', runapp , name='runapp'),
    url(r'^listCMTools(.*)$', listCMTools , name='listCMTools'),
    url(r'^showJob/(?P<task_id>.*)/$', showJob, name='showJob'),
    url(r'^jobStatus(.*)$', jobStatus , name='jobStatus'),
    url(r'^getConverter(.*)$', getConverter , name='getConverter'),
    url(r'^jobResult(.*)$', jobResult , name='jobResult'),
    url(r'^toolDetails(.*)$', toolDetails , name='toolDetails'),
    url(r'^launchCMTool(.*)$', launchCMTool , name='launchCMTool')]
