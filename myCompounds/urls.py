#!/usr/bin/python
# -*- coding: utf-8 -*-

from .views import *
from django.urls import path, re_path

#app_name = 'myCompounds'
urlpatterns = [re_path(r'^addCompounds/?(?P<resource>\w*)/?(?P<job_id>\d*)/?$', uploadCompound, name='uploadCompound'),
               re_path(r'^download/(?P<outputFormat>(sdf|smi))$', downloadCompounds, name='downloadCompounds'),
               path('ajax/<action>', ajax, name='myCompounds-ajax'),
               re_path(r'^(?P<resource>\S*)$', showCompounds, name='showCompounds')]
