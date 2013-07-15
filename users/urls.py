#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from views import * 

urlpatterns = patterns('',
                       url(r'^register/$',
                       'users.views.register')
                       )
