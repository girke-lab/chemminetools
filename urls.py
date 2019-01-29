#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
#from django.conf.urls.defaults import *
from django.conf.urls import *
from django.contrib import admin
from django.conf import settings
#from django.views.generic.simple import direct_to_template
from django.views.generic import TemplateView
from django.views.generic import RedirectView
from django.contrib.auth.views import login, logout


admin.autodiscover()
print('\n before urls in general')
urlpatterns =[
    url(r'^admin/', include(admin.site.urls)),
    url(r'^eisearch/', include('eisearch.urls')),
    url(r'^accounts/', include('userena.urls')),
    url(r'^compounds/', include('compounddb.urls')),
    url(r'^my[Cc]ompounds/', include('myCompounds.urls')),
    url(r'^tools/', include('tools.urls')),
    url(r'^drugbank/', include('drugbank.urls')),
    url(r'^search/?',  RedirectView.as_view(url='/eisearch/query/')), 
    #url(r'^similarity/', include('similarityworkbench.urls')),
    url(r'^similarity/', include('similarity.urls')),
    url(r'^ChemmineR/', include('ChemmineR.urls')),
    url(r'^robots\.txt/?$', TemplateView.as_view(template_name='robots.txt')),
    url(r'^ei/?',  RedirectView.as_view(url='/downloads/')), 
    url(r'^', include('cms.urls')),
    ]
print('\n after urls in general')
if settings.DEBUG:
    import django
    urlpatterns = [url(r'^working/(?P<path>.*)$',
                           django.views.static.serve,
                           {'document_root': settings.MEDIA_ROOT,
                           'show_indexes': False}), 
                   url(r'', include('django.contrib.staticfiles.urls'))
                    ] + urlpatterns
