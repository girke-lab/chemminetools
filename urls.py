#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.conf.urls.defaults import *
from django.contrib import admin
from django.conf import settings
from django.views.generic.simple import direct_to_template
from django.views.generic import RedirectView
from django.contrib.auth.views import login, logout

# enable cron

import django_cron
django_cron.autodiscover()

admin.autodiscover()
print '\n before urls in general'
urlpatterns = patterns(
    r'',
    url(r'^admin/', include(admin.site.urls)),
    (r'^eisearch/', include('eisearch.urls')),
    (r'^accounts/', include('userena.urls')),
    (r'^compounds/', include(r'compounddb.urls')),
    (r'^my[Cc]ompounds/', include('myCompounds.urls')),
    (r'^tools/', include('tools.urls')),
    (r'^drugbank/', include('drugbank.urls')),
    (r'^search/?',  RedirectView.as_view(url='/eisearch/query/')), 
    (r'^similarity/', include('similarityworkbench.urls')),
    (r'^ChemmineR/', include('ChemmineR.urls')),
    (r'^robots\.txt/?$', direct_to_template, {'template': 'robots.txt',
     'mimetype': 'text/plain'}),
    (r'^ei/?',  RedirectView.as_view(url='/downloads/')), 
    url(r'^', include('cms.urls')),
    )
print '\n after urls in general'
if settings.DEBUG:
    urlpatterns = patterns(r'', url(r'^working/(?P<path>.*)$',
                           r'django.views.static.serve',
                           {'document_root': settings.MEDIA_ROOT,
                           'show_indexes': False}), url(r'',
                           include(r'django.contrib.staticfiles.urls'
                           ))) + urlpatterns
