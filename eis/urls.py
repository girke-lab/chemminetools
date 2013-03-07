from django.conf.urls.defaults import *
from django.conf import settings
import logging


urlpatterns = patterns('',
    url(r'^$', 'eis.views.search', name="compound_search"),
	url(r'^wait/(?P<sid>[a-z0-9]+)/$', 'eis.views.wait', name="search_wait"),
	url(r'^result/(?P<sid>[a-z0-9]+)/$', 'eis.views.result', name="search_result"),
	url(r'^download/(?P<sid>[a-z0-9]+)/(?P<format>[a-z]*)/$', 'eis.views.download', name="download_search_result"),
)
