from django.conf.urls.defaults import *
from django.contrib import admin
from django.conf import settings

# enable cron
import django_cron
django_cron.autodiscover()

admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    (r'^compounds/', include('compounddb.urls')),
    (r'^my[Cc]ompounds/', include('myCompounds.urls')),
    (r'^tools/', include('tools.urls')),
    (r'^search/pug/', include('pugsearch.urls')),
    (r'^search/', include('eis.urls')),
    (r'^similarity/', include('similarityworkbench.urls')),
    url(r'^', include('cms.urls')),
)

if settings.DEBUG:
    urlpatterns = patterns('',
        url(r'^working/(?P<path>.*)$', 'django.views.static.serve',
        {'document_root': settings.MEDIA_ROOT, 'show_indexes': False}),
        url(r'', include('django.contrib.staticfiles.urls')),
    ) + urlpatterns
