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
    url(r'^', include('cms.urls')),
)

if settings.DEBUG:
    urlpatterns = patterns('',
        url(r'^media/(?P<path>.*)$', 'django.views.static.serve',
        {'document_root': settings.MEDIA_ROOT, 'show_indexes': True}),
        url(r'', include('django.contrib.staticfiles.urls')),
    ) + urlpatterns
