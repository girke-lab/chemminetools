from django.conf.urls.defaults import *
from views import uploadCompound

urlpatterns = patterns('',
	(r'^jscripts/', 'django.views.static.serve', { 'document_root': 'jscripts' }),
	url(r'',  uploadCompound),
)
