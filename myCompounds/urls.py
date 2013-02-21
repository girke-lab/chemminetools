from django.conf.urls.defaults import *
from views import uploadCompound

urlpatterns = patterns('',
	url(r'',  uploadCompound),
)
