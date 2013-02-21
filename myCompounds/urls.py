from django.conf.urls.defaults import *
from views import * 

urlpatterns = patterns('',
	url(r'^addCompounds/?$', uploadCompound),
	url(r'',  showCompounds),
)
