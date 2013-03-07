from django.conf.urls.defaults import *
urlpatterns = patterns('',
	url(r'^$', 'similarityworkbench.views.ui', name="similarity-workbench-ui"),
)
