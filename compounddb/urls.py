from django.conf.urls.defaults import *
from views import compound_detail

urlpatterns = patterns('',
	url(r'^(?P<id>\d+)/png/?(?P<filename>\S*)$',
		'compounddb.views.render_image'),
	url(r'^(?P<id>\d+)/(?P<resource>\w*)/?(?P<filename>\S*)$',
		'compounddb.views.compound_detail', name='compound_detail'),
	# url(r'^_ajax$',
	#	'compounddb.views.ajax', name='ajax'),
)
