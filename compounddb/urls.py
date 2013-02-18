from django.conf.urls import patterns, include, url

urlpatterns = patterns('',
	# url(r'^compounds/(?P<cid>[a-zA-Z0-9 _-]+)/png',
	# 	compounddb.views.compound_image', name='compound_image'),
	url(r'^compounds/(?P<cid>[a-zA-Z0-9 _-]+)/$',
		'compounddb.views.compound_detail', name='compound_detail'),
	# url(r'^_ajax$',
	#	'compounddb.views.ajax', name='ajax'),
)
