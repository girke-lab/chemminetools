from django.conf.urls.defaults import *
from views import * 

urlpatterns = patterns('',
	url(r'^manage_application/(?P<chooseForm>\w*)/$', manage_application),
	url(r'^launch_job/$', launch_job),
	url(r'^list_jobs/$', list_jobs),
	url(r'^view_job/(?P<job_id>\d+)/?(?P<resource>\w*)/(?P<filename>\S*)$', view_job),
	url(r'^', list_jobs),
)
