from django.conf.urls.defaults import *
from views import * 

urlpatterns = patterns('',
	url(r'^manage_application/?$', manage_application),
	url(r'^launch_job/?$', launch_job),
	url(r'^view_job/(?P<job_id>\d+)$', view_job),
)
