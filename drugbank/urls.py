from __future__ import print_function

from django.conf.urls import *
from drugbank.views import *
from django.views.decorators.csrf import csrf_exempt

app_name = 'drugbank'
urlpatterns = [ url(r'^$',
                       csrf_exempt(drugbank_lookup),name = 'drugbank-lookup')]


