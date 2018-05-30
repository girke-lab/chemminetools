
from django.conf.urls.defaults import *
from drugbank.views import *
from django.views.decorators.csrf import csrf_exempt

print '\n before url'
app_name = 'drugbank'
urlpatterns = patterns('', url(r'^$',
                       csrf_exempt(drugbank_lookup),name = 'drugbank-lookup'))
print '\n after url'


