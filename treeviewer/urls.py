from django.conf.urls.defaults import *
urlpatterns = patterns('',
    url(r'^([a-fA-F0-9]+)/$', 'treeviewer.views.viewer', name='tree-viewer'),
    url(r'^([a-fA-F0-9]+)/data/$', 'treeviewer.views.userdata',
        name='userdata'),
    url(r'^([a-fA-F0-9]+)/tile/$', 'treeviewer.views.serv_tile',
        name='tile'),
    url(r'^([a-fA-F0-9]+)/full/$', 'treeviewer.views.serv_fullsize',
        name='fullsize'),
    )
