#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.urls import path, re_path
from .views import *

urlpatterns = [
        path('', newTS, name='newTS'),
        path('drugIndTable', drugIndTable, name='drugIndTable'),
        path('compoundNames/<query>', compoundNames, name='compoundNames'),
        path('targetNames/<query>', targetNames, name='targetNames'),
        path('ajax/<action>', ajax, name='targetsearch-ajax'),
        path('detail/<id>', detailPage),
        #path('extanno/<id_type>/<cid>/', extAnno, name='extAnno'),
        path('extannobychembl/<chembl_id>/', extAnnoByChembl, name='extAnnoByChembl'),
        path('extannobychembl/<chembl_id>/<db>/', extAnnoByChembl, name='extAnnoByChembl'),
        path('filter/', tsFilter, name='tsFilter'),
        path('filter2/', tsFilter2, name='tsFilter2'),
        path('annofilter1/', tsAnnoFilter1, name='tsAnnoFilter1'),
        path('annofilter2/', tsAnnoFilter2, name='tsAnnoFilter2'),

        path('<initial_ids>',newTS), # this needs to be last
        ]
