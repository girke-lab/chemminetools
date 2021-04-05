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
        path('tsgraph/<id_type>/<table_name>/<ids>/', tsGraph, name='tsGraph'),

        path('<initial_ids>',newTS), # this needs to be last
        ]
