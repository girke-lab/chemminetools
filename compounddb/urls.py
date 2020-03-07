#!/usr/bin/python
# -*- coding: utf-8 -*-

# URL routes for /compounds/

from django.urls import path, re_path
from .views import *

urlpatterns = [path('cid_lookup/', cid_lookup),
               path('<int:id>/png/', render_image),
               path('cid/<cid>/png/', render_image),
               path('<int:id>/svg/', render_svg, name='render_svg'),
               path('chembl/<chembl_id>/svg/', render_chembl_svg, name='render_chembl_svg'),
               path('tagCompounds/', tagCompounds),
               path('tagCompounds/<action>/', tagCompounds),
               path('batch/<action>/', batchOperation),
               path('withTags/<tags>/count/', countCompoundsWithTags),
               path('cid/<cid>/', compound_detail, name='compound_detail'),
               path('<int:id>/', compound_detail, name='compound_detail'),
               path('cid/<cid>/<resource>/', compound_detail, name='compound_detail'),
               path('<int:id>/<resource>/', compound_detail, name='compound_detail'),
               ]
