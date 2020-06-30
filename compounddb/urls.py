#!/usr/bin/python
# -*- coding: utf-8 -*-

# URL routes for /compounds/

from django.urls import path, re_path
from .views import *

urlpatterns = [path('cid_lookup/', cid_lookup),

               path('<int:id>/png', render_image),
               path('<int:id>/png/<filename>', render_image),
               path('cid/<cid>/png', render_image),
               path('cid/<cid>/png/<filename>', render_image),

               path('<int:id>/svg', render_svg, name='render_svg'),
               path('<int:id>/svg/<filename>', render_svg, name='render_svg'),

               path('chembl/<chembl_id>/svg', render_chembl_svg, name='render_chembl_svg'),
               path('chembl/<chembl_id>/svg/<filename>', render_chembl_svg, name='render_chembl_svg'),

               path('tagCompounds/', tagCompounds),
               path('tagCompounds/<action>/', tagCompounds),
               path('tagDuplicates/', tagDuplicateCompounds),
               path('batch/<action>/', batchOperation),
               path('withTags/<tags>/count/', countCompoundsWithTags),

               path('cid/<cid>/', compound_detail, name='compound_detail'),
               path('cid/<cid>/<filename>', compound_detail, name='compound_detail'),
               path('<int:id>/', compound_detail, name='compound_detail'),
               path('<int:id>/<filename>', compound_detail, name='compound_detail'),
               path('cid/<cid>/<resource>/', compound_detail, name='compound_detail'),
               path('cid/<cid>/<resource>/<filename>', compound_detail, name='compound_detail'),
               path('<int:id>/<resource>/', compound_detail, name='compound_detail'),
               path('<int:id>/<resource>/<filename>', compound_detail, name='compound_detail'),
               ]
