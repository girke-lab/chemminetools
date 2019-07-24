#!/usr/bin/python
# -*- coding: utf-8 -*-

from builtins import object
from django.db import models
from django.contrib.auth.models import User


class Tag(models.Model):
    name = models.CharField(max_length=256)
    user = models.ForeignKey(User, db_index=True)

    def allUserTagNames(user):
        return [tag.name for tag in Tag.objects.filter(user = user)]

class Compound(models.Model):

    cid = models.CharField(max_length=256, db_index=True)
    name = models.CharField(max_length=256)
    formula = models.CharField(max_length=1024)
    weight = models.DecimalField(max_digits=10, decimal_places=2)
    inchi = models.TextField()
    smiles = models.TextField()
    user = models.ForeignKey(User, db_index=True)
    tags = models.ManyToManyField(Tag, db_index=True)

    def byTagNames(tagNames,emptyMeansAll=False):
        tags = Tag.objects.filter(name__in=tagNames)
        if "all" in tags or (emptyMeansAll and len(tags)==0):
            return compound.objects.all()
        else:
            return Compound.objects.filter(tags__in=tags)

    class Meta(object):

        ordering = ['id']

    def __unicode__(self):
        return '%s_%s' % (self.id, self.cid)



class SDFFile(models.Model):

    sdffile = models.TextField()
    compound = models.ForeignKey(Compound)

    def __unicode__(self):
        return '%s' % self.id


