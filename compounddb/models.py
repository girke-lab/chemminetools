#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.db import models
from django.contrib.auth.models import User


class Compound(models.Model):

    cid = models.CharField(max_length=256, db_index=True)
    name = models.CharField(max_length=256)
    formula = models.CharField(max_length=1024)
    weight = models.DecimalField(max_digits=10, decimal_places=2)
    inchi = models.TextField()
    smiles = models.TextField()
    user = models.ForeignKey(User, db_index=True)

    class Meta:

        ordering = ['id']

    def __unicode__(self):
        return '%s_%s' % (self.id, self.cid)


class SDFFile(models.Model):

    sdffile = models.TextField()
    compound = models.ForeignKey(Compound)

    def __unicode__(self):
        return '%s' % self.id


