#!/usr/bin/python
# -*- coding: utf-8 -*-

from builtins import object
from django.db import models
from django.contrib.auth.models import User,AnonymousUser


class Tag(models.Model):
    name = models.CharField(max_length=256)
    user = models.ForeignKey(User, db_index=True,on_delete=models.CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields= ['name','user'],name="unique_name_user_tag") ]

    def allUserTagNames(user):
        if user.id == None:  # user is AnonymousUser
            return set()
        return {tag.name for tag in Tag.objects.filter(user = user)}

    def ensureAllExist(tagNames,user):
        existingTags = Tag.allUserTagNames(user)
        givenTags = set(tagNames)

        for newTag in givenTags.difference(existingTags):
            #print("creating new tag: "+newTag+" for user "+user.username)
            Tag.objects.create(name = newTag, user=user)



class Compound(models.Model):

    cid = models.CharField(max_length=256, db_index=True)
    name = models.CharField(max_length=1024)
    formula = models.CharField(max_length=1024)
    weight = models.DecimalField(max_digits=10, decimal_places=2)
    inchi = models.TextField()
    smiles = models.TextField()
    user = models.ForeignKey(User, db_index=True,on_delete=models.CASCADE)
    tags = models.ManyToManyField(Tag, db_index=True)


    def byTagNames(tagNames,user):
        if "all" in tagNames:
            return Compound.objects.filter(user=user)
        else:
            tags = Tag.objects.filter(name__in=tagNames)
            return Compound.objects.filter(user=user,tags__in=tags).distinct("id")

    class Meta(object):

        ordering = ['id']

    def __unicode__(self):
        return '%s_%s' % (self.id, self.cid)

    def inWorkbench(user, ids):
        workbench = Compound.objects.filter(user=user)
        common_ids = { w.cid for w in workbench } & set(ids)
        return common_ids


class SDFFile(models.Model):

    sdffile = models.TextField()
    compound = models.ForeignKey(Compound,on_delete=models.CASCADE)

    def __unicode__(self):
        return '%s' % self.id


