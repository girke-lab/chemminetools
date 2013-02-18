from django.db import models
from django.contrib.auth.models import User

class Compound(models.Model):
        cid = models.CharField(max_length=256)
        name = models.CharField(max_length=256)
        formula = models.CharField(max_length=256)
        weight = models.DecimalField(max_digits=10, decimal_places=2)
        inchi = models.TextField()
        smiles = models.CharField(max_length=1024)
        username = models.CharField(max_length=30, default='', blank=True)

        class Meta:
                #ordering = ['id']
                #get_latest_by = 'pub_date'
                pass

        def __unicode__(self):
                return "%s_%s" % (self.id, self.cid)

        @models.permalink
        def get_absolute_url(self):
                return ('compound_detail', (),
                        dict( library=self.library,
                                        cid=self.header.cid ) )

class SDFFile(models.Model):
        sdffile = models.TextField()
        compound = models.ForeignKey(Compound)

        class Meta:
                #ordering = ['id']
                pass

        def __unicode__(self):
                return "%s" % self.id
