# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Drug_targets_org'
        db.create_table('drugbank_drug_targets_org', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('target_drugs', self.gf('django.db.models.fields.TextField')()),
            ('uniprot_id', self.gf('django.db.models.fields.TextField')()),
            ('Organism', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal('drugbank', ['Drug_targets_org'])

        # Adding model 'Ensb_uniport'
        db.create_table('drugbank_ensb_uniport', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ensemble_id', self.gf('django.db.models.fields.TextField')()),
            ('uniprot_id', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal('drugbank', ['Ensb_uniport'])


    def backwards(self, orm):
        # Deleting model 'Drug_targets_org'
        db.delete_table('drugbank_drug_targets_org')

        # Deleting model 'Ensb_uniport'
        db.delete_table('drugbank_ensb_uniport')


    models = {
        'drugbank.drug_targets_org': {
            'Meta': {'object_name': 'Drug_targets_org'},
            'Organism': ('django.db.models.fields.TextField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'target_drugs': ('django.db.models.fields.TextField', [], {}),
            'uniprot_id': ('django.db.models.fields.TextField', [], {})
        },
        'drugbank.ensb_uniport': {
            'Meta': {'object_name': 'Ensb_uniport'},
            'ensemble_id': ('django.db.models.fields.TextField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uniprot_id': ('django.db.models.fields.TextField', [], {})
        }
    }

    complete_apps = ['drugbank']