from django.db import models

# Create your models here.
#on aws csv
#DrugBank.ID,Name,Type,UniProt.ID,UniProt.Name
#
#DB00001,Lepirudin,BiotechDrug,P00734,Prothrombin

#drugbank_drug_targets_org on virtualbox
# id | target_drugs | uniprot_id | organism
#----+--------------+------------+----------
  # 1 | DP001        | PI001      | Humans
  # 2 | DP002        | PI002      | Humans
  # 3 | DP003        | PI003      | Humans


class Drug_targets_org(models.Model):

    target_drugs = models.TextField()
    uniprot_id = models.TextField()
    organism = models.TextField()

    class Meta:
        managed = True
        app_label = 'drugbank'



    def get_fields(self):
        return [(field.name, field.value_to_string(self)) for field in Drug_targets_org._meta.fields]

    def get_field(self, name):
        return getattr(self, name)

    def __unicode__(self):
        return '%s', (self.target_drugs)


class Drugtargets_org(models.Model):

    drugbank_id = models.TextField()
    drug_name = models.TextField()
    drug_type = models.TextField()
    uniprot_id = models.TextField()
    uniprot_name = models.TextField()

    class Meta:
        managed = False
        app_label = 'drugbank'
        db_table = 'drugtargets_org'



    def get_fields(self):
        return [(field.name, field.value_to_string(self)) for field in Drugtargets_org._meta.fields]

    def get_field(self, name):
        return getattr(self, name)

    def __unicode__(self):
        return '%s', (self.drugbank_id)

class Drugtargets_distinct_org(models.Model):

    drugbank_id = models.TextField()
    drug_name = models.TextField()
    drug_type = models.TextField()
    #uniprot_id = models.TextField()
    #uniprot_name = models.TextField()

    class Meta:
        managed = False
        app_label = 'drugbank'
        db_table = 'drugtargets_distinct_org'



    def get_fields(self):
        return [(field.name, field.value_to_string(self)) for field in Drugtargets_distinct_org._meta.fields]

    def get_field(self, name):
        return getattr(self, name)

    def __unicode__(self):
        return '%s', (self.drugbank_id)


class Ensb_uniport(models.Model):

    ensemble_id = models.TextField()
    uniprot_id = models.TextField()

    class Meta:
        managed = True
        app_label = 'drugbank'

    def get_fields(self):
        return [(field.name, field.value_to_string(self)) for field in Ensb_uniport._meta.fields]

    def get_field(self, name):
        return getattr(self, name)


    def __unicode__(self):
        return '%s' ,self.ensemble_id

#SQL:
# CREATE TABLE drugtargets_org (
#     drugbank_id    text,
#     drug_name      text,
#     drug_type    text,
#     uniprot_id   text,
#     uniprot_name   text
#
# );
# COPY drugtargets_org(drugbank_id, drug_name,drug_type, uniprot_id, uniprot_name) from '/vagrant/drug_target_uniprot.csv'  WITH NULL AS '' DELIMITER ',' CSV;