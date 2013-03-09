#!/bin/bash
# These are notes for configuring and installing the chemminetools container

# 1st:
# installed custom screen & vim plugins


# Install Django
apt-get install -y git
apt-get install -y python-django
apt-get install -y postgresql-9.1
apt-get install -y python-psycopg2
apt-get install -y python-skimage
apt-get install -y python-pip
apt-get install -y python-django-reversion
apt-get install -y python-openbabel
apt-get install -y libcurl4-openssl-dev
apt-get install -y rabbitmq-server
apt-get install -y python-pylibmc
apt-get install -y libxml2-dev
apt-get install -y openjdk-7-jre
apt-get install -y libreadline5

pip install django-bootstrap-toolkit
pip install django-cms south
pip install django-appmedia
pip install django-cron
pip install ZSI
pip install django-celery
pip install simplejson
pip install ghostscript

# pull down project
su tbackman
cd /srv
git clone git@github.com:TylerBackman/chemminetools.git

# copy config file
scp biocluster.ucr.edu:/home_girkelab/tbackman/Projects/chemminetools/settings.py.devel .
ln -s settings.py.devel settings.py
# made changes to get it running

# create postgresql database
su postgres
createuser -U postgres chemminetools -P
createdb -E utf8 -O chemminetools chemminetools -T template0

# manually install packages in /usr/local/lib/python2.7/dist-packages:
svn checkout http://django-guest.googlecode.com/svn/trunk/ django-guest-read-only
svn checkout http://django-cron.googlecode.com/svn/trunk/ django-cron-read-only
mv django-guest-read-only/guest /usr/local/lib/python2.7/dist-packages/
mv django-guest-read-only/gyroid_utils /usr/local/lib/python2.7/dist-packages/
mv django-cron-read-only/django_cron /usr/local/lib/python2.7/dist-packages/

# add sql commands to blank database
python manage.py syncdb
python manage.py migrate

# exit if you're running this as a bash script
exit 0

# launch test page on local port
python manage.py runserver 8020

# launch celery worker on local port
python manage.py celery worker --loglevel=info

#####################
# Django Console    #
#####################

# test compound upload from django console
python manage.py shell
from myCompounds import views
from compounddb import tools
sdf = open('/home/tbackman/example_db.sdf', 'r')
sdf = sdf.read()
views.addMyCompounds(sdf, 'tbackman')

# get a users compounds
from compounddb.models import Compound
base_queryset = Compound.objects
base_queryset.filter(username='tbackman')

# list all users
from django.contrib.auth.models import User
base_queryset = User.objects
base_queryset.filter()

# test compound clustering and celery 
from tools.runapp import launch
from compounddb import tools
sdf = open('/home/tbackman/example_db.sdf', 'r')
sdf = sdf.read()
result = launch.delay("apcluster.R", "", sdf)
result.ready()
result.get()
result.result
result.successful()

# show results from a job
from tools.models import Job
base_queryset = Job.objects
from tools.runapp import launch
result = launch.AsyncResult(base_queryset.filter()[1].task_id)
result.result

# list all applications
from tools.models import Application
base_queryset = Application.objects
base_queryset.filter()

# make SDF from a users compounds
from compounddb.models import Compound, SDFFile
base_queryset = Compound.objects
compoundList = base_queryset.filter(username='tbackman')
# get sdf:
compound.sdffile_set.all()[0].sdffile

#####################
# adding a new app  #
#####################

python manage.py startapp <newAppName>

# add app to settings.py:
'<newAppName>',

# add app to urls.py:
(r'^<newAppName>/', include('<newAppName>.urls')),

# edit <newAppName>/models.py to add base models

# check new database entries
python manage.py sql <newAppName>

# load new tables into database
python manage.py syncdb

#####################
# Database use      #
#####################
psql chemminetools -U chemminetools -h localhost -W
\d # lists tables
\d+ <tablename> # show schema
drop table <tablename> # show schema
\q # quit
