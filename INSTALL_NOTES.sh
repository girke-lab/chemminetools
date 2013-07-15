#!/bin/bash
# These are notes for configuring a chemminetools web server 
# These instructions were tested on Debian Linux 7.0 aka "Wheezy"

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
apt-get install -y r-base-core
apt-get install -y subversion
apt-get install -y apache2
apt-get install -y libapache2-mod-wsgi
apt-get install -y memcached

pip install django-bootstrap-toolkit
pip install django-cms south
pip install django-appmedia
pip install ZSI
pip install django-celery
pip install simplejson
pip install ghostscript
pip install pyyaml
pip install django-userena

# exit now if a bash script, the rest should be run interactively 
exit 0

# create a user for chemmine tools, and su to that user
cd ~
git clone git@github.com:TylerBackman/chemminetools.git

# as root move to /srv
mv chemminetools /srv/

# rename config file
mv /srv/chemminetools/chemminetools/settings_sample.py /srv/chemminetools/chemminetools/settings.py 
# made changes to get it running: add in database settings and secret key

# create postgresql database as postgres user
su postgres
createuser -U postgres chemminetools -P
createdb -E utf8 -O chemminetools chemminetools -T template0

# as root
# manually install packages in /usr/local/lib/python2.7/dist-packages:
svn checkout http://django-guest.googlecode.com/svn/trunk/ django-guest-read-only
svn checkout http://django-cron.googlecode.com/svn/trunk/ django-cron-read-only
mv django-guest-read-only/guest /usr/local/lib/python2.7/dist-packages/
mv django-guest-read-only/gyroid_utils /usr/local/lib/python2.7/dist-packages/
mv django-cron-read-only/django_cron /usr/local/lib/python2.7/dist-packages/

# unset LC_ALL
export LC_ALL=

# add sql commands to blank database
cd /srv/chemminetools
python manage.py syncdb
python manage.py migrate
python manage.py collectstatic
python manage.py check_permissions

# install R packages
R
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("ChemmineR", "ctc", "rjson", "R.utils", "eiR"))
q()

# create working directory and set permissions
mkdir /srv/chemminetools/working
sudo chown www-data /srv/chemminetools/working

# register all applications in database
cd /srv/chemminetools/tools/tool_scripts
./loader.py -i <appname>.yaml

# setup Apache: 
# add to /etc/apache2/mods-available/wsgi.conf:
    Alias /static/ /srv/chemminetools/static_production/
    WSGIScriptAlias / /srv/chemminetools/chemminetools/wsgi.py
    <Location />
        Order Allow,Deny
        Allow from all
    </Location>
    <Location /admin>
        Order Deny,Allow
        Deny from all
        Allow from .ucr.edu
    </Location>

# add apache module
sudo a2enmod wsgi

# daemonize celery
cd /etc/init.d
wget https://raw.github.com/celery/celery/3.0/extra/generic-init.d/celeryd
chmod ugo+x celeryd
sudo update-rc.d celeryd defaults

# put the following in /etc/default/celeryd:
###############################################
CELERYD_NODES="w1"

# Where to chdir at start.
CELERYD_CHDIR="/srv/chemminetools"

# How to call "manage.py celeryd_multi"
CELERYD_MULTI="$CELERYD_CHDIR/manage.py celeryd_multi"

# How to call "manage.py celeryctl"
CELERYCTL="$CELERYD_CHDIR/manage.py celeryctl"

# Extra arguments to celeryd
CELERYD_OPTS="--time-limit=172800 --concurrency=4"

# %n will be replaced with the nodename.
CELERYD_LOG_FILE="/var/log/celery/%n.log"
CELERYD_PID_FILE="/var/run/celery/%n.pid"
CELERY_CREATE_DIRS=1

# Workers should run as an unprivileged user.
CELERYD_USER="www-data"
CELERYD_GROUP="www-data"

# Name of the projects settings module.
export DJANGO_SETTINGS_MODULE="chemminetools.settings"
############################################################

# to get joelib working:
# -install java jre1.7.0_17 in /opt/jre
# - install JOELib2-alpha-20070303 in /opt/JOELib2-alpha-20070303

# log into /admin and add users and static content

# collect static files:
./manage.py collectstatic

### OPTIONAL FOR DEV SITE: ###
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
drop table <tablename> # delete table 
\q # quit
