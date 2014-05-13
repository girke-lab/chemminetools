#!/bin/bash
# These are notes for configuring a chemminetools web server 
# These instructions were tested on Debian Linux 7.0 aka "Wheezy"

# 1st:
# installed custom screen & vim plugins


# Install Required Debian Packages (including django) 
apt-get update
apt-get upgrade -y
apt-get install -y git
apt-get install -y python-django
apt-get install -y postgresql-9.1
apt-get install -y python-psycopg2
apt-get install -y python-skimage
apt-get install -y python-pip
apt-get install -y python-django-reversion
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
apt-get install -y libpq-dev

# install these w/ sid allowed via apt pinning
# apt-get install -y -t=sid libopenbabel4 openbabel
# update: this is commented out because the sid openbabel wasn't working....
#       instead we are now compiling OB from source

pip install django-bootstrap-toolkit
pip install django-cms==2.3.5 south
pip install django-appmedia
pip install ZSI
pip install django-celery
pip install simplejson
pip install ghostscript
pip install pyyaml
pip install django-userena
pip install beautifulsoup4 

# clean up package install 
apt-get clean

# create symbolic link for /srv/chemminetools
ln -s /vagrant /srv/chemminetools

# copy config file
cd /srv/chemminetools
cp chemminetools/settings_sample.py chemminetools/settings.py 
# made changes to get it running: add in database settings and secret key

# create postgresql database as postgres user
sudo -u postgres createuser -U postgres cmt -w -S -R -d
sudo -u postgres psql -U postgres -d postgres -c "alter user cmt with password 'cmt';"
sudo -u postgres createdb -E utf8 -O cmt chemminetools -T template0 --locale=C.UTF-8

# manually install packages in /usr/local/lib/python2.7/dist-packages:
cd /tmp
svn checkout http://django-guest.googlecode.com/svn/trunk/ django-guest-read-only
svn checkout http://django-cron.googlecode.com/svn/trunk/ django-cron-read-only
mv django-guest-read-only/guest /usr/local/lib/python2.7/dist-packages/
mv django-guest-read-only/gyroid_utils /usr/local/lib/python2.7/dist-packages/
mv django-cron-read-only/django_cron /usr/local/lib/python2.7/dist-packages/
rm -rf django-cron-read-only django-guest-read-only

# add sql commands to blank database
cd /srv/chemminetools
mkdir static_production
python manage.py syncdb --noinput
python manage.py migrate --noinput
python manage.py collectstatic --noinput
python manage.py check_permissions

# install R packages
printf "source(\"http://bioconductor.org/biocLite.R\")
biocLite()
biocLite(c(\"ChemmineR\", \"ctc\", \"rjson\", \"R.utils\", \"eiR\", \"RPostgreSQL\"))
" | R --slave

# create working directory and set permissions
cd /srv/chemminetools
mkdir /srv/working
chown www-data /srv/working
ln -s /srv/working working

# register all applications in database
cd /srv/chemminetools/tools/tool_scripts
find *.yaml -print | xargs -I {} ./loader.py -i {}

# append correct settings to apache config
cat /srv/chemminetools/apacheconfig >> /etc/apache2/mods-available/wsgi.conf
# echo "export PYTHONPATH=\"/usr/local/lib/:/usr/local/lib/python2.7/dist-packages/\"" >> /etc/profile

# add apache module
a2enmod wsgi

# daemonize celery
cd /etc/init.d
cp /srv/chemminetools/celeryd .
chmod ugo+x celeryd
update-rc.d celeryd defaults 98 02
cp /srv/chemminetools/celery_config /etc/default/celeryd

# cause celery to launch AFTER /vagrant is mounted
cp /srv/chemminetools/udev_rule /etc/udev/rules.d/50-vagrant-mount.rules

# start celery now
/etc/init.d/celeryd start

# restart apache
/etc/init.d/apache2 restart

# exit now if a bash script, the rest should be run interactively 
exit 0

# The following steps are only necessary for creating a real production environment

# create a user for chemmine tools, and su to that user
cd ~
git clone git@github.com:TylerBackman/chemminetools.git

# as root move to /srv
mv chemminetools /srv/

# to get joelib working:
# -install java jre1.7.0_17 in /opt/jre
# - install JOELib2-alpha-20070303 in /opt/JOELib2-alpha-20070303

# log into /admin and add users and static content

### OPTIONAL FOR DEV SITE: ###
# launch test page on local port
python manage.py runserver 8020

# launch celery worker on local port
python manage.py celery worker --loglevel=info

########################
# Creating a superuser #
########################
python manage.py createsuperuser

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
