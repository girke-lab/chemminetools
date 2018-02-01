#!/bin/bash
# These are notes for configuring a chemminetools web server 
# These instructions were tested on Debian Linux 7.0 aka "Wheezy"

# 1st:
# installed custom screen & vim plugins


# Install Required Debian Packages (including django) 
apt-get update
apt-get install -y git
apt-get install -y postgresql-9.1
apt-get install -y python-psycopg2
apt-get install -y python-skimage
apt-get install -y python-pip
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
apt-get install -y nginx

# clean up package install 
apt-get clean

# install PyXML-0.8.4 (missing from pip)
cd /tmp
wget http://biocluster.ucr.edu/~tbackman/vagrantImages/PyXML-0.8.4.tar.gz
tar xvfz PyXML-0.8.4.tar.gz
cd PyXML-0.8.4
# There is a bug in newer Ubuntu systems that prevents this from building
# Solution: echo '#define HAVE_MEMMOVE 1' >> /usr/include/python2.7/pyconfig.h
python setup.py build
sudo python setup.py install
cd ..
rm -rf PyXML-0.8.4.tar.gz PyXML-0.8.4

# For correct dependency resolution pip needs everything on a single line
pip install Django==1.4.5 django-guardian==1.1.1 ZSI==2.0-rc3 django-bootstrap-toolkit==2.8.0 \
django-cms==2.3.5 South==0.7.5 django-appmedia==1.0.1 django-celery==3.0.11 simplejson==3.1.0 \
ghostscript==0.4.1 PyYAML==3.10 django-userena==1.2.1 beautifulsoup4==4.3.2 celery==3.0.16 subprocess32

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
wget http://biocluster.ucr.edu/~tbackman/vagrantImages/django_cron.tgz
wget http://biocluster.ucr.edu/~tbackman/vagrantImages/django_guest.tgz
wget http://biocluster.ucr.edu/~tbackman/vagrantImages/gyroid_utils.tgz
tar xvfz django_cron.tgz
tar xvfz django_guest.tgz
tar xvfz gyroid_utils.tgz
mv guest /usr/local/lib/python2.7/dist-packages/
mv django_cron /usr/local/lib/python2.7/dist-packages/
mv gyroid_utils /usr/local/lib/python2.7/dist-packages/
rm -rf guest django_cron django_cron.tgz django_guest.tgz gyroid_utils.tgz gyroid_utils

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

# setup and start nginx pubchem proxy
cp /srv/chemminetools/nginx_config /etc/nginx/sites-enabled/default
/etc/init.d/nginx restart

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
