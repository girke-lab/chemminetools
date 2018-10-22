#!/bin/bash
# These are notes for configuring a chemminetools web server 
# These instructions were tested on Debian Linux 7.0 aka "Wheezy"

# 1st:
# installed custom screen & vim plugins

set -o xtrace

# for debian 9, both python 2.7 and python 3.5 are installed by default. 
# this will make python3 the default for cli 
update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1
update-alternatives --install /usr/bin/python python /usr/bin/python3.5 2

#use python 2.7 for install debian packages
update-alternatives --set python /usr/bin/python2.7

# Install Required Debian Packages (including django) 
apt-get update
apt-get install -y git
apt-get install -y postgresql-9.6
apt-get install -y python3-psycopg2
apt-get install -y python-skimage
apt-get install -y python3-pip
apt-get install -y libcurl4-openssl-dev
apt-get install -y rabbitmq-server
apt-get install -y python-pylibmc
apt-get install -y libxml2-dev
apt-get install -y openjdk-8-jre
apt-get install -y libreadline5
apt-get install -y r-base-core
apt-get install -y subversion
apt-get install -y apache2
apt-get install -y libapache2-mod-wsgi-py3
apt-get install -y memcached
apt-get install -y libpq-dev
apt-get install -y nginx
apt-get install -y libgsl-dev  # needed for eiR
apt-get install -y openbabel libopenbabel-dev swig
apt-get install -y libgc1c2 # for fmcs

# clean up package install 
apt-get clean

#swithc back to python 3 for the rest
update-alternatives --set python /usr/bin/python3.5

# For correct dependency resolution pip needs everything on a single line
pip3 install Django==1.11.14 django-guardian==1.4.9  django-bootstrap-toolkit==2.15.0 \
django-cms==3.5.2 South==1.0.2  django-appmedia==1.0.1 django-celery==3.2.2 simplejson==3.16.0 \
ghostscript PyYAML==3.12 beautifulsoup4==4.6.1 celery==3.1.26.post2 django-userena==2.0.1 \
'html5lib<0.99999999' subprocess32 openbabel django-mptt==0.9.1 djangocms_text_ckeditor==3.6.0 \
djangocms_picture djangocms_link djangocms_file djangocms_googlemap django_cron==0.5.1 \
python-memcached future configparser django-ckeditor

# create symbolic link for /srv/chemminetools
ln -s /vagrant /srv/chemminetools

# copy config file
cd /srv/chemminetools
cp chemminetools/settings_sample.py chemminetools/settings.py 
# made changes to get it running: add in database settings and secret key

# create postgresql database as postgres user
sudo -u postgres createuser -U postgres chemminetools -w -S -R -d
sudo -u postgres psql -U postgres -d postgres -c "alter user chemminetools with password 'chemminetools';"
sudo -u postgres createdb -E utf8 -O chemminetools chemminetools -T template0 --locale=C.UTF-8

# manually install packages in /usr/local/lib/python2.7/dist-packages:
cd /tmp
wget http://biocluster.ucr.edu/~khoran/guest.tgz
wget http://biocluster.ucr.edu/~khoran/gyroid_utils.tgz
tar xfz gyroid_utils.tgz -C /usr/local/lib/python3.5/dist-packages/
tar xfz guest.tgz -C /usr/local/lib/python3.5/dist-packages/
rm -rf guest.tgz gyroid_utils.tgz

# add sql commands to blank database
cd /srv/chemminetools
mkdir static_production
python manage.py syncdb --noinput
python manage.py migrate --noinput
python manage.py collectstatic --noinput
python manage.py check_permissions

# install R packages

printf "install.packages(c(\"amap\",\"bitops\",\"R.oo\",\"gridExtra\"),repos=\"https://cloud.r-project.org\")" | R --slave

printf "source(\"http://bioconductor.org/biocLite.R\")
biocLite()
biocLite(c(\"ChemmineR\", \"ctc\", \"rjson\", \"R.utils\", \"eiR\", \"RPostgreSQL\"),dependencies=c(\"Imports\"))
" | R --slave

# create working directory and set permissions
cd /srv/chemminetools
mkdir /srv/working
chown www-data /srv/working
ln -s /srv/working working

# register all applications in database
cd /srv/chemminetools/tools/tool_scripts
find *.yaml -print | xargs -I {} ./loader.py -i {}

# enable wsgi settings in apache
cp /srv/chemminetools/apacheconfig  /etc/apache2/conf-available/chemminetools.conf
a2enconf chemminetools
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
