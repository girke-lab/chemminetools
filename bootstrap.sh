#!/bin/bash
# This script is for configuring a ChemmineTools web server
# It is known to work on Debian 10.2 (Buster)

# In its current form, this script should be run as root on a fresh, minimal,
# updated install of Debian Buster. It is assumed that ChemmineTools has been
# extracted/cloned to /srv/chemminetools, and that this script is running from
# there. It has not been tested with Vagrant yet.

# Echo commands, and exit on error
set -xe

# Install GPG early
apt update
apt install -y gnupg

# Add CRAN's APT repo to sources
echo 'deb https://cloud.r-project.org/bin/linux/debian buster-cran35/' >> /etc/apt/sources.list
# Add GPG key to APT keyring
apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
## See https://cran.r-project.org/bin/linux/debian/ for guidance when this
## inevitably becomes out-of-date.


# Install Debian packages
apt update
## If this is a development installation, these packages may be useful
#apt install -y git man rsync ssh sudo vim wget
## Direct ChemmineTools dependencies
apt install -y apache2 curl libapache2-mod-wsgi-py3 memcached openbabel postgresql \
python3 python3-future python3-memcache python3-pip python3-psycopg2 \
python3-venv rabbitmq-server r-base-core
## Indirect ChemmineTools dependencies (probably to compile pip/R packages)
apt install -y libcurl4-openssl-dev libopenbabel-dev libpq-dev librsvg2-dev \
libssl-dev libxml2-dev libzmq3-dev swig
## We won't need Apache until the very end. Might as well stop it early.
systemctl stop apache2


# Setup Python venv
mkdir -p /srv/venvs
pushd /srv/venvs
python3 -m venv --system-site-packages chemminetools
source chemminetools/bin/activate
popd
## Install Python packages in venv
pip install -r django_req.txt beautifulsoup4 celery html5lib openbabel PyYAML \
requests simplejson
## Manually install gyroid_utils
curl http://cluster.hpcc.ucr.edu/~khoran/gyroid_utils.tgz | tar -xzf- --no-same-owner --no-same-permissions -C /srv/venvs/chemminetools/lib/python3.7/site-packages/


set +x
# Install R packages
for pkg in amap bitops R.oo gridExtra rzmq RPostgreSQL BiocManager shiny; do
    echo "Installing R package [$pkg]. Check /tmp/$pkg.log for errors."
    echo 'install.packages(c("'$pkg'"), repos="https://cloud.r-project.org")' | R --slave &> /tmp/$pkg.log
done
# Install Bioconductor packages
for pkg in ChemmineOB ChemmineR fmcsR ctc rjson R.utils eiR; do
    echo "Installing Bioconductor package [$pkg]. Check /tmp/$pkg.log for errors."
    echo 'BiocManager::install(c("'$pkg'"))' | R --slave &> /tmp/$pkg.log
done
set -x


# Setup PostgreSQL user for ChemmineTools
runuser -u postgres -- createuser -U postgres chemminetools -w -S -R -d
runuser -u postgres -- psql -U postgres -d postgres -c "alter user chemminetools with password 'chemminetools';"
runuser -u postgres -- createdb -E utf8 -O chemminetools chemminetools -T template0 --locale=C.UTF-8

# Setup Celery
install -m 644 celery_config /etc/default/celeryd
install -m 644 celery_tmpfiles.conf /etc/tmpfiles.d/celery.conf
install -m 644 celery.service /etc/systemd/system/celery.service

# Setup RabbitMQ for remote connections
rabbitmqctl add_user remote_worker askfj4l3nbb43
rabbitmqctl set_permissions -p '/' remote_worker ".*" ".*" ".*"

# Setup Apache HTTP Server
install -m 644 apacheconfig /etc/apache2/conf-available/chemminetools.conf
a2enmod wsgi
a2enconf chemminetools

# Create working directory and set permissions
mkdir -p /srv/working
chown www-data:www-data /srv/working
ln -s /srv/working /srv/chemminetools/working

# Copy sample Django settings file and generate a secret key
# Note: you can't use sed's s/// command here because the secret key may
# contain an '&', which will freak sed out...
pushd /srv/chemminetools/chemminetools
SECRET_KEY=$(python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())')
sed -e "/^SECRET_KEY/cSECRET_KEY = '${SECRET_KEY}'" settings_sample.py > settings.py
popd

# Initialize ChemmineTools static files and DB schema
pushd /srv/chemminetools
mkdir static_production
python manage.py migrate --noinput
python manage.py collectstatic --noinput
python manage.py check_permissions
popd

# Register all applications in database
pushd /srv/chemminetools/tools/tool_scripts
for tool in *.yaml; do
    python loader.py -i $tool
done
popd

# Create temp directories for Celery
systemd-tmpfiles --create

# Start remaining services
systemctl daemon-reload
systemctl enable celery
systemctl start celery
systemctl start apache2

# Present user with installation success message and further instructions.
set +x
cat << EOF
ChemmineTools installation *ALMOST* complete!
Before using ChemmineTools, you'll need to make a few more configuration
changes:

- Review /srv/chemminetools/chemminetools/settings.py to make sure the settings
  are appropriate to your installation. Minimally, change the ALLOWED_HOSTS
  list to include your site's domain name(s) and/or IP address. Restart Apache
  after making changes (systemctl restart apache2)
- You'll need to prepare two additional databases, one with ChEMBL data, and
  another with eiR index data. See the following file for more information:
  /srv/chemminetools/chembl_update/README
- Create a Django superuser with:
  $ source /srv/venvs/chemminetools/bin/activate
  $ cd /srv/chemminetools/
  $ python manage.py createsuperuser
- Consider using a firewall to limit access to your server. This will be
  specific to your installation (AWS, VirtualBox, iptables, etc.) and is beyond
  the scope of this script.
EOF
# Exit now if running as a Bash script
exit 0

################
# Useful notes #
################

# Activate venv before doing anything Django-related
source /srv/venvs/chemminetools/bin/activate

# Create a superuser for Django site
python manage.py createsuperuser

# Start a shell with Django environment loaded
python manage.py shell

# Adding a new app
python manage.py startapp newAppName
## Add newAppName to INSTALLED_APPS set in settings.py
## Add app to urlpatterns list in urls.py with something like:
## url(r'^newAppName/', include('newAppName.urls')),

# Create database migrations (changes to DB schema over time)
python manage.py makemigrations newAppName

# Apply database migrations
python manage.py migrate

###############################
# Older notes - for reference #
###############################

# The following steps are only necessary for creating a real production environment

# create a user for chemmine tools, and su to that user
cd ~
git clone git@github.com/girke-lab/chemminetools.git

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
# Database use      #
#####################
psql chemminetools -U chemminetools -h localhost -W
\d # lists tables
\d+ <tablename> # show schema
drop table <tablename> # delete table 
\q # quit
