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
apt-get install -y python3-psycopg2
apt-get install -y python-skimage
apt-get install -y python3-pip
apt-get install -y libcurl4-openssl-dev
apt-get install -y python-pylibmc
apt-get install -y libxml2-dev
apt-get install -y libreadline5
apt-get install -y openjdk-8-jre-headless
apt-get install -y r-base-core
apt-get install -y memcached
apt-get install -y libpq-dev
apt-get install -y libgsl-dev  # needed for eiR
apt-get install -y libssl-dev  # needed for eiR
apt-get install -y librsvg2-dev  # needed for ChemmineR
apt-get install -y openbabel libopenbabel-dev swig
apt-get install -y libgc1c2 # for fmcs
apt-get install -y libzmq3-dev #for rzmq server used by eiSearch

# clean up package install 
apt-get clean

#swithc back to python 3 for the rest
update-alternatives --set python /usr/bin/python3.5

# For correct dependency resolution pip needs everything on a single line
pip3 install Django==1.11.14 django-guardian==1.4.9  django-bootstrap-toolkit==2.15.0 \
django-cms==3.5.2 South==1.0.2  django-appmedia==1.0.1 django-celery==3.2.2 simplejson==3.16.0 \
ghostscript PyYAML==3.12 beautifulsoup4==4.6.1 celery django-userena==2.0.1 \
'html5lib<0.99999999' subprocess32 openbabel django-mptt djangocms_text_ckeditor==3.6.0 \
djangocms_picture djangocms_link djangocms_file djangocms_googlemap django_cron==0.5.1 \
python-memcached future configparser django-ckeditor 


# copy config file
cd /srv/chemminetools
cp chemminetools/settings_sample.py chemminetools/settings.py 
# made changes to get it running: add in database settings and secret key

# manually install packages in /usr/local/lib/python2.7/dist-packages:
cd /tmp
wget http://cluster.hpcc.ucr.edu/~khoran/guest.tgz
wget http://cluster.hpcc.ucr.edu/~khoran/gyroid_utils.tgz
tar xfz gyroid_utils.tgz -C /usr/local/lib/python3.5/dist-packages/
tar xfz guest.tgz -C /usr/local/lib/python3.5/dist-packages/
rm -rf guest.tgz gyroid_utils.tgz

# install R packages

printf "install.packages(c(\"amap\",\"bitops\",\"R.oo\",\"gridExtra\",\"rzmq\",\"RPostgreSQL\",\"BiocManager\"),repos=\"https://cloud.r-project.org\")" | R --slave

#printf "source(\"http://bioconductor.org/biocLite.R\")
#biocLite()
#printf "BiocManager::install(c(\"ChemmineR\", \"fmcsR\",\"ChemmineOB\",\"ctc\", \"rjson\", \"R.utils\", \"eiR\", \"RPostgreSQL\"),dependencies=c(\"Imports\"))
printf "BiocManager::install(c(\"ChemmineR\", \"fmcsR\",\"ChemmineOB\",\"ctc\", \"rjson\", \"R.utils\", \"eiR\", \"RPostgreSQL\"))
" | R --slave

# create working directory and set permissions
cd /srv/chemminetools
mkdir /srv/working
chown www-data /srv/working
ln -s /srv/working working


# daemonize celery
cd /etc/init.d
cp /srv/chemminetools/celeryd .
chmod ugo+x celeryd
update-rc.d celeryd defaults 98 02
cp /srv/chemminetools/celery_config /etc/default/celeryd

# start celery now
systemctl start celeryd

