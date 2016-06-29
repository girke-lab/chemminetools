#

# Change this to point to your Python2 exectuable
PYTHON = python

# This should be the Python that generated RPMs depend on and install
# into; this may be different from the Python used to build the RPM.
RPM_TARGET_PYTHON=/usr/bin/python


all:	rpm

ins:
	soup $(PYTHON) setup.py install

kit:	doc rpm	# make a kit

rpm:	ZSI/version.py
	rm -f dist/*
	$(PYTHON) setup.py bdist_rpm --python=$(RPM_TARGET_PYTHON)

doc:	doc/version.tex		# build the docs
	$(MAKE) -C doc

ver:			# update the build number
	$(PYTHON) newver.py --incr

.PHONY: all ins kit rpm doc ver

ZSI/version.py doc/version.tex:	setup.cfg
	$(PYTHON) newver.py
