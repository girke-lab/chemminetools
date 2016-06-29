#! /usr/bin/env python
# $Header$
import sys
try:
    from setuptools import setup
    hasSetuptools = True
except ImportError:
    from distutils.core import setup
    hasSetuptools = False

_url = "http://pywebsvcs.sf.net/"

import ConfigParser
cf = ConfigParser.ConfigParser()
cf.read('setup.cfg')
major = cf.getint('version', 'major')
minor = cf.getint('version', 'minor')
patchlevel = cf.getint('version', 'patchlevel')
candidate = cf.getint('version', 'candidate')

_version = "%d.%d" % ( major, minor )
if patchlevel:
    _version += '.%d' % patchlevel
if candidate:
    _version += '_rc%d' % candidate

try:
    open('ZSI/version.py', 'r').close()
except:
    print 'ZSI/version.py not found; run "make"'
    sys.exit(1)

_packages = [ "ZSI", "ZSI.generate", "ZSI.wstools"]
if sys.version_info[0:2] >= (2, 4):
    _packages.append("ZSI.twisted")
    

# setuptools specific logic
additional_params = {}
if hasSetuptools:
    additional_params['entry_points'] = {
        'console_scripts': [
            'wsdl2py = ZSI.generate.commands:wsdl2py',
            'wsdl2dispatch = ZSI.generate.commands:wsdl2dispatch',
        ],
    }
    additional_params['install_requires'] = [ "PyXML >= 0.8.3", ]
    additional_params['setup_requires'] = [ "setuptools >= 0.6c3", ]
    additional_params['dependency_links'] = [
        "http://sourceforge.net/project/showfiles.php?group_id=6473&package_id=6541&release_id=286213",
    ]
else:
    additional_params['scripts'] = ["scripts/wsdl2py", "scripts/wsdl2dispatch"]

setup(
    name="ZSI",
    version=_version,
    license="Python",
    packages=_packages,
    description="Zolera SOAP Infrastructure",
    author="Rich Salz, et al",
    author_email="rsalz@datapower.com",
    maintainer="Rich Salz, et al",
    maintainer_email="pywebsvcs-talk@lists.sf.net",
    url=_url,
    long_description="For additional information, please see " + _url,
    **additional_params
)
