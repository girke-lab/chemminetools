#!/usr/bin/python
# -*- coding: utf-8 -*-

from yaml import load, Loader, Dumper

# make it possible to run as standalone program

import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/tools')  # allow tools.py to be imported
import os
from django.core.management import setup_environ
import chemminetools.settings
import argparse
setup_environ(chemminetools.settings)

from tools.models import *
from tools.runapp import deleteApp

parser = argparse.ArgumentParser(description='Register tool in database'
                                 )
parser.add_argument('-i', '--infile', help='input YAML file',
                    required=True)
args = vars(parser.parse_args())


def main():
    f = open(args['infile'], 'r')
    data = load(f, Loader=Loader)
    f.close()
    try:
        category = \
            ApplicationCategories.objects.get(name=data['category'])
    except:
        category = ApplicationCategories(name=data['category'])
        category.save()

    # if app exists, delete it before re-adding it

    deleteApp(name=data['name'])

    newApp = Application(
        name=data['name'],
        script=data['script'],
        input_type=data['input_type'],
        output_type=data['output_type'],
        description=data['description'],
        category=category,
        )
    newApp.save()
    try:
        for option in data['ApplicationOptions'].keys():
            optionData = data['ApplicationOptions'][option]
            newOption = ApplicationOptions(name=option,
                    realName=optionData['realName'], application=newApp)
            newOption.save()
            for optionValue in optionData['options']:
                newValue = ApplicationOptionsList(category=newOption,
                        name=optionValue[0], realName=optionValue[1])
                newValue.save()
    except:
        pass
    print 'Registered application ' + newApp.name + '\n'
    return True


if __name__ == '__main__':
    main()
