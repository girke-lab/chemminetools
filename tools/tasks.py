from __future__ import absolute_import, unicode_literals
from celery import shared_task

import subprocess
from tempfile import NamedTemporaryFile
from django.conf import settings
from compounddb.models import Compound,Tag

outputPath = settings.TOOLS_RESULTS
projectDir = settings.PROJECT_DIR

@shared_task
def launch(
    appname,
    commandOptions,
    input,
    job_id,
    user,
    tagNames = [],
    ):

    if input == 'chemical/x-mdl-sdfile':
        input = makeSDF(user,tagNames)
    outputFileName = outputPath + '/job_' + str(job_id)
    command = [projectDir + '/tools/tool_scripts/' + appname,'--outfile=' + outputFileName] + commandOptions
    print('Running: ' + str(command) + '\n')
    runningTask = subprocess.Popen(command, shell=False,
                                   stdin=subprocess.PIPE,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    try:
        print("input type: "+str(type(input)))
        #outs, errs = runningTask.communicate(input, timeout=(60 * 60 * 24 * 7)) # wait a week
        outs, errs = runningTask.communicate(input.encode())#, timeout=(60 * 60 * 24 * 7))  # wait a week
        print('\n outs --', outs, ' errs --', errs)
        return outputFileName
    except Exception as e:
        print("An exception occured while running "+str(command))
        print (e)
        runningTask.kill()
        return False

def makeSDF(user,tagNames):
    compoundList = None
    print(" in makeSDF, given tag names: "+str(tagNames))
    compoundList = Compound.byTagNames(tagNames,user)

    sdf = u''
    for compound in compoundList:
        sdf = sdf + compound.sdffile_set.all()[0].sdffile.rstrip() \
            + '\n'
    return sdf

