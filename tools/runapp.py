#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import time 
import socket
import errno
from celery import task
from tempfile import NamedTemporaryFile
from django.conf import settings
from django.contrib.auth.models import User
import subprocess32 as subprocess
import random
from django.forms import Form, FileField, ModelChoiceField, \
    IntegerField, HiddenInput, CharField, TextInput
from django import forms
from models import *
from compounddb.models import Compound
from types import NoneType
outputPath = settings.TOOLS_RESULTS
projectDir = settings.PROJECT_DIR


def createJob(
    user,
    applicationName,
    optionsList,
    commandOptions,
    input,
    inputvar='',
    async=True,
    ):

    application = Application.objects.get(name=applicationName)
    newJob = Job(
        user=user,
        application=application,
        options=optionsList,
        input=inputvar,
        output='',
        task_id='',
        )
    newJob.save()
    result = launch.delay(application.script, commandOptions, input,
                          newJob.id, user)
    if not async:
        result.wait()
    newJob.task_id = result.id
    newJob.save()
    return newJob


@task()
def launch(
    appname,
    commandOptions,
    input,
    job_id,
    user,
    ):

    if input == 'chemical/x-mdl-sdfile':
        input = makeSDF(user)
    outputFileName = outputPath + '/job_' + str(job_id)
    command = [projectDir + '/tools/tool_scripts/' + appname,'--outfile=' + outputFileName] + commandOptions
    print 'Running: ' + str(command) + '\n'
    runningTask = subprocess.Popen(command, shell=False,
                                   stdin=subprocess.PIPE)
    try:
        outs, errs = runningTask.communicate(input, timeout=(60 * 60 * 24 * 7)) # wait a week
        return outputFileName
    except:
        runningTask.kill()
        return False

def getAppForm(application_id, user):
    fields = {}
    application = Application.objects.get(id=application_id)
    if application.input_type == 'upload':
        fields['File Upload'] = FileField()
    for option in \
        ApplicationOptions.objects.filter(application=application).order_by('id'
            ):
        try: # check for a special fileInput option
            value = ApplicationOptionsList.objects.get(name='fileInput'
                    , category=option)
            getJobList(user)
            fields[option.name] = \
                ModelChoiceField(queryset=Job.objects.filter(application__output_type=value.realName,
                                 user=user,
                                 status=Job.FINISHED).order_by('-id'),
                                 empty_label='None', required=False)
            continue
        except:
            pass
        try: # check for a special stringInput option
            value = ApplicationOptionsList.objects.get(name='stringInput', category=option)
            fields[option.name] = CharField(max_length=2056, widget=forms.Textarea(attrs={'cols': 200, 'rows': 1,}))
            continue
        except:
            pass
        fields[option.name] = \
            ModelChoiceField(queryset=ApplicationOptionsList.objects.filter(category=option).order_by('id'
                             ), empty_label=None)
    fields['application'] = IntegerField(initial=application.id,
            widget=HiddenInput())
    return type('%sForm' % str(application.name), (Form, ), fields)

def parseToolForm(form):

    # parses a form created by getAppForm to return command line options

    commandOptions = []
    optionsList = u''
    application = Application.objects.get(id=form.cleaned_data['application'])
    for question in form.cleaned_data.keys():
        if question != 'application' and question != 'File Upload':
            questionObject = \
                ApplicationOptions.objects.get(application=application,
                    name=question)
            try:
                job = form.cleaned_data[question]
                option = job.output
                optionName = str(job)
            except:
                try:
                    answerObject = form.cleaned_data[question]
                    optionName = answerObject.name
                    option = answerObject.realName
                except:
                    answerObject = form.cleaned_data[question]
                    optionName = answerObject
                    option = answerObject 
            if option == None:
                option = u'None'
            if optionName == None:
                optionName = u'None'
            commandOptions = commandOptions + ['--'+questionObject.realName + '=' + option]
            optionsList = optionsList + questionObject.name + ': ' \
                + optionName + ', '
    optionsList = re.sub(",\s$", '', optionsList, count=0)
    return commandOptions, optionsList

def getJobList(user):

    # get all running jobs and update those that need updating

    try:
        runningJobs = Job.objects.filter(user=user, status=Job.RUNNING)
        for job in runningJobs:
            updateJob(user, job.id)
    except:
        return False
    return Job.objects.filter(user=user).order_by('-start_time')


def updateJob(user, job_id):

    # checks if a job is done, updates it's status, and then returns it

    # wait in case another process is still creating the new job
    for i in range(5):
        try:
            job = Job.objects.get(id=job_id, user=user)
        except Job.DoesNotExist:
            if i == 4:
                return False
            time.sleep(1)
            continue
        break
    if job.status == Job.RUNNING:
        result = launch.AsyncResult(job.task_id)
        if result.ready():
            if result.result:
                job.status = Job.FINISHED
            else:
                job.status = Job.FAILED
            job.output = result.result
            result.forget()
            job.save()
    return job


def deleteJob(user, job_id):
    job = updateJob(user, job_id)
    if isinstance(job.output, unicode):
        if os.path.isfile(job.output):
            os.remove(job.output)
    try:
        result = launch.AsyncResult(job.task_id)
        result.forget()
    except:
        pass
    job.delete()
    return True


def deleteApp(name):
    try:
        app = Application.objects.get(name=name)
    except:
        return False
    jobs = Job.objects.filter(application=app)
    for job in jobs:
        deleteJob(job.user, job.id)
    app.delete()
    return True


def deleteOrphanJobs():
    try:
        jobs = Job.objects.filter(user__isnull=True)
        for job in jobs:
            deleteJob(job.user, job.id)
        return True
    except:
        return False


def makeSDF(user):
    compoundList = Compound.objects.filter(user=user)
    sdf = u''
    for compound in compoundList:
        sdf = sdf + compound.sdffile_set.all()[0].sdffile.rstrip() \
            + '\n'
    return sdf


