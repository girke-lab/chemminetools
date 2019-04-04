#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django_cron import CronJobBase, Schedule
from .runapp import (deleteOrphanJobs,deleteJob)
from .models import Job
from datetime import (datetime,timedelta)
from django.contrib.auth.models import User
import traceback
import os


class deleteOrphans(CronJobBase):

    RUN_EVERY_MINS = 4000

    schedule = Schedule(run_every_mins=RUN_EVERY_MINS)
    code = 'tools.cron.deleteOrphans'

    def do(self):
        deleteOrphanJobs()


class deleteOldJobs(CronJobBase):
    RUN_EVERY_MINS = 4000
    schedule = Schedule(run_every_mins=RUN_EVERY_MINS)
    code = 'tools.cron.deleteOldJobs'

    def do(self):
        try:
            #print("checking for old jobs")
            how_many_days = 365
            oldJobs = Job.objects.filter(start_time__lte=datetime.now()-timedelta(days=how_many_days))
            for thisjob in oldJobs:
                #print("deleting job for user "+str(thisjob.user)+", job id: "+str(thisjob.id)+", started on "+str(thisjob.start_time)+
                        #", file: "+str(thisjob.output))
                try:
                    deleteJob(thisjob.user, thisjob.id)
                except OSError:
                    print("failed to remove output file for job "+str(thisjob.id)+": "+str(thisjob.output))
        except:
            print("error in tools.cron.deleteOldJobs")
            print(traceback.format_exc())

class deleteOldUsers(CronJobBase):
    #print("checking for old users")
    RUN_EVERY_MINS = 4000
    schedule = Schedule(run_every_mins=RUN_EVERY_MINS)
    code = 'tools.cron.deleteOldUsers'

    def do(self):
        try:
            how_many_days = 365
            oldUsers = User.objects.filter(
                    email='', 
                    is_staff=False,
                    is_active=False, 
                    is_superuser=False, 
                    last_login__lte=datetime.now()-timedelta(days=how_many_days))

            for user in oldUsers:
                print("deleting user "+str(user.username)+", joined "+str(user.date_joined)+", last login: "+str(user.last_login)+", email: "+str(user.email))
                for userJob in Job.objects.filter(user_id=user.id):
                    deleteJob(user,userJob.id)
                user.delete()
        except:
            print("error in tools.cron.deleteOldJusers")
            print(traceback.format_exc())


