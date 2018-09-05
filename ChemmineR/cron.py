#!/usr/bin/python
# -*- coding: utf-8 -*-

from django_cron import CronJobBase, Schedule
from datetime import datetime, timedelta
from django.contrib.auth.models import User
from tools.runapp import deleteJob
from tools.models import Job as toolJob

class deleteOldChemmmineR(CronJobBase):

    RUN_EVERY_MINS = 4000

    schedule = Schedule(run_every_mins = RUN_EVERY_MINS)
    code = 'ChemmineR.cron.deleteOldChemmineR'

    def do(self):
        how_many_days = 60
        user = User.objects.get(username='ChemmineR')
        oldJobs = toolJob.objects.filter(user=user,start_time__lte=datetime.now()-timedelta(days=how_many_days))

        for thisjob in oldJobs:
            deleteJob(thisjob.user, thisjob.id)

