#!/usr/bin/python
# -*- coding: utf-8 -*-

from django_cron import cronScheduler, Job
from runapp import deleteOrphanJobs


class deleteOrphans(Job):

    run_every = 4000

    def job(self):
        deleteOrphanJobs()


cronScheduler.register(deleteOrphans)
