#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from django_cron import CronJobBase, Schedule
from .runapp import deleteOrphanJobs


class deleteOrphans(CronJobBase):

    RUN_EVERY_MINS = 4000

    schedule = Schedule(run_every_mins=RUN_EVERY_MINS)
    code = 'tools.cron.deleteOrphans'

    def do(self):
        deleteOrphanJobs()


