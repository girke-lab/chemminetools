from django_cron import CronJobBase, Schedule

from .utils import cleanup_guests
from . import settings

class DeleteOldGuests(CronJobBase):
    """
        Cron Job that deletes expired Sessions
    """
    RUN_EVERY_MINS = settings.GUEST_DELETE_FREQUENCY
    schedule = Schedule(run_every_mins = RUN_EVERY_MINS)
    code = 'guest.cron.DeleteOldGuests'
    def do(self):
        cleanup_guests()

