"""
The default settings used by the quizzes application.

These values will be overriden by those in the DJANGO_SETTINGS_MODULE.
To make sure that this happens these settings should always be
accessed through guest.settings.

"""

import datetime

# The username of a guest user for display purposes.
GUEST_USER_NAME = 'Guest'
# A dummy password for guests.
GUEST_PASSWORD = 'whatever'
# A dummy email for guests.
GUEST_EMAIL = ''

# The amount of time after which an unused guest user can be deleted.
GUEST_DELETE_TIME = datetime.timedelta(hours = 72)
# Frequency with which to delete old guests in seconds if using the
# django-cron application to do this.
GUEST_DELETE_FREQUENCY = 86400
