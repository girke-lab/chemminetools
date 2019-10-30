import datetime

from django.db import models
from django.contrib.auth.models import User
from django.contrib.sessions.models import Session

class Guest(models.Model):
    
    """
    A temporary user.

    Fields:
    ``user`` - The temporary user.
    ``last_used`` - The last time we noted this user doing something.
    
    All users with a record in this model are temporary and should be
    deleted after GUEST_DELETE_TIME.
    
    """

    user = models.ForeignKey(User,on_delete=models.CASCADE)
    last_used = models.DateTimeField(User)

    @classmethod
    def create_guest(self, user):
        guest = Guest(user=user, last_used=datetime.datetime.now())
        return guest
