"""
Contains various functions to aide the creation of guest users.

It is assumed that deleting a user, automatically takes care of deleting all the
other associated objects that you would want to delete with it.  This can
be taken care of using pre-delete signals.

"""

import datetime
import os
import random
import time

from django.contrib.auth.models import User
from django.contrib.auth import authenticate, login
#from django.utils.hashcompat import md5_constructor
from hashlib import md5
from django.conf import settings as django_settings

from .models import Guest
from .exceptions import NotAGuest
from . import settings as guest_settings

# Use the system (hardware-based) random number generator if it exists.
if hasattr(random, 'SystemRandom'):
    randrange = random.SystemRandom().randrange
else:
    randrange = random.randrange
MAX_SESSION_KEY = 18446744073709551616     # 2 << 63

# This function is modified from django.contrib.sessions.base._get_new_session_key
def _get_new_username():
    "Returns a username that isn't being used."
    # The random module is seeded when this Apache child is created.
    # Use settings.SECRET_KEY as added salt.
    try:
        pid = os.getpid()
    except AttributeError:
        # No getpid() in Jython, for example
        pid = 1
    while 1:
        #username = md5_constructor(
        username = md5(
            ("%s%s%s%s" % (randrange(0, MAX_SESSION_KEY), pid, time.time(),
                          django_settings.SECRET_KEY)).encode()).hexdigest()
        try:
            User.objects.get(username=username)
        except User.DoesNotExist:
            break
    return username[:30]

def create_guest():
    """Creates a guest user."""
    user = User.objects.create_user(_get_new_username(),
        guest_settings.GUEST_EMAIL, guest_settings.GUEST_PASSWORD)
    user.save()
    guest = Guest.create_guest(user)
    guest.save()
    return user

def assign_guest(request):
    """Creates a guest user, authenticates it and logs it in."""
    user = create_guest()
    user = authenticate(username=user.username,
                             password=guest_settings.GUEST_PASSWORD)
    login(request, user)

def get_guest(user):
    try:
        return Guest.objects.get(user=user)
    except Guest.DoesNotExist:
        raise NotAGuest()

def is_a_guest(user):
    try:
        get_guest(user)
        return True
    except NotAGuest:
        return False

def display_username(user):
    """Returns a username suitable for display (the guest usernames are not)."""
    if is_a_guest(user):
        return guest_settings.GUEST_USER_NAME
    else:
        return user.username

def make_guest_permanent(request, user, username, email, password):
    """Turns a guest user into a permanent user."""
    get_guest(user).delete()
    user.username = username
    user.email = email
    user.set_password(password)
    user.save()

def delete_guest(user, deletable_classes):
    get_guest(user).delete()
    user.delete()

def cleanup_guests():
    """Delete guest users who have not interacted with the site for at least
    GUEST_DELETE_TIME"""
    cut_off_time = datetime.datetime.now() - guest_settings.GUEST_DELETE_TIME
    old_guests = Guest.objects.filter(last_used__lt = cut_off_time)
    for guest in old_guests:
        guest.user.delete()
        guest.delete()
