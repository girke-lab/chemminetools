"""
Various modified authentication routines to allow for guest users.

These are for users who haven't yet signed up, but for whom temporary
user accounts are created so that they can test the site out.
"""

# Set up settings for this application.
# They default to those in guest.default_settings.py but are
# overridden by DJANGO_SETTINGS_MODULE.
from gyroid_utils.conf import Settings, module_from_env_var
from . import default_settings
settings = Settings([module_from_env_var(), default_settings])
