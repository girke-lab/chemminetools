try:
    from functools import wraps
except ImportError:
    from django.utils.functional import wraps  # Python 2.3, 2.4 fallback.

from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth import REDIRECT_FIELD_NAME

from .utils import assign_guest, is_a_guest

def guest_allowed(view_func):
    """
    Decorator for views where if the user is not logged in, then a temporary
    guest account is created.
    """
    @wraps(view_func)
    def wrapper(request, *args, **kwargs):
        if not request.user.is_authenticated:
            assign_guest(request)
        return view_func(request, *args, **kwargs)
    return wrapper

# Modified from django.contrib.auth.decorators.login_required
def login_required(function=None, redirect_field_name=REDIRECT_FIELD_NAME):
    """
    Decorator for views that checks that the user is logged in, redirecting
    to the log-in page if necessary.  Does not allow guest users.
    """
    actual_decorator = user_passes_test(
        lambda u: u.is_authenticated and not is_a_guest(u),
        redirect_field_name=redirect_field_name
    )
    if function:
        return actual_decorator(function)
    
