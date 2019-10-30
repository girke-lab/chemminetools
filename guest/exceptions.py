class CannotLogOutGuest(Exception):
    """A guest cannot be logged out. Only deleted."""
    pass

class NotAGuest(Exception):
    """This user is not a guest."""
    pass
