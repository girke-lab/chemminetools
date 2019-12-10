#from django.template.defaulttags import register
from django.utils.html import format_html

from django import template

register = template.Library()

@register.filter
@register.simple_tag
def dict_lookup(d, k):
    return d.get(k)

@register.simple_tag
def url_format(format_string, arg):
    return format_string % arg

@register.simple_tag
def html_format(format_string, arg):
    return format_html(format_string, arg)

@register.filter
def order_by(l, orderByIndex):
    l=list(l)
    l.sort(key=lambda x : x[orderByIndex].lower())
    return l

@register.simple_tag
def find_col_by_key(info, key, value):
    """Given a list of dicts 'info', return the index for the first instance in
    which info[key] = value."""
    for i in range(0, len(info)):
        if info[i].get(key) == value:
            return i
    return None
