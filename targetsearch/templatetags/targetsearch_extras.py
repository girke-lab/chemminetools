#from django.template.defaulttags import register

from django import template

register = template.Library()

@register.filter
def dict_lookup(d, k):
    return d.get(k)

@register.simple_tag
def dict_lookup2(d, k):
    return d.get(k)

def dict_lookup3(d, k):
    return d.get(k)

register.filter(dict_lookup3)
register.simple_tag(dict_lookup3)

@register.simple_tag
def url_format(format_string, arg):
    return format_string % arg
