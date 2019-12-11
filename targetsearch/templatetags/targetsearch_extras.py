#from django.template.defaulttags import register

from django import template

register = template.Library()

def dict_lookup(d, k):
    #print("d: "+str(d))
    #print("k: "+str(k))
    return d.get(k)

register.filter(dict_lookup)
register.simple_tag(dict_lookup)

@register.simple_tag
def url_format(format_string, arg):
    return format_string % arg

@register.filter
def order_by(l, orderByIndex):
    l=list(l)
    l.sort(key=lambda x : x[orderByIndex].lower())
    return l
