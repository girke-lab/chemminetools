from django.template.defaulttags import register

@register.filter
def dict_lookup(d, k):
    return d.get(k)
