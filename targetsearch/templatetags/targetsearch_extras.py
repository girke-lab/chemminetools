from django.conf import settings
#from django.template.defaulttags import register
from django.utils.html import format_html, mark_safe

from django import template

register = template.Library()

@register.filter
@register.simple_tag
def dict_lookup(d, k):
    #print("d: "+str(d))
    #print("k: "+str(k))
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
def find_col_by_key(info, key, value_prefix, value, default=None):
    """Given a list of dicts 'info', return the index for the first instance in
    which info[key] == value_prefix + value."""
    #print("find_col_by_key key: "+str(key)+", value: "+str(value_prefix+value));
    if info != None:
        for i in range(0, len(info)):
            if info[i].get(key) == value_prefix+value:
                return i
    return default

@register.simple_tag
def col_index_list(info, key, value):
    """Given a list of dicts 'info', return a list of indices corresponding to
    columns in which info[key] == value. Use to build lists of default columns,
    non-exportable columns, etc."""
    index_list = list()
    if info != None:
        for i in range(0, len(info)):
            if info[i].get(key) == value:
                index_list.append(i)
    return index_list

#TODO: Move this to its own app after local_settings.py is implemented
@register.simple_tag
def analytics():
    try:
        return settings.ANALYTICS
    except AttributeError as e:
        return True

@register.simple_tag
def extanno_format(extanno_dict):
    chembl_id = extanno_dict['chembl_id']
    output_html = ''

    if extanno_dict['extannos'] is None:
        return 'None'

    for db_name, display_name in [('drugbank', 'Drugbank'), ('drugage', 'Drugage'), ('lincs', 'LINCS'), ('cmap', 'CMAP')]:
        if db_name in extanno_dict['extannos']:
            output_html += format_html('<a class="btn btn-primary" role="button" target="_blank" href="/targetsearch/extannobychembl/{}/{}/">{}</a> ', chembl_id, db_name, display_name)
    return mark_safe(output_html)
