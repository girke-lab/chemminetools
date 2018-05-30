# Create your views here.
from django.shortcuts import redirect, render_to_response, render
from django.http  import HttpResponse
from django.template import RequestContext
from guest.decorators import guest_allowed, login_required
from drugbank.models import Drugtargets_distinct_org,Drugtargets_org
import json
import operator
from django.db.models import Q


def drugbank_lookup(request):
    dict_page = {}

    try:
        drugs = request.GET.get('druglist')
    except:
        drugs = request.POST.get('druglist')
    druglist = []
    if drugs:
        str(drugs).split(",")
        druglist = drugs.split(",")

    if druglist:
        #queryset = Drugtargets_distinct_org.objects.filter(drugbank_id__in=druglist)
        queryset = Drugtargets_org.objects.filter(drugbank_id__in=druglist)

        results = get_results(queryset, 'drugbank_id', request)
    else:
        #queryset = Drugtargets_distinct_org.objects.all()
        queryset = Drugtargets_org.objects.all()
        results = get_results(queryset, 'drugbank_id', request)

    #showFields = ['drugbank_id','drug_name','drug_type']#,'uniprot_id','uniprot_name']
    showFields = ['drugbank_id','drug_name','drug_type','uniprot_id','uniprot_name']


    #dict_page['showFields'] = showFields
    url = '/drugbank/?isajax=true&druglist=' + drugs#DB02512,DB02511'

    c = {'showFields':showFields, 'url':url}

    if 'isajax' in request.REQUEST:
        print '\n in ajax in process_table_data_for_json'#, ' and total is ', results['total']
        dict_page['total'] = results['total']
        dict_page['rows'] = get_rows(showFields, results['results'], int(results['page_number']), int(results['rows_number']))
        print '\n this is rows for json', dict_page["rows"]
        result_json = json.dumps(dict_page)
        print '\n result_json is ', result_json

       # result_json =  {"rows": [{"uniprot_id": "PI001", "organism": "Humans", "target_drugs": "DP001"}, {"uniprot_id": "PI002", "organism": "Humans", "target_drugs": "DP002"}, {"uniprot_id": "PI003", "organism": "Humans", "target_drugs": "DP003"}], "total": 3, "theFields": ["id", "target_drugs", "uniprot_id", "organism"]}

        # result_json = {"rows": [{"uniprot_id": "PI001", "organism": "Humans", "target_drugs": "DP001"},
        #           {"uniprot_id": "PI002", "organism": "Humans", "target_drugs": "DP002"},
        #           {"uniprot_id": "PI003", "organism": "Humans", "target_drugs": "DP003"}], "total": 3,
        #  "theFields": ["id", "target_drugs", "uniprot_id", "organism"]}
        #
        #

        return HttpResponse(result_json)
    return render_to_response('annotation_5.html',c, context_instance=RequestContext(request) )


def get_results(queryset, sorting_order, request):
    """ Gets data based on the search type and passed parameters

    @param queryset: result query set
    @param sorting_order: string
    @param request: http request
    @return: results

    """

    results = {}
    page_number = 1
    ordering = 'asc'
    rows_number = 10
    offset_number = (page_number - 1) * rows_number
    # sorting
    print '\n in sorting order'
    sorting_order = get_sorting_order(request, sorting_order)
    print '\n finished sorting ', sorting_order

    # filter
    try:
        results['results'] = queryset
    #     results['results'] = self.get_all_filters(request.REQUEST['filterRules'],queryset, sorting_order)
    #     print '\n got all_filters ', results['results']
    except:
        #     results['results'] = queryset
        print '\n did not get queryset '
    # paging
    try:
        page_number = request.REQUEST['page']
        ordering = request.REQUEST['order']
        rows_number = request.REQUEST['rows']
        offset_number = (page_number - 1) * rows_number
    except:
        pass
    results['results'] = results['results'].order_by(sorting_order)
    results['ordering'] = ordering
    results['rows_number'] = rows_number
    results['page_number'] = page_number
    results['sorting_order'] = sorting_order
    results['offset_number'] = offset_number
    results['total'] = results['results'].count()
    print '\n what is total? ', results['total']

    return results

def get_results_page(results, page_number, rows_number):

        offset_number = int((int(page_number) - 1)) * int(rows_number)
        beg = int(offset_number)
        end = int(rows_number) + beg
        print '\n\n\n get beginning and end of the page ', beg, '->', end
        results_page = results[beg:end]
        return results_page

def get_rows(showFields, results, page_number, rows_number):

        r_list = []
        results_page = get_results_page(results, page_number, rows_number)
        print '\n rows number, page_number ', rows_number, page_number

        n = 0
        for r in results_page:
            n = n + 1
            r_dict = {}

            for f in showFields:
                print '\n field f is ',
                if str(f) != 'id':  # todo:add more fields here
                    r_dict[f] = r.get_field(f)
            #r_dict['structure'] = '<img src="https://www.drugbank.ca/structures/DB02512/image.png" width="160" height="160">'
            drugname = r_dict['drugbank_id']
            #drugname =  'DB02512'
            r_dict['structure'] = '<img src="https://www.drugbank.ca/structures/' + drugname + '/image.png"  width="160" height="160">'
            print '\n\n\n\n structure is for  ', drugname
            r_list.append(r_dict)
        print '\n list is ', r_list
        return r_list

def get_range(string, separator):
        return [x.strip() for x in string.split(separator)]

def get_all_filters(filtering, queryset, sorting_order):
        kwargs = {}
        result = queryset
        try:
            filters = json.loads(filtering)

            print '\n filters are now ', filters
            for filter in filters:
                print '\n for filter in filters ', filter

                argument_list = []
                if str(filter["op"]) == 'similar':
                    value_list = get_range(str(filter["value"]), ",")
                    for l in value_list:
                        argument_list.append(Q(**{str(filter['field']) + '__icontains': l}))

                elif str(filter["op"]) == 'equal':
                    # kwargs[str(filter['field'])] = str(filter["value"])
                    value_list = get_range(str(filter["value"]), ",")
                    for l in value_list:
                        print '\n in l is (', l, ')'
                        argument_list.append(Q(**{str(filter['field']) + '__iexact': l}))

                elif str(filter["op"]) == 'notequal':
                    argument_list.append(~Q(**{str(filter['field']) + '__iexact': filter["value"]}))

                elif str(filter["op"]) == 'between':
                    range = get_range(str(filter["value"]), ":")
                    kwargs = {"%s__range" % (str(filter['field'])): (str(range[0]), str(range[1]))}

                elif str(filter["op"]) == 'less':
                    kwargs = {"%s__lt" % (str(filter['field'])): (filter["value"])}
                elif str(filter["op"]) == 'greater':
                    kwargs = {"%s__gt" % (str(filter['field'])): (filter["value"])}
                elif str(filter["op"]) == 'nofilter':
                    try:
                        del kwargs[str(filter['field'])]
                    except:
                        try:
                            del argument_list[str(filter['field'])]
                        except:
                            print '\n was not able to remove filter'

                result = result.filter(**kwargs).order_by(sorting_order)
                print '\n before argument_list ', argument_list
                result = result.filter(reduce(operator.or_, argument_list))

        except:
            pass
        return result


def get_sorting_order(request,sorting_order):
        ordering = 'asc'

        try:
            sorting_order = request.REQUEST['sort']
            ordering = request.REQUEST['order']
            print '\n try sorting_order is ', sorting_order
            print '\n try ordering is ', ordering

        except:
            pass

        if ordering == 'desc':
            sorting_order = '-'+str(sorting_order)
        print '\n sorting_order is ', sorting_order
        return sorting_order



#=======================================================
def drugbank_lookup_1(request):
    dict_page = {}  # for json and ajax


    dict_page['theFields'] = ['DID','PID','organism' ]
    dict_page['rows'] = ['hello1', 'hello2', 'hello3']
    dict_page['total'] = 5
    result_json = json.dumps(dict_page)
    print 'result_json ', result_json

    url = '/drugbank/?isajax=true/' # not sure
    print '\n request is ', request.REQUEST
    if 'isajax' in request.REQUEST:
        print '\n drugbank urls and is ajax ', url
    else:
        print '\n not in request, ', url, '==', request.REQUEST
    if 'isajax' in request.REQUEST:
        print '\n\n\n\n\n\n\n +++++++ in AJAX Drugbank', request.REQUEST
        # +++++++ in AJAX MALI with checking filter {u'sort': u'id', u'startdate': u'', u'rows': u'20', u'enddate': u'', u'searchtype': u'emergence', u'order': u'asc', u'displaytype': u'table', u'location': u'mali', u'isajax': u'true', u'amp': u'', u'showFields': u'village', u'filterRules': u'[{"field":"morphosp_name","op":"notequal","value":"hello"}]', u'page': u'1'}


        dict_page['rows'] =['hello1','hello2','hello3']
        dict_page['total'] = 5
        result_json = json.dumps(dict_page)
        print '\n do resutl_json', result_json
        return HttpResponse(result_json)

    return render_to_response('annotation_3.html', {'result': result_json, 'url': url}, context_instance=RequestContext(request))
                        #      'annotation_test.html',{'result': result_json, 'url': url})
   # return render_to_response('annotation_test.html', {},
                                   #context_instance=RequestContext(request))


#
# #from django.http import JsonResponse
# from django.db.models import Q
# from .models import Drug_targets_org
#
# import math
#
#
# def drugbank(request):
#     if request.method == 'POST':
#         pagination_content = ""
#         page_number = request.POST['data[page]'] if request.POST['data[page]'] else 1
#         page = int(page_number)
#         name = request.POST['data[th_name]']
#         sort = '-' if request.POST['data[th_sort]'] == 'DESC' else ''
#         search = request.POST['data[search]']
#         max = int(request.POST['data[max]'])
#
#         cur_page = page
#         page -= 1
#         per_page = max  # Set the number of results to display
#         start = page * per_page
#
#         # If search keyword is not empty, we include a query for searching
#         # the "content" or "name" fields in the database for any matched strings.
#         if search:
#             all_posts = Drug_targets_org.objects.filter(Q(content__contains=search) | Q(name__contains=search)).order_by(
#                 sort + name)[start:per_page]
#             count = Drug_targets_org.objects.filter(Q(content__contains=search) | Q(name__contains=search)).count()
#
#         else:
#             all_posts = Drug_targets_org.objects.order_by(sort + name)[start:cur_page * max]
#             count = Drug_targets_org.objects.count()
#
#         if all_posts:
#             for post in all_posts:
#                 pagination_content += '''
# 					<tr>
# 						<td><img src='%s' width='100' /></td>
# 						<td>%s</td>
# 						<td>$%s</td>
# 						<td>%s</td>
# 						<td>%s</td>
# 						<td>%s</td>
# 						<td>
# 							<a href='#' class='text-success'>
# 								<span class='glyphicon glyphicon-pencil' title='Edit'></span>
# 							</a>
# 							<a href='#' class='text-danger delete-product' item_id='%s'>
# 								<span class='glyphicon glyphicon-remove' title='Delete'></span>
# 							</a>
# 						</td>
# 					</tr>
# 				''' % (post.featured_image, post.name, post.price, post.status, post.date, post.quantity, post.id)
#         else:
#             pagination_content += "<tr><td colspan='7' class='bg-danger p-d'>No results</td></tr>"
#
#         # return JsonResponse({
#         #     'content': pagination_content,
#         #     'navigation': nagivation_list(count, per_page, cur_page)
#         # })
#     else:
#         return render(request, 'annotation_1.html')
#
#
# def nagivation_list(count, per_page, cur_page):
#     pagination_nav = ""
#     previous_btn = True
#     next_btn = True
#     first_btn = True
#     last_btn = True
#
#     no_of_paginations = int(math.ceil(count / per_page))
#
#     if cur_page >= 7:
#         start_loop = cur_page - 3
#         if no_of_paginations > cur_page + 3:
#             end_loop = cur_page + 3
#         elif cur_page <= no_of_paginations and cur_page > no_of_paginations - 6:
#             start_loop = no_of_paginations - 6
#             end_loop = no_of_paginations
#         else:
#             end_loop = no_of_paginations
#     else:
#         start_loop = 1
#         if no_of_paginations > 7:
#             end_loop = 7
#         else:
#             end_loop = no_of_paginations
#
#     # Pagination Buttons logic
#     pagination_nav += "<div class='pagination-container'><ul>"
#
#     if first_btn and cur_page > 1:
#         pagination_nav += "<li p='1' class='active'>First</li>"
#     elif first_btn:
#         pagination_nav += "<li p='1' class='inactive'>First</li>"
#
#     if previous_btn and cur_page > 1:
#         pre = cur_page - 1
#         pagination_nav += "<li p='" + str(pre) + "' class='active'>Previous</li>"
#     elif previous_btn:
#         pagination_nav += "<li class='inactive'>Previous</li>"
#
#     for i in range(start_loop, end_loop + 1):
#         if cur_page == i:
#             pagination_nav += "<li p='" + str(i) + "' class = 'selected'>" + str(i) + "</li>"
#         else:
#             pagination_nav += "<li p='" + str(i) + "' class='active'>" + str(i) + "</li>"
#
#     if next_btn and cur_page < no_of_paginations:
#         nex = cur_page + 1
#         pagination_nav += "<li p='" + str(nex) + "' class='active'>Next</li>"
#     elif next_btn:
#         pagination_nav += "<li class='inactive'>Next</li>"
#
#     if last_btn and cur_page < no_of_paginations:
#         pagination_nav += "<li p='" + str(no_of_paginations) + "' class='active'>Last</li>"
#     elif last_btn:
#         pagination_nav += "<li p='" + str(no_of_paginations) + "' class='inactive'>Last</li>"
#
#     pagination_nav = pagination_nav + "</ul></div>"
#
  #  return pagination_nav

# If
# you
# have ?students = 23, 24, 25, 26 in the
# URL, then
# you
# want
# students = request.GET.get('students')
# instead
# of
# getlist.Then
# convert
# to
# a
# list
# of
# integers
# with sd_list =[int(x) for x in students.split(',')