class ForceHttpHeadersMiddleware:
    """A middleware class for overriding HTTP headers that may have been
    placed by other middleware. To use, add an attribute named
    'cmt_force_http_headers' to the response object with a dict object:

    hr = HttpResponse(...)
    d = dict()
    d['Cache-Control'] = 'max-age=9000, immutable'
    d['Vary'] = None
    hr.cmt_force_http_headers = d
    return hr
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)

        try:
            for header, value in response.cmt_force_http_headers.items():
                if value is None:
                    del response[header]
                else:
                    response[header] = value
        except AttributeError:
            pass

        return response
