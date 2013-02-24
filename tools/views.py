from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.forms import ModelForm
from django.contrib import messages
from guest.decorators import guest_allowed, login_required
from myCompounds.views import makeSDF
from runapp import launch
from models import *

class applicationForm(ModelForm):
	class Meta:
		model = Application

class jobForm(ModelForm):
	class Meta:
		model = Job
		fields = ('application',)

@user_passes_test(lambda u: u.is_superuser)
def manage_application(request):
	if request.method == 'POST': # If the form has been submitted...
		form = applicationForm(request.POST) # A form bound to the POST data
		if form.is_valid(): # All validation rules pass
		    # Process the data in form.cleaned_data
			form.save()
			messages.success(request, 'Success: application added.')				
			return render_to_response('genericForm.html', dict(
				title='Add Application',
				form=form,
			),
			context_instance=RequestContext(request)) 
		else:
			messages.error(request, 'error: invalid form data.')
			return render_to_response('genericForm.html', dict(
				title='Add Application',
				form=form,
			),
			context_instance=RequestContext(request)) 
	else:
		form = applicationForm() # An unbound form
		return render_to_response('genericForm.html', dict(
			title='Add Application',
			form=form,
		),
		context_instance=RequestContext(request))

@guest_allowed
def launch_job(request):
	username = request.user.username
	if request.method == 'POST':
		form = jobForm(request.POST)
		if form.is_valid():
			try:
				application = Application.objects.get(id__iexact=form['application'].value())
			except Application.DoesNotExist:
				raise Http404
		else:
			raise Http404

		sdf = makeSDF(username)
		result = launch.delay(application.script, "", sdf)
		newJob = Job(
			username=username,
			application=application,
			options='',
			input='myCompounds sdf',
			output='',
			task_id=result.id,
		)
		newJob.save()
		messages.success(request, 'Success: job launched.')
		return render_to_response('genericForm.html', context_instance=RequestContext(request))
	else:
		form = jobForm()
		return render_to_response('genericForm.html', dict(
			title='Start Job',
			form = form,
		),
		context_instance=RequestContext(request)) 

