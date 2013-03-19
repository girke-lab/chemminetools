import re
import os
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import redirect, render_to_response
from django.template import RequestContext
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.forms import ModelForm
from django.contrib import messages
from guest.decorators import guest_allowed, login_required
from myCompounds.views import makeSDF, uploadCompound
from compounddb.models import Compound
from tools.runapp import * 
from models import *
from simplejson import dumps

class applicationForm(ModelForm):
	class Meta:
		model = Application

class ApplicationOptionsForm(ModelForm):
	class Meta:
		model = ApplicationOptions

class ApplicationOptionsListForm(ModelForm):
	class Meta:
		model = ApplicationOptionsList

@user_passes_test(lambda u: u.is_superuser)
def manage_application(request, chooseForm):
	if request.method == 'POST': # If the form has been submitted...
		if chooseForm == 'applicationForm':
			form = applicationForm(request.POST) # A form bound to the POST data
			title = 'Add Application'
		elif chooseForm == 'ApplicationOptionsForm':
			title = 'Add Option Type'
			form = ApplicationOptionsForm(request.POST)
		else:
			title = 'Add Option Value'
			form = ApplicationOptionsListForm(request.POST)
		if form.is_valid(): # All validation rules pass
		    # Process the data in form.cleaned_data
			form.save()
			messages.success(request, 'Success: application added.')				
			return render_to_response('genericForm.html', dict(
				title=title,
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
		if chooseForm == 'applicationForm':
			form = applicationForm()
			title = 'Add Application'
		elif chooseForm == 'ApplicationOptionsForm':
			form = ApplicationOptionsForm()
			title = 'Add Option Type'
		else:
			title = 'Add Option Value'
			form = ApplicationOptionsListForm()
		return render_to_response('genericForm.html', dict(
			title=title,
			form=form,
		),
		context_instance=RequestContext(request))

@guest_allowed
def launch_job(request, category=None):
	if request.is_ajax():
		# for ajax requests, return HTML form for each app
		currentApp = request.GET['currentApp']
		try:
			application = Application.objects.get(id=currentApp)
			AppFormSet = getAppForm(application.id, request.user)	
			form = AppFormSet()
			form = str(form)		
			response = dict(form=form, desc=application.description)
		except:
			response = dict(form="ERROR")
		return HttpResponse(dumps(response),'text/json')
	if request.method == 'POST':
		appForm = getAppForm(request.POST['application'], request.user)
		form = appForm(request.POST, request.FILES)
		if form.is_valid():
			try:
				appid = int(form.cleaned_data['application'])
				application = Application.objects.get(id=str(appid))
			except Application.DoesNotExist:
				messages.error(request, 'Application does not exist')
				return redirect(launch_job, category=category)
		else:
			messages.error(request, str(form.errors))
			return redirect(launch_job, category=category)
		commandOptions = u''
		optionsList = u''
		for question in form.cleaned_data.keys():
			if question != 'application' and question != 'File Upload':
				questionObject = ApplicationOptions.objects.get(application=application, name=question) 
				try:
					job = form.cleaned_data[question]
					option = job.output
					optionName = str(job) 
				except:
					try:
						answerObject = form.cleaned_data[question]
						optionName = answerObject.name
						option = answerObject.realName
					except:
						option = 'None'
						optionName = 'None'
				commandOptions = commandOptions + " --" + questionObject.realName + "=" + option 
				optionsList = optionsList + questionObject.name + ": " + optionName + ", "
		optionsList = re.sub(",\s$", "", optionsList, count=0)

		# setup input
		if application.input_type == 'chemical/x-mdl-sdfile':
			input = makeSDF(request.user)
		elif application.input_type == 'upload':
			input = request.FILES['File Upload'].read()	
			if not isinstance(input, unicode):
				input = unicode(input, 'utf-8')
			input = input.encode('ascii', 'ignore')
		else:
			input = ''
		newJob = Job(
			user=request.user,
			application=application,
			options=optionsList,
			input='myCompounds sdf',
			output='',
			task_id='',
		)
		newJob.save()
		result = launch.delay(application.script, commandOptions, input, newJob.id)
		newJob.task_id = result.id
		newJob.save()
		messages.success(request, 'Success: job launched.')
		return redirect(view_job, job_id=newJob.id, resource='', filename='')
	else:
		if category:
			try:
				category = ApplicationCategories.objects.get(name=category)
				compoundCount = Compound.objects.filter(user=request.user).count()
				if category.name == 'Clustering' and compoundCount < 3:
					messages.info(request, 'Notice: you must have at least 3 compounds to perform clustering. Please use this form to add more compounds and then click "Cluster" again.')
					return redirect(uploadCompound)
				if category.name == 'Properties' and compoundCount < 1:
					messages.info(request, 'Notice: you must have at least one compound to compute properties. Please use this form to add compounds and then click "Properties" again.')
					return redirect(uploadCompound)
				title = "Launch " + category.name + " Job"
				apps = Application.objects.filter(category=category)
			except:
				raise Http404	
		else:
			title = "Launch Job"
			apps = Application.objects.filter()
		fields = {}
		fields['application'] = ModelChoiceField(queryset=apps, empty_label="")
		form = type('%sForm' % 'choose application', (Form,), fields)
		return render_to_response('submitForm.html', dict(
			title=title,
			form = form,
		),
		context_instance=RequestContext(request)) 

@guest_allowed
def view_job(request, job_id, resource, filename):
	try:
		job = updateJob(request.user, job_id)
	except Job.DoesNotExist:
		raise Http404
	if resource:
		if resource == 'delete':
			deleteJob(request.user, job.id)
			return HttpResponse("deleted", mimetype='text/plain')
		if resource == 'download':
			f = open(job.output, 'r')
			result = f.read()
			f.close()
			return HttpResponse(result, mimetype=job.application.output_type)
	if request.is_ajax():
		if job.status == Job.RUNNING:
			response = dict(status="RUNNING")
		else:
			response = dict(status="DONE")
		return HttpResponse(dumps(response),'text/json')
	if job.status == Job.FINISHED:
		finalResult = job.output 
		finalResult = re.sub(".*/", "", finalResult, count=0)
		finalResult = '/working/' + finalResult
		# select correct viewer here based on output type
		if(job.application.output_type == 'application/json.canvasxpress'):	
			f = open(job.output, 'r')
			plotJSON = f.read()
			f.close()
			return render_to_response('view_job.html', dict(
				title = str(job.application) + " Results",
				result = finalResult,
				job = job,
				plotJSON = plotJSON,
			),
			context_instance=RequestContext(request))
		elif(job.application.output_type == 'text/properties.table'):
			f = open(job.output, 'r')
			textfile = f.read().rstrip()
			f.close()
			csv = []
			for line in textfile.split('\n'):
				csv.append(line.rstrip().split(','))
			return render_to_response('view_csv.html', dict(
				title = str(job.application) + " Results",
				result = finalResult,
				job = job,
				csv = csv,
			),
			context_instance=RequestContext(request))
		else:
			# if mimetype is unknown, just send the file to the user
			return redirect(view_job, job_id=job.id, resource='download', filename='output')	
	elif job.status == Job.RUNNING:
		return render_to_response('wait.html', dict(
			title = 'Job Running <img src="/static/images/spinner.gif"/>',
		),
		context_instance=RequestContext(request))
	elif job.status == Job.FAILED:
		return render_to_response('view_job.html', dict(
			title = "Error: " + str(job.application) + " Job Failed",
			job = job,
		),
		context_instance=RequestContext(request))

@guest_allowed
def list_jobs(request):
	matches = getJobList(request.user)
	return render_to_response('list_jobs.html', dict(matches=matches,), context_instance=RequestContext(request))
