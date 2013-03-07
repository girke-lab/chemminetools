import os
from celery import task
from tempfile import NamedTemporaryFile
from django.conf import settings
import subprocess
import random
from django.forms import Form, ModelChoiceField, IntegerField, HiddenInput 
from models import *
outputPath = settings.TOOLS_RESULTS
projectDir = settings.PROJECT_DIR

@task()
def launch(appname, commandOptions, input):
	outputFile = NamedTemporaryFile(dir=outputPath, delete=True)
	outputFileName = outputFile.name
	outputFile.close()
	command = projectDir + '/tools/tool_scripts/' + appname + ' --outfile=' + outputFileName + " " + commandOptions
	print("Running: " + command + "\n")
	runningTask = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
	runningTask.stdin.write(input)
	runningTask.stdin.close()
	runningTask.wait()
	if runningTask.returncode != 0:
		return False 
	else:
		return outputFileName 

def getAppForm(application_id, user):
	fields = {}
	application = Application.objects.get(id=application_id)
	for option in ApplicationOptions.objects.filter(application=application).order_by('id'):
		try:
			value = ApplicationOptionsList.objects.get(name='fileInput', category=option)
			getJobList(user)
			fields[option.name] = ModelChoiceField(queryset=Job.objects.filter(application__output_type=value.realName, user=user, status=Job.FINISHED).order_by('-id'), empty_label="None", required=False)	
		except:
			fields[option.name] = ModelChoiceField(queryset=ApplicationOptionsList.objects.filter(category=option).order_by('id'), empty_label=None)		
	fields['application'] = IntegerField(initial=application.id,widget=HiddenInput())
	return type('%sForm' % str(application.name), (Form,), fields)

def getJobList(user):
	# get all running jobs and update those that need updating
	try:
		runningJobs = Job.objects.filter(user=user, status=Job.RUNNING)
		for job in runningJobs:
			updateJob(user,job.id)
	except:
		return False
	return Job.objects.filter(user=user).order_by('-start_time')


def updateJob(user, job_id):
	# checks if a job is done, updates it's status, and then returns it
	try:
		job = Job.objects.get(id=job_id, user=user)
	except:
		return False
	if job.status == Job.RUNNING:
		result = launch.AsyncResult(job.task_id)
		if result.ready():
			if result.result:
				job.status = Job.FINISHED
			else:
				job.status = Job.FAILED 
			job.output = result.result
			result.forget()
			job.save()
	return job

def deleteJob(user, job_id):
	job = updateJob(user, job_id)
	job.delete()
	return True

def deleteApp(name):
	try:
		app = Application.objects.get(name=name)
	except:
		return False
	app.delete()
	return True
