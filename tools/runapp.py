from celery import task
from tempfile import NamedTemporaryFile
from os import unlink
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

def getAppForm(application_id):
	fields = {}
	application = Application.objects.get(id=application_id)
	for option in ApplicationOptions.objects.filter(application=application):
		fields[option.name] = ModelChoiceField(queryset=ApplicationOptionsList.objects.filter(category=option).order_by('id'), empty_label=None)		
	if fields == {}:
		return False
	else:
		fields['application'] = IntegerField(initial=application.id,widget=HiddenInput())
		return type('%sForm' % str(application.name), (Form,), fields)
