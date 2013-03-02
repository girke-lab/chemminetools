from celery import task
from tempfile import NamedTemporaryFile
from os import unlink
from django.conf import settings
import subprocess
import random
outputPath = settings.TOOLS_RESULTS
projectDir = settings.PROJECT_DIR

@task()
def launch(appname, commandOptions, input):
	outputFile = NamedTemporaryFile(dir=outputPath, delete=True)
	outputFileName = outputFile.name
	outputFile.close()
	command = projectDir + '/tools/tool_scripts/' + appname + ' --outfile=' + outputFileName
	runningTask = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
	runningTask.stdin.write(input)
	runningTask.stdin.close()
	runningTask.wait()
	if runningTask.returncode != 0:
		return False 
	else:
		return outputFileName 
