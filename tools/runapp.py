from celery import task
from django.conf import settings
import subprocess
import random
outputPath = settings.TOOLS_RESULTS
projectDir = settings.PROJECT_DIR

@task()
def add(x, y):
    return x + y

@task()
def launch(appname, commandOptions, input):
	outputFile = outputPath + '/random' + str(random.randint(10000000, 90000000))
	command = projectDir + '/tools/tool_scripts/' + appname + ' > ' + outputFile
	runningTask = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
	runningTask.stdin.write(input)
	runningTask.stdin.close()
	runningTask.wait()
	if runningTask.returncode != 0:
		return False 
	else:
		return outputFile
