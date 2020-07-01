import sys
import datetime
import os
import gzip 
import subprocess
from datetime import datetime

#Check if a file exists and trigger proper actions 
#Operation can be open or write based on how you want to interact with the file
#Mode can be rename, overwrite or stop and trigger diffent actions
#When stop is selected invoke sys.exit if check fails
def checkfile(filename, operation="open", mode="rename", logger=None):
	now = datetime.now()
	
	if operation == "write":
		if os.path.isfile(filename) == True:
			if mode == "rename":
				filebase, extension = os.path.splitext(filename)
				newfilename = filebase + '.' + now.strftime("%Y%m%d_%H%M") + '.' + extension
				if logger is None:
					print(filename, "already exists")
					print("Saving to", newfilename )
				else:
					logger.warning("%s already exists", filename)
					logger.warning("Saving to: %s", newfilename)
				return newfilename
			elif mode == "overwrite":
				if logger is None:
					print(filename, "already exists and will be removed!")
				else:
					logger.warning("%s already exists and will be overwritten!", filename)
				os.remove(filename)
				return filename	
			elif mode == "stop":
				if logger is None:
					logger.critical("%s already exists! Exiting now!", filename)
					sys.exit()
				else:
					sys.exit(filename, "already exists! Exiting now!")			
		else:
			return filename
	elif operation == "open":
		if os.path.isfile(filename) == True:
			if logger is not None:
				logger.info("Open file for reading: %s", filename)
			return True
		else:
			if mode == "stop":
				if logger is None:
					sys.exit("Unable to open file", filename)
				else:
					logger.critical("Unable to open file %s", filename)
			else:
				return False

#Send email from the server
def sendEmail(address, subject, message, sender="pipeline", attachment=""):
	if type(address) is str: address = address
	if type(address) is list: address = ",".join(address)

	if os.path.isfile(message) == True:
		command='cat {0} | mailx -s "{1}" -r "{2}" {3}'.format(message, subject, sender, address) 
	else:
		command='echo "{0}" | mailx -s "{1}" -r "{2}" {3}'.format(message, subject, sender, address) 

	try:
   		subprocess.check_output(command, shell=True)
	except subprocess.CalledProcessError as grepexc:
		return grepexc.returncode
	else:
		return 0

#Return list with files having a given extension
def listFiles(mydir, myext):
	filelist = []
	for file in os.listdir(mydir):
		if file.endswith(myext):
			filelist.append(os.path.join(mydir, file))
	return filelist

#Return absolute path for a directory
def absolutePath(mydir):
	if mydir.startswith("/"):
		return mydir
	else: 
		fulldir = os.getcwd() + "/" + mydir
		return fulldir

#Create directories, eventually as subfolder of a given parent folder
#return absolute paths
def createDir(folders, parent=None):
	if parent is None:
		parent_dir = ""
	else:
		parent_dir = parent + "/"
	if type(folders) is list:
		dir_list = []
		for f in folders:
			myfolder = os.getcwd() + "/" + parent_dir + f
			dir_list.append(myfolder)
			os.makedirs(myfolder,exist_ok=True)
		return dir_list
	elif type(folders) is str:
		myfolder = os.getcwd() + "/" + parent_dir + f
		os.makedirs(myfolder,exist_ok=True)
		return myfolder
	else:
		raise TypeError('folders must be string or list')

#Return present time (HH,MM,SS) with desired separator
def now(sep=":"):
    now = datetime.now()
    current_time = now.strftime("%H{sep}%M{sep}%S".format(sep=sep))
    return current_time

#Run external command and return exitcode, stdout and stderr
def run(cmd, shell=False):
    proc = subprocess.Popen(cmd, shell=shell,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
 
    return proc.returncode, stdout.decode('utf-8'), stderr.decode('utf-8')

#Run external command and iterate over stdout
def get_stdout(cmd, shell=False):
	proc = subprocess.Popen(cmd, shell=shell,
		stdout = subprocess.PIPE,
		stderr = subprocess.DEVNULL)
    
	for line in iter(proc.stdout.readline, b''):
		yield line.decode('utf-8')

#split line into list
def tokenize(line,sep):
	line = line.rstrip('\n')
	line = line.split(sep)
	return line

#simple reader for a file plan text or gzip
def reader(file, sep, header=True, skipper=None):
	if file.endswith('gz'):
		f = gzip.open(file, 'rt')
	else:
		f = open(file, errors='replace')
	
	line = f.readline()
	if isinstance(skipper, int):
		for _ in range(skipper):
			line = f.readline()
	else:
		while skipper:
			line = f.readline()
	
	if header==True:
		colnames = tokenize(line, sep)
		line = f.readline()

	while line:
		line = tokenize(line, sep)
		if header == True:
			line = dict(zip(colnames, line))
		yield line
		line = f.readline()
	f.close()

#A modified reader for VCF, will output header to out_file
#and then yield the variants lines
#out_file is a output stream generated with open(file, "w+")
def VCFreader(file, out_file, new_anno=None):
	if file.endswith('gz'):
		f = gzip.open(file, 'rt')
	else:
		f = open(file, errors='replace')
	
	line = f.readline()
	while line.startswith("##"):
		out_file.write(line)
		line = f.readline()
	
	if new_anno: 
		for anno in new_anno: out_file.write(anno + "\n")
	
	out_file.write(line)
	line=f.readline()
	while line:
		line = tokenize(line, "\t")
		yield line
		line = f.readline()
	f.close()

#Return elements shared across a list of lists, skip non list elements
def shared_all(input_list):
	shared_set = set(input_list[0])
	for s in input_list[1:]: 
		if isinstance(s, list): 
			shared_set.intersection_update(s)
	return shared_set

#Return a flat list from a list of lists, skip non list elements
def flat(input_list):
    flat_list = [item for sublist in input_list if isinstance(sublist,list) for item in sublist]
    return flat_list