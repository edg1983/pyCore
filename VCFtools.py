import pyCore.constants as constants
import re
from collections import OrderedDict
import pandas as pd
import gzip

'''
parse data from VCF sample columns
samplesID is a list containing samples in the same order as in VCF header
getSample is ALL to return data from all sample or you can specify a single sample
tag is ALLTAGS to process all information or you can specificy a single tag (like GQ)
miss is the value returned if the tag is not found
'''
def parseFORMAT(vcfline, samplesID, getSample="ALL", tag="ALLTAGS", miss="NA"):
	values = {}
	fields = vcfline.split("\t")
	vcfFormat = fields[8].split(":")
	for i in range(0, len(vcfFormat)): values[vcfFormat[i]] = {}
	
	for i in range(9, len(fields)):
		myID = samplesID[i-9]
		mysample = fields[i].split(":")
		for n in range(0, len(mysample)):
			if mysample[n] == ".": mysample[n] = miss
			values[vcfFormat[n]][myID] = mysample[n]
	
	if tag == "ALLTAGS" and getSample == "ALL": #return nested dictionary [tag][sampleID]
		return values 
	elif tag == "ALLTAGS" and getSample != "ALL": #return dictionary [tag]
		output = {}
		for key in values.keys():
			output[key] = values[key][getSample]
		return output 
	elif tag != "ALLTAGS" and getSample == "ALL": #return list of tag values across samples
		output = []
		for key, value in values[tag].items():
			output.append(value)
		return output 
	elif tag != "ALLTAGS" and getSample != "ALL": # return single value
		return values[tag][getSample] 

def parseINFO(vcfline, tag="ALLTAGS", miss="NA"):
	infos = {}
	fields = vcfline.split("\t")
	info_tags = fields[7].split(";")
	for value in info_tags:
		p = re.compile('(.+)=(.+)')
		m = p.match(value)
		if m:
			infos[m.group(1)] = m.group(2)
		else:
			infos[value] = True
	
	if tag != "ALLTAGS": 
		if tag in infos.keys():
			myinfo = infos[tag]
		else:
			myinfo = miss
	else:
		myinfo = infos
	
	return myinfo

'''
A useful class to deal with gene-level consequence annotations
Right now works just for snpEff format
worst property contains the worst consequence
'''
class Consequence:
	def __init__(self, annfield, tool):
		self.annotation = []
		self.worst = {}
		bestrank = 99

		if tool is "snpEff":
			fieldnames = ['Allele','Annotation','Annotation_Impact','Gene_Name','Gene_ID', \
				'Feature_Type','Feature_ID','Transcript_BioType','Rank','HGVS.c', 'HGVS.p', \
				'cDNA.pos','cDNA.length','CDS.pos','AA.pos','Distance','ERRORS']
			sep1 = ","
			sep2 = "|"

		records = annfield.split(sep1)
		for r in records:
			annotations = {}
			annotations['ANN'] = r
			fields = r.split(sep2)
			for n,val in enumerate(fields):
				annotations[fieldnames[n]] = val
			
			terms = annotations['Annotation'].split('&')
			terms.sort(key=lambda x: constants.CONSEQUENCE_RANKING.index(x))
			annotations['Annotation'] = terms[0]
			if constants.CONSEQUENCE_RANKING.index(annotations['Annotation']) < bestrank:
				myworst = annotations
				bestrank = constants.CONSEQUENCE_RANKING.index(annotations['Annotation'])
			self.annotation.append(annotations)
		self.worst = myworst 

	#Extract the selected field from all the annotations and return a list	
	def getField(self, field):
		output = []
		for a in self.annotation:
			output.append(a[field])
		output = list(set(output))
		return output

	#Sort the annotations from worst to low impact
	def sortByAnnotation(self):
		self.annotation.sort(key=lambda x: constants.CONSEQUENCE_RANKING.index(x['Annotation']))
	
	#Remove any annotation with consequence below the specified threshold
	def cleanByImpact(self, impact):
		threshold = constants.CONSEQUENCE_RANKING.index(impact)
		self.annotation = [a for a in self.annotation if not constants.CONSEQUENCE_RANKING.index(a['Annotation']) > threshold]

	#Evaluate if there is at least one annotation with consequence more severe than threshold
	#Return True/False	
	def evalConsequence(self, impact):
		threshold = constants.CONSEQUENCE_RANKING.index(impact)
		result = False
		for a in self.annotation:
			if constants.CONSEQUENCE_RANKING.index(a['Annotation']) <= threshold:
				result = True
				break
		return result

	#Evaluate if there is at least one annotation with IMPACT more severe than threshold
	#Return True/False
	def evalImpact(self, impact):
		threshold = constants.IMPACT_RANKING.index(impact)
		result = False
		for a in self.annotation:
			if constants.IMPACT_RANKING.index(a['Annotation_Impact']) <= threshold:
				result = True
				break
		return result

	#create a worstByGene list that contains only the worst consequence per gene
	def getWorstByGene(self):
		self.worstByGene = []
		bygene = {}
		for a in self.annotation:
			if bygene[a['Gene_ID']] is None: bygene[a['Gene_ID']] = []
			bygene[a['Gene_ID']].append(a)
		for g in bygene.keys():
			bygene[g].sort(key=lambda x: constants.CONSEQUENCE_RANKING.index(x['Annotation']))
			self.worstByGene.append(bygene[g][0])

	#Return a string formatted as snpEff ANN field. 
	#Output can be set to all (return all annotations) or worst (for the worst only)
	def printANN(self, output="all"):
		if output == "all":
			ANNfield = []
			for a in self.annotation: 
				ANNfield.append(a['ANN'])
			return ",".join(ANNfield)
		elif output == "pergene":
			self.getWorstByGene()
			for a in self.worstByGene: 
				ANNfield.append(a['ANN'])
			return ",".join(ANNfield)
		elif output == "worst":
			return self.worst['ANN']

'''
Apply filters based on INFO column. 
Use parseINFO to parse INFO tags and values
conditions is a list of conditions to evaluate (like DP>20)
logic set the logic for multiple conditions evaluation (AND / OR)
filtertag set the name of the filter to return in case the variant fails the filtering
'''
def evaluateINFO(line, conditions, logic="AND"):

	if logic == "AND": 
		threshold = len(conditions)
	elif logic == "OR":
		threshold = 1

	failed = 0
	for c in conditions:
		m = re.match("(.+)([<>=!]+)(.+)", c)
		value = parseINFO(line, m[1], False)
		if value != False:
			exp = value + m[2] + m[3]
			if eval(exp):
				passed_conditions += 1

	if passed_conditions >= threshold:
		return True
	else:
		return False

'''
Apply filters based on sample information. 
Use parseFORMAT to parse information from the single sample columns
conditions is a list of conditions to evaluate (like GQ>20)
samplesID is a list containing the sampleID from the header line
applyTo decide if the filter is applied to all sample or only a specific one
cond_logic set the logic for condition evaluation (AND / OR)
sample_logic set the logico for multiple sample evaluation (ALL: filter if all samples fail the conditions, ANY: filter if one sample fails)
filtertag set the name of the filter to return in case the variant fails the filtering
'''
def evaluateSAMPLE(line, conditions, samplesID, cond_logic="AND"):

	if cond_logic == "AND": 
		thr_cond = len(conditions)
	elif cond_logic == "OR":
		thr_cond = 1

	sample_pass = {}
	for s in samplesID:
		sample_pass[s] = 0
		pass_cond = 0
		for c in conditions:
			m = re.match(r"([^<>=!]+)([<>=!]+)(.+)", c)
			value = parseFORMAT(line,samplesID,s,m[1],False)
			if value != False:
				exp = value + m[2] + m[3]
				if eval(exp):
					pass_cond += 1
		if pass_cond >= thr_cond:
			sample_pass[s] = 1
	
	return sample_pass

'''
Class to manage file headers
Can process generic headers, single line headers as well as VCF headers
For VCF header property is dictionary containing the different header sections
You have also specific method to add FILTER, INFO or FORMAT lines
header['DICT'] is a nested dictionary of sequence name and lengths
You have also a specific method changeDICT to change/append/update the header
'''
class Header():
	def __init__(self, file, preset="VCF", headchar="#", sep="\t"):
		line = file.readline()
		self.header = []
		self.__type = preset

		#Generic header: collect all lines starting with headchar
		#header is a list of lines
		if preset == "generic":
			while line.startswith(headchar):
				line = line.rstrip("\n")
				self.header.append(line)
				line = file.readline()
			self.firstline = line

		#Parse a single line header and return list of column names
		#header is a list of col names
		if preset == "singleline":
			line.rstrip("\n")
			self.header = line.split(sep)
			line = file.readline()
			self.firstline = line

		#Parse a VCF header
		#header is a dictionary for the different header sections
		if preset == "VCF":
			self.samples = []
			self.header = {}
			self.header['TOP'] = []
			self.header['DICT'] = {}
			self.header['INFO'] = []
			self.header['FILTER'] = []
			self.header['FORMAT'] = []
			self.header['OTHER'] = []
			while line.startswith(headchar):
				line.rstrip("\n")
				if line.startswith(('##fileformat', '##fileDate', '##source', '##reference')):
					self.header['TOP'].append(line)
				if line.startswith('##contig'):
					m = re.search(r"ID=([^,]+),length=([0-9]+)", line)
					seq = m.group(1), 
					if m.group(2) is not None: length = m.group(2) 
					self.header['DICT'][seq] = length
				if line.startswith('##INFO'):
					self.header['INFO'].append(line)
				if line.startswith('##FILTER'):
					self.header['FILTER'].append(line)
				if line.startswith('##FORMAT'):
					self.header['FORMAT'].append(line)
				if line.startswith('#CHROM'):
					col_names = line.split("\t")
					self.headline = line
					self.samples = col_names[9:]
				line = file.readline()
			self.firstline = line

	#Substitute header with a new one. Can accept header class, a list or a string
	def change(self, newheader):
		if isinstance(newheader, Header):
			self.header = newheader.header
		elif isinstance(newheader, list):
			for h in newheader: 
				self.header.clear()
				self.header.append(h)
		elif isinstance(newheader, str):
			self.header.clear()
			self.header.append(newheader)
	
	#Add a new line to an exisisting header
	def appendLine(self, newline):
		if isinstance(newline, list):
			for h in newline: self.header.append(h)
		elif isinstance(newline, str):
			self.header.append(newline)
	
	#Create a new INFO line for a VCF header
	def appendINFO(self, tag, datatype, description="new annotation", number=1):
		newline = '##INFO=<ID={newtag},Number={newnumber},Type={newtype},Description="{newdescription}">'.format(newtag=tag,newnumber=number,newtype=datatype,newdescription=description)
		self.header['INFO'].append(newline)

	#Create a new FILTER line for a VCF header
	def appendFILTER(self, filter_name, description="new filter"):
		newline = '##FILTER=<ID={newfilter},Description="{newdescription}">'.format(newfilter=filter_name,newdescription=description)
		self.header['FILTER'].append(newline)

	#Create a new FORMAT line for a VCF header
	def appendFORMAT(self, tag, datatype, description="new annotation", number=1):
		newline = '##FORMAT=<ID={newtag},Number={newnumber},Type={newtype},Description="{newdescription}">'.format(newtag=tag,newnumber=number,newtype=datatype,newdescription=description)
		self.header['FORMAT'].append(newline)

	#Manipulate sequence dictionary, dict as input
	def changeDICT(self, newdict):
		for seq, length in newdict.items():
			self.header['DICT'][seq] = length

	#Remove a sequence / list of sequences from the seq dictionary 
	def removeDICTseq(self, seqname):
		if isinstance(seqname, list):
			for s in seqname:
				del self.header['DICT'][s]
		elif isinstance(seqname, str):
			del self.header['DICT'][seqname]

	#Return a string of text containing the whole header for printing/writing
	def printHeader(self):
		if self.__type == "generic":
			fullheader = "\n".join(self.header)
		elif self.__type == "singleline":
			fullheader = "\t".join(self.header)
		elif self.__type == "VCF":
			fullheader = "\n".join(self.header['TOP']) + "\n" + \
				"\n".join(self.header['DICT']) + "\n" + \
				"\n".join(self.header['INFO']) + "\n" + \
				"\n".join(self.header['FILTER']) + "\n" + \
				"\n".join(self.header['FORMAT']) + "\n" + \
				self.headline
		return fullheader

class INFO():
	def __init__(self, info_field):
		exp = re.compile('(.+)=(.+)')
		self.infos = OrderedDict()
		info_tags = info_field.split(";")
		for value in info_tags:
			m = exp.match(value)
			if m:
				self.infos[m.group(1)] = m.group(2)
			else:
				self.infos[value] = True
	
	def selectINFO(self, info_keys):
		if isinstance(info_keys, list):
			self.infos = { key: self.infos[key] for key in info_keys }
		else:
			self.infos = { info_keys : self.infos[info_keys]}
	
	def dropINFO(self, info_keys):
		if isinstance(info_keys, list):
			for key in info_keys: self.infos.pop(key)
		else:
			self.infos.pop(key)
    
	def addINFO(self, info_key, value):
		self.infos[info_key] = value
    
	def updateINFO(self, info_key, value):
		if info_key in self.infos.keys():
			self.infos[info_key] = value
		else:
			self.addINFO(info_key, value)
	
	def getINFO(self,format):
		if format == "dict":
			return self.infos
		elif format == "string":
			output = [] 
			for key, value in self.infos.items():
				output.append(key + "=" + str(value))
		return ';'.join(output)

class Variant():
	def __init__(self,line,samples):
		if len(line) < 8:
			raise Exception("Malformed VCF")
		self.CHROM = line[0]
		self.POS = line[1]
		self.ID = line[2]
		self.REF = line[3]
		if re.match(',', line[4]):
			self.TYPE = "MNP"
			self.ALT = line[4].split(',') 
		else:
			self.ALT = line[4]
			if len(line[3]) == len(line[4]):
				self.TYPE = "SNV"
			else:
				self.TYPE = "INSDEL"
			
		self.QUAL = line[5]
		self.FILTER = line[6]
		self.INFO = INFO(line[7])
				
		try:
			self.ANN = Consequence(self.INFO.infos['ANN'], "snpEff")
		except:
			self.ANN = 0

		try:
			self.FORMAT = line[8].split(':')
		except:
			self.FORMAT = "GT"
				
		if len(line) >= 10:
			self.SAMPLES = OrderedDict()
			i = -1
			for sample_field in line[9:]:
				i += 1
				sample_field = sample_field.split(':')
				self.SAMPLES[samples[i]] = OrderedDict(zip(self.FORMAT,sample_field))
	
	def updatePosition(self, location):
		self.CHROM = location[0]
		self.POS = location[1]

	def writeVariant(self):
		if isinstance(self.ALT, list):
			ALT = ','.join(self.ALT)
		else:
			ALT = self.ALT

		SAMPLES = [':'.join(sample.values()) for sample in self.SAMPLES.values()]
		
		outline = [
			self.CHROM,
			self.POS,
			self.ID,
			self.REF,
			ALT,
			self.FILTER,
			self.QUAL,
			self.INFO.getINFO("string"),
			':'.join(self.FORMAT),
			'\t'.join(SAMPLES)
		]

		return '\t'.join(outline)

class VCF():
	def __init__(self, VCFfile, interval=None):
		if VCFfile.endswith('.gz'): 
			VCFfile = gzip.open(VCFfile, 'rt')
		elif VCFfile.endswith('.vcf'): 
			VCFfile = open(VCFfile)
		else:
			raise Exception("VCF or VCF.GZ extension expected")
		
		chrom = 'ALL'
		start = 0
		end = 99999999999
		
		if interval is not None:
			m = re.match(r'([^:]+):([0-9]+)-([0-9]+)', interval)
			if m is not None:
				chrom =  m.group(1)
				start = int(m.group(2))
				end = int(m.group(3))
			else:
				chrom = interval

		with VCFfile as f:
			self.header = Header(f)
			line = self.header.firstline
			self.variants = []
			self.stats = {
				'tot_variants' : [0, 0],
				'SNV' : [0, 0],
				'INSDEL' : [0, 0],
				'MNP' : [0, 0]}
			n = 0
			while line:
				n += 1
				#print("Line", n, end='\r')
				line = line.rstrip('\n')
				line = line.split('\t')
				
				if chrom != 'ALL':
					if line[0] != chrom: 
						line = f.readline()
						continue
				if int(line[1]) > end: break
				if not start < int(line[1]) < end: 
					line = f.readline()
					continue
				
				var = Variant(line,self.header.samples)
				self.variants.append(var)
				for var_type in ['SNV', 'INSDEL', 'MNV']:
					if var.TYPE == var_type: 
						self.stats[var_type][0] += 1
						self.stats['tot_variants'][0] += 1
						if var.FILTER == "PASS": 
							self.stats[var_type][1] += 1
							self.stats['tot_variants'][1] += 1
				line = f.readline()

	def readVariants(self):
		for var in self.variants:
			yield var
	
	def removeSamples(self, samples):
		for var in self.variants:
			var.SAMPLES_DF.drop(var.SAMPLES_DF.index[[samples]])

	def selectSamples(self, samples):
		for var in self.variants:
			var.SAMPLES_DF = var.SAMPLES_DF.loc[samples]
