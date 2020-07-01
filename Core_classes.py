import sys, os
sys.path.append('/well/gel/HICF2/software/BRC_tools/pyCore/')
import re
import pandas as pd
import json, requests, obonet, networkx
import subprocess
from datetime import datetime
from Core import run, flat, shared_all

#Default file locations on rescomp
GADO_JAR="/well/gel/HICF2/software/GadoCommandline-1.0.1/GADO.jar"
GADO_DATA="/well/gel/HICF2/ref/GADO_resources"
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1"
OBO_FILE = "/well/gel/HICF2/ref/HPO/hp_20191120.obo"
HPO2GENE = "/well/gel/HICF2/ref/HPO/phenotype_to_genes_20191120.txt"
DISEASE2HPO = "/well/gel/HICF2/ref/HPO/phenotype_annotation_20191120.tab"
CLINVAR_PATH_GENES="/well/gel/HICF2/ref/ClinVar/pathogenic_genes_conditions_20190704.std.txt"

#Interface to HPO ontology. Allow extracting genes for given HPOs or disease terms
class HPO():
    def __init__(self, obo_file, hpo2gene, disease2hpo):
        self._obofile = obo_file
        self._disease_df = pd.read_csv(disease2hpo, sep="\t", usecols=[0,1,2,4], names=['source','disease_id','disease','HPO_id'],comment="#")
        self._gene_df = pd.read_csv(hpo2gene, sep="\t", usecols=[0,3], names=['HPO_id','gene'],comment="#")
        self._HPO2gene = self._gene_df.groupby(by='HPO_id')['gene'].apply(list).reset_index(name='genes')
        self._disease2genes = self._disease_df.merge(self._HPO2gene, on='HPO_id')
        self._disease2genes = self._disease2genes.groupby(by=['source','disease_id'])['genes'].agg(genes=pd.NamedAgg(column='genes', aggfunc='sum'))

        self._disease_df.set_index(['source','disease_id'], inplace=True)
        self._gene_df.set_index('gene', inplace=True)
        self._HPO2gene.set_index('HPO_id', inplace=True)
        
        self._ontology, self._obsoletes = obonet.read_obo(obo_file)
        self.id_to_name = {id_: data.get('name') for id_, data in self._ontology.nodes(data=True)}
        self.name_to_id = {data['name']: id_ for id_, data in self._ontology.nodes(data=True) if 'name' in data}
        self.n_terms = len(self._ontology)

    #Get genes from HPO ids
    #If shared=True return only terms shared between provided HPOs    
    def genes_by_hpo(self, hpo_ids, shared=False):
        genes = list(self._HPO2gene.loc[hpo_ids,'genes'])
        if shared:
            out_list = shared_all(genes)
        else:
            out_list = set(flat(genes))
        return out_list
    
    #Get genes from disease IDs (default OMIM id)
    #Other sources can be specified with source option (DECIPHER, ecc)
    #If shared=True return only terms shared between provided diseases 
    def genes_by_diseaseID(self, id_list, shared=False, source='OMIM'):
        id_list = [int(x) for x in id_list]
        genes = list(self._disease2genes.loc[(source, id_list),'genes'])
        if shared:
            out_list = shared_all(genes)
        else:
            out_list = set(flat(genes))
        return out_list

    #Get subterms for a list of HPO ids
    def get_subterms(self, hpo_ids):
        output_ids = []
        output_names = []
        for hpo in hpo_ids:
            output_names.extend([self.id_to_name[subterm] for subterm in networkx.ancestors(self._ontology, hpo)])
            output_ids.extend(networkx.ancestors(self._ontology, hpo))
        return output_ids, output_names
    
    #Get disease ids from a list of keywords.
    #Returns list if single id is specified as source (like OMIM)
    #Returns dict of source, ids if ALL is specified as source
    def get_disease_ids(self, terms_list, exact_match=False, source="OMIM"):
        if source == "ALL":
            sources = [x[0] for x in self._disease_df.index]
            disease_ids = dict.fromkeys(sources,set())
            if exact_match:
                for term in terms_list:
                    idx = self._disease_df[self._disease_df.disease == term].index
                    for key, value in idx:
                        disease_ids[key].add(value)
            else:
                for term in terms_list:
                    selection = pd.DataFrame(self._disease_df['disease'].str.contains(term, case=False))
                    idx = selection.index[selection.disease == True]
                    for key, value in idx:
                        disease_ids[key].add(value)                
        else:
            disease_ids = set()
            if exact_match:
                for term in terms_list:
                    idx = self._disease_df[self._disease_df.disease == term].index
                    disease_ids.update([x[1] for x in idx if x[0]==source])                
            else:
                for term in terms_list:
                    selection = pd.DataFrame(self._disease_df['disease'].str.contains(term, case=False))
                    idx = selection.index[selection.disease == True]
                    disease_ids.update([x[1] for x in idx if x[0]==source])

        return disease_ids
    
    #For a list of HPO ids update obsolete ones according to ontology obo
    def update_obsolete(self, hpo_ids):
        new_ids = []
        for hpo_id in hpo_ids:
            new_ids.append(self._obsoletes.get(hpo_id, hpo_id))
        return new_ids

#Interface to PanelApp API. Allow extraction of genes using panel IDs or disease terms
class PanelApp():
    def __init__(self, url):
        self._url = url
        panels_list = list()
        uri = PANELAPP + "/panels"
        while ( uri is not None):
            print(uri)
            get_request = requests.get(uri)
            if get_request.status_code != 200: raise Exception("Unable to get panels from", uri) 
            ppage = get_request.json()
            panels_list = panels_list + ppage['results']
            uri=ppage["next"]
        self._panels = pd.DataFrame.from_dict(panels_list)
        self._panels['n_genes'] = [x['number_of_genes'] for x in list(self._panels['stats'])]
        self._panels = self._panels[['id','name','disease_group','disease_sub_group','version','version_created','n_genes','relevant_disorders']]
        self._panels.set_index('id', inplace=True)
        now = datetime.now()
        self.time = now.strftime("%Y{sep}%m{sep}%d".format(sep=""))

    def listPanels(self, panel_id=False, name=False, disease=False):
        panels = []
        
        if not name and not disease and not panel_id:
            name='.*'
        
        if disease:
            r = re.compile(disease, re.IGNORECASE)
            mask = []
            for line in list(self._panels.relevant_disorders):
                mask.append(any(r.search(x) for x in line))
            selected_panels_disease = self._panels[mask]
            panels.append(selected_panels_disease)
            
        if name:
            selected_panels_name = self._panels[self._panels.name.str.contains(name,case=False,regex=True)]
            panels.append(selected_panels_name)
            
        if panel_id:
            selected_panels_id = pd.DataFrame(self._panels.loc[panel_id,]).transpose()
            panels.append(selected_panels_id)
        
        return pd.concat(panels)
              
    def getPanelId(self, name='.*', disease='.*'):
        selected_panels = self.listPanels(name, disease)
        return list(selected_panels.index)
    
    def getGenes(self, pid=False, name=False, disease=False, level=3, out_format="df", build="GRCh38"):
        genes = pd.DataFrame(columns = ['entity_name','confidence_level','panel_name','penetrance','mode_of_inheritance','GRCh37','GRCh38'])
        panels = self.listPanels(pid,name,disease)
        ids = list(panels.index)
        for panel_id in ids:
            get_request = requests.get(self._url + "/panels/" + str(panel_id))
            if get_request.status_code != 200: 
                print("Unable to get genes for panel ", panel_id)
                continue
            
            if len(get_request.json()['genes']) == 0:
                continue
            else:
                panel_genes = pd.DataFrame.from_dict(get_request.json()['genes'])
                GRCh37_coords = list()
                GRCh38_coords = list()
                for x in panel_genes['gene_data']:
                    try:
                        GRCh37_coords.append(x['ensembl_genes']['GRch37']['82']['location'])
                    except:
                        GRCh37_coords.append("0:0-0")
                    try:
                        GRCh38_coords.append(x['ensembl_genes']['GRch38']['90']['location'])
                    except:
                        GRCh38_coords.append("0:0-0")
                panel_genes['GRCh37'] = GRCh37_coords
                panel_genes['GRCh38'] = GRCh38_coords
                panel_genes['panel_name'] = panels.name[panel_id]
                panel_genes['confidence_level'] = pd.to_numeric(panel_genes['confidence_level'])
                panel_genes = panel_genes[panel_genes.confidence_level >= level]
                panel_genes = panel_genes[['entity_name','confidence_level','panel_name','penetrance','mode_of_inheritance','GRCh37','GRCh38']]
                genes = pd.concat([genes, panel_genes])
            
                if out_format == "bed":
                    chrs_order = [str(x) for x in list(range(0,23)) + ['X','Y','M']]
                    coordinates = genes[build].str.extractall(r'([0-9XYM]+):(\d+)-(\d+)')
                    genes['chrom'] = pd.Categorical(coordinates[0].values, chrs_order)
                    genes['start'] = pd.to_numeric(coordinates[1].values)
                    genes['stop'] = pd.to_numeric(coordinates[2].values)
                    genes = genes[['chrom', 'start', 'stop','entity_name','panel_name']]
                    genes.sort_values(by=['chrom','start'], inplace=True)
            
                if out_format == "detailed_bed":
                    chrs_order = [str(x) for x in list(range(0,23)) + ['X','Y','M','MT']]
                    coordinates = genes[build].str.extractall(r'([0-9XYMT]+):(\d+)-(\d+)')
                    genes['chrom'] = pd.Categorical(coordinates[0].values, chrs_order)
                    genes['start'] = pd.to_numeric(coordinates[1].values)
                    genes['stop'] = pd.to_numeric(coordinates[2].values)
                    genes = genes[['chrom', 'start', 'stop','entity_name','confidence_level','penetrance','mode_of_inheritance']]
                    genes.sort_values(by=['chrom','start'], inplace=True)
        
        if len(genes.index) > 0:
            return True, genes
        else:
            return False, genes
    
    def dumpPanels(self, output_dir, level=3, build="GRCh38"):
        index_file = output_dir+"/Index_table_"+self.time+".tsv"
        index_df = self._panels
        index_df['n_green'] = 0
        for panel_id in list(self._panels.index):
            print("Saving panel ", panel_id)
            genes_file = output_dir+"/"+build+"_Panel_"+str(panel_id)+".bed"
            success, genes = self.getGenes(pid = panel_id, level=0, build=build, out_format="detailed_bed")
            if success:
                try:
                    n_green = int(genes.groupby('confidence_level').count().loc[3,'entity_name'])
                except:
                    n_green = 0
                genes = genes[genes.confidence_level >= level]
                header = list(genes.columns)
                header[0] = "#" + header[0]
                genes.to_csv(genes_file,sep="\t",index=False,header=header)
                index_df.loc[panel_id,'n_green'] = n_green
        index_df = index_df[['name', 'disease_group', 'disease_sub_group', 'version', 'version_created', 'n_genes', 'n_green', 'relevant_disorders']]
        index_df.sort_values(by=['name'], inplace=True)
        index_df.to_csv(index_file, sep="\t")

#Simply store and search in ClinVar pathogenic genes
class ClinVar():
    def __init__(self, pathogenic_genes):
        to_list = lambda x: x.split(";")
        self._genes = pd.read_csv(pathogenic_genes, sep="\t", names = ['gene', 'conditions'], converters={'conditions': to_list})
 
    def getGenes(self, term, exact=False):
        if exact:
            mask = self._genes.conditions.apply(lambda x: term in x)       
            return self._genes[mask]
        else:
            r = re.compile(term, re.IGNORECASE)
            mask = []
            for line in list(self._genes.conditions):
                mask.append(any(r.search(x) for x in line))

            return self._genes[mask]

#Wrapper around GADO prioritizer tool. Process both list of HPOs or a GADO input file
class GADO():
    def __init__(self, gado_jar=GADO_JAR, gado_data=GADO_DATA):
        self._jar = gado_jar
        self._predictions_info = gado_data + '/hpo_predictions_info_01_02_2018.txt'
        self._genes = gado_data + '/hpo_predictions_genes_01_02_2018.txt'
        self._predictions = gado_data + '/hpo_predictions_sigOnly_spiked_01_02_2018'
    
    def process(self, hpo_ids, output_dir, hpo_ontology, sampleID="mysample"):
        obofile = hpo_ontology._obofile
        outputfile = output_dir + "/gado_tmp_process"
        inputfile = output_dir + "/gado_tmp_input"
        
        if isinstance(hpo_ids, list):        
            new_HPOs = hpo_ontology.update_obsolete(hpo_ids)
            with open(inputfile, 'w+') as f:
                f.write(sampleID + "\t" + "\t".join(new_HPOs) + "\n")
            f.close()
        
        elif os.path.isfile(hpo_ids):
            tmp_out = open(inputfile, 'w+')
            with open(hpo_ids) as f:
                line = f.readline()
                while line:
                    line = line.rstrip('\n')
                    line = line.split('\t')
                    old_HPOs = line[1:]
                    new_HPOs = hpo_ontology.update_obsolete(old_HPOs)
                    tmp_out.write(line[0] + '\t' + "\t".join(new_HPOs) + "\n")
                    line = f.readline()
            tmp_out.close()
        
        exitcode, unused_stdout, stderr = run(['java','-jar',self._jar,
                                '--mode', 'PROCESS',
                                '--output', outputfile,
                                '--caseHpo', inputfile,
                                '--hpoOntology', obofile,
                                '--hpoPredictionsInfo', self._predictions_info])
        
        return exitcode, stderr.decode('utf-8'), outputfile
    
    def prioritize(self, input_file, output_dir):
        exitcode, unused_stdout, stderr = run(['java','-jar',self._jar,
                                '--mode', 'PRIORITIZE',
                                '--output', output_dir,
                                '--caseHpoProcessed', input_file,
                                '--genes', self._genes,
                                '--hpoPredictions', self._predictions])
        
        return exitcode, stderr.decode('utf-8')
