"""
A set of classes for aggregation of TERA data sources into common formats.
"""

from rdflib import Graph, Namespace, Literal, URIRef, BNode
from rdflib.namespace import RDF, OWL, RDFS
UNIT = Namespace('http://qudt.org/vocab/unit#')
import pandas as pd
import validators
import glob
import math
from tqdm import tqdm
import warnings

nan_values = ['nan', float('nan'),'--','-X','NA','NC',-1,'','sp.', -1,'sp,','var.','variant','NR']

class DataObject:
    def __init__(self, namespace = 'http://www.example.org/', verbose = True, name = 'Data Object'):
        """
        Base class for aggregation of data.
        
        Parameters
        ----------
        namespace : str 
            Base URI for the data set.
            
        verbose : bool
        """
        self.graph = Graph()
        self.namespace = Namespace(namespace)
        self.name = name 
        self.verbose = verbose
        
    def __str__(self):
        return self.name 
        
    def __dict__(self):
        return {
                'namespace':self.namespace,
                'num_triples':len(self.graph)
            }
    
    def __del__(self):
        self.graph = Graph()
    
    def save(self, path):
        """Save graph to file.
        
        Parameters
        ----------
        path : str 
            ex: file.nt
        """
        self.graph.serialize(path, format=path.split('.').pop(-1))
        
    def replace(self, converted):
        """Replace old entities with new in data object. 
        Usefull after converting between datasets.
        
        Parameters
        ----------
        converted : list
            list of (old, new) tuples. 
            
        """
        if len(converted) < 1:
            warnings.warn('Empty mapping list.')
            return
    
        tmp = set()
        for old, new in converted:
            triples = self.graph.triples((old,None,None))
            tmp |= set([(new,p,o) for _,p,o in triples])
            triples = self.graph.triples((None, None, old))
            tmp |= set([(s,p,new) for s,p,_ in triples])
            self.graph.remove((old,None,None))
            self.graph.remove((None,None,old))
            
        for t in tmp:
            self.graph.add(t)
            
    def apply_func(self, func, dataframe, cols, sub_bar=False):
        pbar = None
        if self.verbose and not sub_bar:
            pbar = tqdm(total=len(dataframe.index),desc=self.name)
            
        for row in zip(*[dataframe[c] for c in cols]):
            func(row)
            if pbar: pbar.update(1)
            

class Taxonomy(DataObject):
    def __init__(self, 
                 namespace = 'https://www.ncbi.nlm.nih.gov/taxonomy/',
                 name = 'NCBI Taxonomy',
                 verbose = True, 
                 directory = None):
        """
        Aggregation of the NCBI Taxonomy. 
        
        Parameters
        ---------- 
        directory : str 
            Path to data set. Downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
            
        """
        super(Taxonomy, self).__init__(namespace, verbose, name)
        
        if directory:
            self._load_ncbi_taxonomy(directory)
        
        self.verbose = verbose
    
    def _load_ncbi_taxonomy(self, directory):
        self._load_hierarchy(directory+'nodes.dmp')
        self._load_divisions(directory+'division.dmp')
        self._load_names(directory+'names.dmp')

    def _load_hierarchy(self, path):
        df = pd.read_csv(path, sep='|', usecols=[0,1,2,4], names=['child','parent','rank','division'], na_values = nan_values, dtype = str)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())

        def func(row):
            c,p,r,d = row
            c = self.namespace['taxon/'+str(c)]
            self.graph.add((c,RDF.type,self.namespace['Taxon']))
            rc = r
            r = r.replace(' ','_')
            self.graph.add((c, self.namespace['rank'], self.namespace['rank/'+r]))
            self.graph.add((self.namespace['rank/'+r], RDFS.label, Literal(rc)))
            self.graph.add((self.namespace['rank/'+r], RDF.type, self.namespace['Rank']))
            p = self.namespace['taxon/'+str(p)]
            self.graph.add((c, RDFS.subClassOf, p))
            d = str(d).replace(' ','_')
            d = self.namespace['division/'+str(d)]
            self.graph.add((c, RDFS.subClassOf, d))
        
        self.apply_func(func, df, ['child','parent','rank','division'])
    
    def _load_names(self, path):
        df = pd.read_csv(path, sep='|', usecols=[0,1,2,3], names=['taxon','name','unique_name','name_type'],na_values = nan_values,dtype = str)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())

        def func(row):
            c,n,un,nt = row
            c = self.namespace['taxon/'+str(c)]
            n = Literal(n)
            un = Literal(un)
            self.graph.add((c, RDFS.label, un))
            ntl = Literal(nt)
            nt = self.namespace[nt.replace(' ','_')]
            self.graph.add((c,nt,n))
            self.graph.add((nt,RDFS.label,ntl))
        
        self.apply_func(func, df, ['taxon','name','unique_name','name_type'])
        
    def _load_divisions(self, path):
        df = pd.read_csv(path, sep='|', usecols=[0,1,2], names=['division','acronym','name'], na_values = nan_values, dtype = str)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())

        def func(row):
            d,a,n = row
            d = self.namespace['division/'+str(d)]
            self.graph.add((d,RDF.type,self.namespace['Division']))
            self.graph.add((d,RDFS.label,Literal(n)))
            self.graph.add((d,RDFS.label,Literal(a)))
        
        self.apply_func(func, df, ['division','acronym','name'])
    
class Traits(DataObject):
    def __init__(self, 
                 namespace = 'https://eol.org/pages/',
                 name = 'EOL Traits',
                 verbose = True,
                 directory = None):
        """
        Encyclopedia of Life Traits. 
        
        Parameters
        ----------
        directory : str  
            Path to data set. See https://opendata.eol.org/dataset/all-trait-data-large
            
        """
        super(Traits, self).__init__(namespace, verbose, name)
        
        if directory:
            self._load_eol_traits(directory)
    
    def _load_eol_traits(self, directory):
        self._load_traits(directory+'trait_bank/traits.csv')
        for f in glob.glob(directory+'eol_rels/*.csv'):
            self._load_eol_subclasses(f)
    
    def _load_traits(self, path):
        df = pd.read_csv(path, sep=',', usecols=['page_id','predicate','value_uri'], na_values = nan_values, dtype=str)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())

        def func(row):
            s,p,o = row
            s = self.namespace[s]
            
            try:
                val = validators.url(o)
                o = URIRef(o)
            except TypeError:
                o = Literal(o)
                val = True

            if validators.url(s) and validators.url(p) and val:
                self.graph.add((URIRef(s),URIRef(p),o))

        self.apply_func(func, df, ['page_id','predicate','value_uri'])            
            
    def _load_eol_subclasses(self, path):
        try: 
            try:
                df = pd.read_csv(path,sep=',',usecols=['child','parent'],na_values = nan_values, dtype=str)
                df.dropna(inplace=True)
                df = df.apply(lambda x: x.str.strip())
            except ValueError:
                df = pd.read_csv(path,sep=',',header=None,na_values = nan_values, dtype=str)
                df.columns = ['parent','child']
                df.dropna(inplace=True)
                df = df.apply(lambda x: x.str.strip())
            
        except FileNotFoundError as e:
            print(e,path)
        
        def func(row):
            c,p = row
            if validators.url(c) and validators.url(p):
                c,p = URIRef(c),URIRef(p)
                self.graph.add((c,RDFS.subClassOf,p))
        
        self.apply_func(func, df, ['child','parent'])
        
    
class Effects(DataObject):
    def __init__(self, 
                    namespace = 'https://cfpub.epa.gov/ecotox/',
                    name = 'Ecotox Effects',
                    verbose = True,
                    directory = None):
        """
        Ecotox effects data aggregation.
        
        Parameters
        ---------- 
        directory : str 
            Path to data set. Downloaded from ftp://newftp.epa.gov/ecotox/ecotox_ascii_12_12_2019.exe
        """
        super(Effects, self).__init__(namespace, verbose, name)
        
        self._load_effect_data(directory + 'tests.txt', directory + 'results.txt')
        
    def _load_effect_data(self, tests_path, results_path):
        tests = pd.read_csv(tests_path, sep='|', dtype = str, na_values = nan_values)
        tests.dropna(inplace=True, subset=['test_id',
                        'test_cas',
                        'species_number'])
        tests.fillna(inplace=True, value='missing')
        tests = tests.apply(lambda x: x.str.strip())
        results = pd.read_csv(results_path, sep='|',  dtype = str, na_values = nan_values)
        results.dropna(inplace=True, subset=['test_id','endpoint','conc1_mean','conc1_unit','effect'])
        results.fillna(inplace=True, value='missing')
        results = results.apply(lambda x: x.str.strip())

        def test_func(row):
            test_id, cas_number, species_number, stdm, stdu, habitat, lifestage, age, ageunit, weight, weightunit = row
            
            #must be included
            t = self.namespace['test/'+str(test_id)]
            s = self.namespace['taxon/'+str(species_number)]
            c = self.namespace['cas/'+str(cas_number)]
            self.graph.add((t, RDF.type, self.namespace['Test']))
            self.graph.add((t, self.namespace['species'], s))
            self.graph.add((t, self.namespace['chemical'], c))

            for v,u,p in zip([stdm,age,weight],[stdu,ageunit,weightunit],['studyDuration','organismAge','organismWeight']):
                if v != 'missing':
                    b = BNode()
                    self.graph.add( (b, RDF.value, Literal(v)) )
                    if u != 'missing':
                        self.graph.add( (b, UNIT.units, Literal(u)) )
                    self.graph.add( (t, self.namespace[p], b) )
            
            if habitat != 'missing':
                self.graph.add((t, self.namespace['organismHabitat'],self.namespace['habitat/'+habitat]))
            if lifestage != 'missing':
                self.graph.add((t, self.namespace['organismLifestage'],self.namespace['lifestage/'+lifestage]))

        def results_func(row):
            test_id, endpoint, conc, conc_unit, effect = row
            t = self.namespace['test/'+str(test_id)]
            
            r = BNode()
            ep = self.namespace['endpoint/'+str(endpoint)]
            ef = self.namespace['effect/'+str(effect)]
            
            self.graph.add((r,self.namespace['endpoint'],ep))
            self.graph.add((r,self.namespace['effect'],ef))
            b = BNode()
            self.graph.add( (b, RDF.value, Literal(conc)) )
            if conc_unit != 'missing':
                #TODO creating mapping for units not in UNIT.units
                self.graph.add( (b, UNIT.units, Literal(conc_unit)) )
                
            self.graph.add( (r, self.namespace['concentration'], b) )
            
            self.graph.add((t,self.namespace['hasResult'],r))
            
        self.apply_func(test_func, tests, ['test_id',
                        'test_cas',
                        'species_number',
                        'study_duration_mean',
                        'study_duration_unit',
                        'organism_habitat',
                        'organism_lifestage',
                        'organism_age_mean',
                        'organism_age_unit',
                        'organism_init_wt_mean',
                        'organism_init_wt_unit'])
        
        self.apply_func(results_func, results, ['test_id','endpoint','conc1_mean','conc1_unit','effect'])


class EcotoxTaxonomy(DataObject):
    def __init__(self, 
                    namespace = 'https://cfpub.epa.gov/ecotox/',
                    name = 'Ecotox Taxonomy',
                    verbose = True,
                    directory = None):
        """
        Ecotox taxonomy aggregation. 
        
        Parameters
        ----------
        directory : str 
            Path to dataset. Downloaded from ftp://newftp.epa.gov/ecotox/ecotox_ascii_12_12_2019.exe
        """
        super(EcotoxTaxonomy, self).__init__(namespace, verbose, name)
        
        self._load_species(directory + 'validation/species.txt' )
        self._load_synonyms(directory + 'validation/species_synonyms.txt')
        self._load_hierarchy(directory + 'validation/species.txt')
        
    def _load_species(self, path):
        df = pd.read_csv(path, sep='|',  dtype = str, na_values = nan_values)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())
        
        def func(row):
            s, cn, ln, group = row
            
            s = self.namespace['taxon/'+s]
            names = group.split(',')
            group = group.replace('/','')
            group = group.replace('.','')
            group = group.replace(' ','')
            tmp = group.split(',')
            group_uri = [self.namespace['group/'+gr] for gr in tmp]
            
            for gri,n in zip(group_uri,names):
                self.graph.add((s, RDFS.subClassOf, gri))
                self.graph.add((gri, RDFS.label, Literal(n)))
                self.graph.add((gri, RDF.type, self.onto_namespace['SpeciesGroup']))
                
            self.graph.add((s, RDF.type, self.namespace['Taxon']))
            self.graph.add((s, self.namespace['commonName'], Literal(cn)))
            self.graph.add((s, self.namespace['latinName'], Literal(ln)))
        
        self.apply_func(func, df, ['species_number','common_name','latin_name','ecotox_group'])
            
    def _load_synonyms(self, path):
        df = pd.read_csv(path, sep='|',  dtype = str, na_values = nan_values)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())
        
        def func(row):
            s, ln = row
            s = self.namespace['taxon/'+s]
            self.graph.add((s, self.namespace['latinName'], Literal(ln)))
        
        self.apply_func(func, df, ['species_number','latin_name'])
            
    
    def _load_hierarchy(self, path):
        df = pd.read_csv(path, sep= '|', dtype = str)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())
        
        def func(row):
            sn, lineage = row
            curr = sn
            for l in filter(lambda x: x != '', lineage):
                curr = self.namespace['taxon/'+str(curr)]
                self.graph.add((curr, RDFS.subClassOf, self.namespace['taxon/'+str(l)]))
                curr = l
        
        pbar = None
        if self.verbose: pbar = tqdm(total=len(df.index))
        for row in zip(df['species_number'],
                        zip(df['variety'],
                            df['subspecies'],
                            df['species'],
                            df['genus'],
                            df['family'],
                            df['tax_order'],
                            df['class'],
                            df['superclass'],
                            df['subphylum_div'],
                            df['phylum_division'],
                            df['kingdom'])):
            func(row)
            if pbar: pbar.update(1)
        
        
class EcotoxChemicals(DataObject):
    def __init__(self, 
                    namespace = 'https://cfpub.epa.gov/ecotox/',
                    name = 'Ecotox Chemicals',
                    verbose = True,
                    directory = None):
        """
        Ecotox chemicals aggregation. 
        
        Parameters
        ----------
        directory : str 
            Path to dataset. Downloaded from ftp://newftp.epa.gov/ecotox/ecotox_ascii_12_12_2019.exe
        """
        super(EcotoxChemicals, self).__init__(namespace, verbose, name)
        
        self._load_chemicals(directory + 'validation/chemicals.txt')
        
    def _load_chemicals(self, path):
        df = pd.read_csv(path, sep='|',  dtype = str, na_values = nan_values)
        df.dropna(inplace=True)
        df = df.apply(lambda x: x.str.strip())
        
        def func(row):
            c, n, group = row
            c = self.namespace['cas/'+str(c)]
            self.graph.add((c, RDF.type, self.namespace['Chemical']))
            for a in n.split(', '):
                self.graph.add((c, RDFS.label, Literal(a)))
            
            names = group.split(',')
            group = group.replace('/','')
            group = group.replace('.','')
            group = group.replace(' ','')
            tmp = group.split(',')
            group_uri = [self.namespace['group/'+gr] for gr in tmp]
        
            for gri,n in zip(group_uri,names):
                self.graph.add((c, RDFS.subClassOf, gri))
                self.graph.add((gri, RDFS.label, Literal(n)))
                self.graph.add((gri, RDF.type, self.namespace['ChemicalGroup']))
        
        self.apply_func(func, df, ['cas_number','chemical_name','ecotox_group'])
        
class PubChem(DataObject):
    def __init__(self, 
                    namespace = 'http://rdf.ncbi.nlm.nih.gov/pubchem/compound/',
                    name = 'PubChem',
                    directory = None):
        """
        PubChem data loading.
        
        Parameters
        ----------
        directory : str 
            Path to turtle RDF files. Downloaded from ftp://ftp.ncbi.nlm.nih.gov/pubchem/RDF/ . Used files: compound/pc_compound_type.ttl and compound/pc_compound2parent.ttl .
        """
        super(PubChem, self).__init__(namespace, name)
        
        for f in glob.glob(directory+'*.ttl'):
            self._load_data(f)
        
    def _load_data(self, path):
        self.graph += Graph().parse(path, format = 'ttl')

class ChEBI(DataObject):
    def __init__(self, 
                    namespace = 'http://purl.obolibrary.org/obo/',
                    name = 'ChEBI',
                    directory = None):
        """
        ChEBI data loading. 
        
        Parameters
        ---------- 
        directory : str 
            Path to turtle RDF files. See https://www.ebi.ac.uk/rdf/datasets/
        """
        super(ChEBI, self).__init__(namespace, name)
        
        for f in glob.glob(directory+'*.ttl'):
            self._load_data(f)
        
    def _load_data(self, path):
        self.graph += Graph().parse(path, format = 'ttl')
        
class MeSH(DataObject):
    def __init__(self, 
                    namespace = 'http://id.nlm.nih.gov/mesh/',
                    name = 'MeSH',
                    directory = None):
        """
        MeSH data loading. 
        
        Parameters
        ----------
        directory : str  
            Path to nt RDF files. See https://id.nlm.nih.gov/mesh/
        """
        super(MeSH, self).__init__(namespace, name)
        
        for f in glob.glob(directory+'*.nt'):
            self._load_data(f)
        
    def _load_data(self, path):
        self.graph += Graph().parse(path, format = 'nt')
        
