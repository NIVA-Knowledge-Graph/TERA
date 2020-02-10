"""
A set of classes for aligning data aggregated with tools in DataAggregation.
"""
from rdflib import Graph, Namespace, Literal, URIRef
from rdflib.namespace import RDF, OWL, RDFS
import pandas as pd
import validators
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from collections import defaultdict
import tera.utils as ut

class Alignment:
    def __init__(self, name = 'Alignment'):
        """Base class for alignment of two data sets. 
        
        Parameters
        ----------
        name : str        
        """
        self.name = name
    
    def load(self):
        raise NotImplementedError
    
    def _to_defaultdict(self):
        self.mappings = defaultdict('no mapping',self.mappings)
    
    def _mapping(self, x, reverse = False):
        """
        Maps from one id type to another. 
        
        Parameters
        ----------
        x : rdflib.URIRef or str 
            URI/identifier to map from. 
            
        reverse : bool 
            Reverse the direction of mapping. 
                
        Returns
        -------
        str 
            If no mapping exists, returns 'no mapping'
        """
        if not hasattr(self, 'mappings'):
            self.load()
            
        tmp = self.mappings
        if reverse:
            tmp = defaultdict('no mapping')
            for k,i in self.mappings.items():
                k,i = str(k),str(i)
                if not i == 'no mapping':
                    tmp[i] = k
        x = str(x)
        return tmp[x]
        
    
    def convert(self, id_, reverse=True, strip = False):
        """
        Convert a set of ids into new identifiers.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, set 
            URI(s)/identifier(s)  
        
        reverse : bool 
            Reverse the direction of mapping. 
            
        strip : bool 
            Remove namespace.
                
        Returns
        -------
        str or dict
            Mapped values.
        """
        if strip:
            id_ = ut.strip_namespace(str(id_),['/','#','CID'])
        return self._mapping(id_,reverse)

class EndpointMapping(Alignment):
    def __init__(self, endpoint):
        super(EndpointMapping, self).__init__()
        """Class for loading mappings based on owl:sameAs property.
        
        Parameters
        ----------
        endpoint : str 
            SPARQL endpoint URL.
        """
        self.endpoint = endpoint
    
    def load(self):
        """
        Load mappings from endpoint.
        
        Parameters
        ----------
        endpoint : str 
            SPARQL endpoint URL
                
        Returns
        -------
        dict
            On the form {from:to}
        """
        query = """
        SELECT ?s ?o WHERE {
            ?s <http://www.w3.org/2002/07/owl#sameAs> ?o .
        } 
        """
        res = ut.query_endpoint(self.endpoint, query, var = ['s','o'])
        self.mappings = {str(s):str(o) for s,o in res}
        self._to_defaultdict()

class WikidataMapping(Alignment):
    def __init__(self, query):
        """
        Class for loading mappings from wikidata.
        
        Parameters
        ----------
        query : str 
            Wikidata query with two variables. 
            
            eg. from inchikey to cas: 
            
            SELECT ?from ?to {  
            ?compound wdt:P235 ?from . 
            ?compound wdt:P231 ?to .} 
        """
        super(WikidataMapping, self).__init__()
        self.query = query
        
    def load(self):
        """
        Load mappings from wikidata.
        
        Parameters
        ----------
        query : str 
            wikidata query. 
                
        Returns
        -------
        dict
            On the form {from:to}
        """
        res = ut.query_endpoint('https://query.wikidata.org/sparql', 
                             self.query, 
                             var = ['from', 'to'])
        self.mappings = {str(f):str(t) for f,t in res}
        self._to_defaultdict()

class LogMapMapping(Alignment):
    def __init__(self, filename, threshold=0.95):
        """
        Class for using LogMap (or other system) alignments. 
        
        Parameters
        ----------
        filename : str 
            Path to logmap output file (.rdf) 
        
        threshold : float 
            Alignment threshold.
        """
        super(LogMapMapping, self).__init__()
        
        self.threshold = threshold
        self.filename = filename
    
    def load(self):
        """
        Load mappings from alignment system.
        
        Parameters
        ----------
        filename : str 
            Path to file. 
                
        Returns
        -------
        dict
            On the form {from:to}
        """
        out = {}
        g = Graph()
        g.parse(self.filename)
        o = URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmentCell')
        for s in g.subjects(predicate=RDF.type, object = o):
            e1 = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmententity1'))).pop(0)
            e2 = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmententity2'))).pop(0)
            score = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmentmeasure'))).pop(0)
            
            score = float(score)
            if score >= self.threshold:
                out[str(e1)] = str(e2)
        
        self.mappings = out
        self._to_defaultdict()

        
class StringMatchingMapping(Alignment):
    def __init__(self, dict1, dict2, threshold = 0.95):
        """
        Class for creating mapping between two label dictonaries using string matching. 
        
        Parameters
        ----------
        dict1 : dict 
            Dictonary on the form {entity:list of labels} 
        
        dict2 : dict 
            Same as dict1.
            
        threshold : float 
            Alignment threshold.
        """
        super(StringMatchingMapping, self).__init__()
        
        self.threshold = threshold
        self.dict1 = dict1
        self.dict2 = dict2
    
    def load(self):
        """
        
        Parameters
        ----------
        dict1 : dict 
            Dictonary on the form {entity:list of labels} 
        
        dict2 : dict 
            Same as dict1.
            
        threshold : float 
            Alignment threshold.
                
        Returns
        -------
        dict 
            On the form {from:to}
        """
        tmp = defaultdict(float)
        for k1 in self.dict1:
            for k2 in self.dict2:
                try:
                    _, score = process.extractOne(self.dict1[k1],self.dict2[k2])
                except TypeError:
                    score = 0
                    
                if score >= self.threshold:
                    tmp[k1,k2] = max(tmp[k1,k2],score)
        
        self.mappings = {k1:k2 for k1,k2 in tmp}
        self._to_defaultdict()
    
class StringGraphMapping(Alignment):
    def __init__(self, g1, g2, threshold = 0.95):
        """
        Class for creating mapping between two graph using string matching. 
        
        Parameters
        ----------
        g1 : rdflib.Graph 
        
        g2 : rdflib.Graph 
        
        threshold : float 
            Alignment threshold.
                
        Returns
        -------
        dict 
            On the form {from:to}
        """
        super(StringGraphMapping, self).__init__()
        
        self.threshold = threshold
        self.g1 = g1
        self.g2 = g2
    
    def load(self):
        """
        Parameters
        ----------
        dict1 : dict 
            Dictonary on the form {entity:list of labels} 
        
        dict2 : dict 
            Same as dict1.
        
        threshold : float 
            Alignment threshold.
                
        Returns
        -------
        dict 
            On the form {from:to}
        """
        dict1 = ut.graph_to_dict(self.g1)
        dict2 = ut.graph_to_dict(self.g2)
        
        tmp = defaultdict(float)
        for k1 in dict1:
            for k2 in dict2:
                try:
                    _, score = process.extractOne(dict1[k1],dict2[k2])
                except TypeError:
                    score = 0
                    
                if score >= self.threshold:
                    tmp[k1,k2] = max(tmp[k1,k2],score)
        
        self.mappings = {k1:k2 for k1,k2 in tmp}
        self._to_defaultdict()
        

class InchikeyToCas(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to cas mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?compound wdt:P235 ?to .
            ?compound wdt:P231 ?from .
            }
        """
        super(InchikeyToCas, self).__init__(query=query)
    
class InchikeyToPubChem(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to pubchem mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P662 ?to .
            }
        """
        super(InchikeyToPubChem, self).__init__(query=query)
    
class InchikeyToChEBI(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to chebi mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P683 ?to .
            }
        """
        super(InchikeyToChEBI, self).__init__(query=query)

class InchikeyToChEMBL(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to chemble mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P592 ?to .
            }
        """
        super(InchikeyToChEMBL, self).__init__(query=query)
        
class InchikeyToMeSH(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to mesh mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P486 ?to .
            }
        """
        super(InchikeyToMeSH, self).__init__(query=query)

class NCBIToEOL(WikidataMapping):
    def __init__(self):
        """Class which creates ncbi to eol mapping."""
        query = """
        SELECT ?from ?to WHERE
            { 
            ?taxon wdt:P685 ?from .
            ?taxon wdt:P830 ?to .
            }
        """
        super(NCBIToEOL, self).__init__(query=query)
        
        
#TODO change ncbi -> ecotox mapping to concensus mappings.
class NCBIToEcotox(StringGraphMapping):
    def __init__(self, dataobject1, dataobject2):
        """Class which creates ncbi to ecotox mapping."""
        super(NCBIToEcotox, self).__init__(dataobject1.graph,
                                           dataobject2.graph)
        





