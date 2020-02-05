"""
A set of classes for aligning data aggregated with tools in DataAggregation.
"""


from rdflib import Graph, Namespace, Literal, URIRef
from rdflib.namespace import RDF, OWL, RDFS
import pandas as pd
import validators
from .utils import query_endpoint, strip_namespace, do_recursively_in_class, graph_to_dict
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from collections import defaultdict

class Alignment:
    def __init__(self, name = 'Alignment'):
        """Base class for alignment of two data sets. 
        
        Parameters
        ----------
        name : str        
        """
        self.name = name
    
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
        tmp = self.mappings
        if reverse:
            tmp = {tmp[k]:k for k in tmp}
        x = str(x)
        if x in tmp:
            return tmp[x]
        else:
            return 'no mapping'
    
    @do_recursively_in_class
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
            id_ = strip_namespace(str(id_),['/','#','CID'])
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
        self.mappings = self._load_mapping(endpoint)
    
    def _load_mapping(self, endpoint):
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
        res = query_endpoint(endpoint, query, var = ['s','o'])
        return {str(s):str(o) for s,o in res}

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
        self.mappings = self._load_mapping(query)
        
    def _load_mapping(self, query):
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
        res = query_endpoint('https://query.wikidata.org/sparql', 
                             query, 
                             var = ['from', 'to'])
        return {str(f):str(t) for f,t in res}

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
        self.mappings = self._load_mapping(filename)
    
    def _load_mapping(self, filename):
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
        g.parse(filename)
        o = URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmentCell')
        for s in g.subjects(predicate=RDF.type, object = o):
            e1 = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmententity1'))).pop(0)
            e2 = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmententity2'))).pop(0)
            score = list(g.objects(subject=s,predicate=URIRef('http://knowledgeweb.semanticweb.org/heterogeneity/alignmentmeasure'))).pop(0)
            
            score = float(score)
            if score >= self.threshold:
                out[str(e1)] = str(e2)
        
        return out

        
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
        self.mappings = self._load_mapping(dict1,dict2)
    
    def _load_mapping(self, dict1, dict2):
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
        for k1 in dict1:
            for k2 in dict2:
                try:
                    _, score = process.extractOne(dict1[k1],dict2[k2])
                except TypeError:
                    score = 0
                    
                if score >= self.threshold:
                    tmp[k1,k2] = max(tmp[k1,k2],score)
        
        return {k1:k2 for k1,k2 in tmp}
    
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
        dict1 = graph_to_dict(g1)
        dict2 = graph_to_dict(g2)
        self.mappings = self._load_mapping(dict1, dict2)
    
    def _load_mapping(self, dict1, dict2):
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
        for k1 in dict1:
            for k2 in dict2:
                try:
                    _, score = process.extractOne(dict1[k1],dict2[k2])
                except TypeError:
                    score = 0
                    
                if score >= self.threshold:
                    tmp[k1,k2] = max(tmp[k1,k2],score)
        
        return {k1:k2 for k1,k2 in tmp}
        

class InchikeyToCas(WikidataMapping):
    def __init__(self):
        """Class which creates inchikey to cas mapping."""
        query = """
        SELECT ?from ?to
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
        SELECT ?from ?to
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
        SELECT ?from ?to
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
        SELECT ?from ?to
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
        SELECT ?from ?to
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
        SELECT ?from ?to
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
        





