
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
        self.name = name
    
    def mapping(self, x, reverse = False):
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
        if strip:
            id_ = strip_namespace(id_,['/','#','CID'])
        return self.mapping(id_,reverse)

class EndpointMapping(Alignment):
    def __init__(self, endpoint):
        super(EndpointMapping, self).__init__()
        self.mappings = self.load_mapping(endpoint)
    
    def load_mapping(self, endpoint):
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
        query :: str
            Wikidata query with two variables. 
            ex: from inchikey to cas
            SELECT ?from ?to
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P231 ?to .
            } 
        """
        super(WikidataMapping, self).__init__()
        self.mappings = self.load_mapping(query)
        
    def load_mapping(self, query):
        res = query_endpoint('https://query.wikidata.org/sparql', 
                             query, 
                             var = ['from', 'to'])
        return {str(f):str(t) for f,t in res}

class LogMapMapping(Alignment):
    def __init__(self, filename, threshold=0.95):
        """
        filename :: str
            path to logmap output file (.rdf)
        """
        super(LogMapMapping, self).__init__()
        
        self.threshold = threshold
        self.mappings = self.load_mapping(filename)
    
    def load_mapping(self, filename):
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
    def __init__(self, dict1, dict2):
        """
        g1 :: dict 
        g2 :: dict
            {entity:list of strings}
        """
        super(StringMatchingMapping, self).__init__()
        
        self.threshold = 0.95
        self.mappings = self.load_mapping(dict1,dict2)
    
    def load_mapping(self, dict1, dict2):
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
    def __init__(self, g1, g2):
        """
        g1 :: Graph
        g2 :: Graph
        Uses literal to align entities in g1 and g2.
        """
        super(StringGraphMapping, self).__init__()
        
        self.threshold = 0.95
        dict1 = graph_to_dict(g1)
        dict2 = graph_to_dict(g2)
        self.mappings = self.load_mapping(dict1, dict2)
    
    def load_mapping(self, dict1, dict2):
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
        super(NCBIToEcotox, self).__init__(dataobject1.graph,
                                           dataobject2.graph)
        
        
