
from rdflib import Graph, Namespace, Literal, URIRef
from rdflib.namespace import RDF, OWL, RDFS
import pandas as pd
import validators
from utils import query_endpoint
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
    
    def convert(self, ids, reverse=True):
        return map(lambda x: self.mapping(x,reverse), ids)
    

class WikidataMapping(Alignment):
    def __init__(self, query, name='WikidataMapping'):
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
        super(WikidataMapping, self).__init__(name)
        self.mappings = self.load_mapping(query)
        
    def load_mapping(self, query):
        res = query_endpoint('https://query.wikidata.org/sparql', 
                             query, 
                             var = ['from', 'to'])
        return {str(f):str(t) for f,t in res}

class LogMapMapping(Alignment):
    def __init__(self, filename, threshold=0.95, name='LogMapMapping'):
        """
        filename :: str
            path to logmap output file (.rdf)
        """
        super(LogMapMapping, self).__init__(name)
        
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
    def __init__(self, g1, g2, threshold=0.99, name='StringMatchingMapping'):
        """
        g1 :: dict 
        g2 :: dict
            {entity:list of strings}
        """
        super(StringMatchingMapping, self).__init__(name)
        
        self.threshold = threshold
        self.mappings = self.load_mapping(g1,g2)
    
    def load_mapping(self, dict1, dict2):
        tmp = defaultdict(float)
        for k1 in dict1:
            for k2 in dict2:
                _, score = process.extractOne(dict1[k2],dict2[k2])
                if score >= self.threshold:
                    tmp[k1,k2] = max(tmp[k1,k2],score)
        
        return {k1:k2 for k1,k2 in tmp}
    
class StringGraphMapping(StringMatchingMapping):
    def __init__(self, g1, g2, threshold=0.99, name='StringGraphMapping'):
        """
        g1 :: Graph
        g2 :: Graph
        Uses literal to align entities in g1 and g2.
        """
        dict1 = self.graph_to_dict(g1)
        dict2 = self.graph_to_dict(g2)
        super(StringMatchingMapping, self).__init__(dict1,
                                                    dict2,
                                                    threshold=threshold,
                                                    name=name)
        
    def graph_to_dict(self, graph):
        entities = graph.subjects()
        d = defaultdict(list)
        
        for e in entities1:
            d[e] = [str(o) for g1.objects(subject=e) if isinstance(o,Literal)]
        return d
        
class CasToInchikey(WikidataMapping):
    def __init__(self, name = 'CasToInchikey'):
        query = """
        SELECT ?from ?to
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P231 ?to .
            }
        """
        super(CasToInchikey, self).__init__(query=query, 
                                            name=name)
    
    
class InchikeyToCas(WikidataMapping):
    def __init__(self, name = 'InchikeyToCas'):
        query = """
        SELECT ?from ?to
            { 
            ?compound wdt:P235 ?to .
            ?compound wdt:P231 ?from .
            }
        """
        super(InchikeyToCas, self).__init__(query=query, 
                                            name=name)
    
class InchikeyToPubChem(WikidataMapping):
    def __init__(self, name = 'InchikeyToPubChem'):
        query = """
        SELECT ?from ?to
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P662 ?to .
            }
        """
        super(InchikeyToPubChem, self).__init__(query=query, 
                                            name=name)
    
class InchikeyToChEBI(WikidataMapping):
    def __init__(self, name = 'InchikeyToChEBI'):
        query = """
        SELECT ?from ?to
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P683 ?to .
            }
        """
        super(InchikeyToChEBI, self).__init__(query=query, 
                                            name=name)

class InchikeyToChEMBL(WikidataMapping):
    def __init__(self, name = 'InchikeyToChEMBL'):
        query = """
        SELECT ?from ?to
            { 
            ?compound wdt:P235 ?from .
            ?compound wdt:P592 ?to .
            }
        """
        super(InchikeyToChEMBL, self).__init__(query=query, 
                                            name=name)

class NCBIToEOL(WikidataMapping):
    def __init__(self, name = 'NCBIToEOL'):
        query = """
        SELECT ?from ?to
            { 
            ?taxon wdt:P685 ?from .
            ?taxon wdt:P830 ?to .
            }
        """
        super(NCBIToEOL, self).__init__(query=query, 
                                        name=name)
    

    
    
        
        
