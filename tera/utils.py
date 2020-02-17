"""
Utilities used by other modules.
"""
from SPARQLWrapper import SPARQLWrapper, JSON
from functools import wraps
from rdflib import Literal
from collections import defaultdict
import warnings
from tqdm import tqdm
from quantulum3 import parser 
from itertools import combinations

nan_values = ['nan', float('nan'),'--','-X','NA','NC',-1,'','sp.', -1,'sp,','var.','variant','NR']

unit_lookup = defaultdict(lambda: '')

unit_lookup.update({'mg':'Milligram',
                    'ug': 'Microgram',
                    'kg':'Kilogram',
                    'mM':'Millimol',
                    'ng':'Nanogram',
                    'g':'Gram',
                    'µg':'Microgram',
                    'L':'Litre',
                    '%':'Percent',
                    'cm':'Centimetre',
                    'mm':'Millimetre',
                    'nm':'Nanometre',
                    'deg':'Degree',
                    'C':'Celcius',
                    'K':'Kelvin',
                    'l':'Litre',
                    'psu':'PracticalSalinityUnit',
                    'h':'Hour',
                    'd':'Day',
                    'w':'Week'
                    }
                   )

def unit_parser(string):
    """
    Takes a unit string and converts to UNIT namespace string. 
    eg. mg/L -> MilligramPerLitre
    Filters out assumed missprints, eg. mg%/L -> MilligramPerLitre.
    
    Parameters
    ----------
    string : str
        Unit string.
    
    Returns
    -------
    str 
    """
    if len(string) < 2 and string not in unit_lookup:
        return ''
    
    if '/' in string:
        a,b = string.split('/', 1)
        return unit_parser(a) + 'Per' + unit_parser(b)
    if '-1' in string:
        return unit_parser(string.replace('-1','/'))
    
    if ' ' in string:
        a,b = string.split(' ', 1)
        return unit_parser(a) + '-' + unit_parser(b)
        
    if string in unit_lookup:
        return unit_lookup[string]
    
    else:
        res1 = [string[x:y] for x, y in combinations(range(len(string) + 1), r = 2)]
        res1.remove(string)
        res2 = map(unit_parser, res1)
        res = zip(res2,res1)
        res = [(a,b) for a,b in res if len(a) > 1]
        u,_ = sorted(res, key=lambda x:len(x[1]),reverse=True).pop(0)
        return u
    
    return ''
        

def tanimoto(fp1, fp2):
    """
    Calculate tanimoto similarity between two chemical fingerprints.
    
    Parameters 
    ----------
    fp1 : str
        Chemical fingerprint on binary form.
        
    fp2 : str
        Chemical fingerprint on binary form.
    
    Returns
    -------
    float
    """
    fp1_count = fp1.count('1')
    fp2_count = fp2.count('1')
    both_count = (fp1 & fp2).count('1')
    return float(both_count) / (fp1_count + fp2_count - both_count)


def test_endpoint(endpoint):
    """
    Test SPARQL endpoint.
    
    Parameters 
    ----------
    endpoint : str 
        SPARQL endpoint URL. ex: https://query.wikidata.org/sparql 
    
    Returns
    -------
    bool
    """
    sparql = SPARQLWrapper(endpoint)
    q = """
        SELECT ?s ?p ?o
        WHERE {?s ?p ?o}
        LIMIT 100
    """ 

    sparql.setQuery(q)
    sparql.setReturnFormat(JSON)
    try:
        results = sparql.query().convert()
        return True
    except:
        return False
    
    
def query_endpoint(endpoint, q, var = 'p'):
    """
    Wrapper for quering SPARQL endpoint with SPARQLWrapper.
    
    Parameters 
    ----------
    endpoint : str
        SPARQL endpoint URL. 
    
    q : str 
        SPARQL query. 
        
    var : str or list 
        Query variables to return.
    
    Returns
    -------
    set 
        Set of tuple query results. Tuple in order specified in input var. 
    """
    if not isinstance(var, list):
        var = [var]
        
    sparql = SPARQLWrapper(endpoint)
    
    out = {}
    try:
        sparql.setQuery(q)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for v in var:
            try:
                out[v] = [r[v]['value'] for r in results['results']['bindings']]
            except KeyError:
                out[v] = [None] * len(results['results']['bindings'])
        return set(zip(*[out[k] for k in out]))
    except Exception as e:
        print(e)
        warnings.warn('Query failed:\n' + q, UserWarning)
        return set()


def query_graph(graph, q):
    """
    Query rdflib.Graph. 
    
    Parameters 
    ----------
    graph : rdflib.Graph 
    
    q : str 
        SPARQL query
   
    Returns
    -------
    set 
    """
    try:
        return set(graph.query(q))
    except Exception as e:
        return set()

def prefixes(initNs):
    """
    Format prefixes for SPARQL. 
    
    Parameters 
    ----------
    initNs : dict 
        ex : {'ex':'http://example.org'} 
    
    Returns
    -------
    str 
    """
    q = ''
    for k,i in initNs.items():
        q += "PREFIX\t"+k+':\t' + '<'+str(i)+'>\n'
    return q

def strip_namespace(string, var = ['/']):
    """
    Remove namespace from URI.
    
    Parameters 
    ----------
    string : str 
        URI 
    var : str or list 
        Symbols to split string. ex. / or #.
    
    Returns
    -------
    str
    """
    if not isinstance(var,list):
        var = [var]
    tmp1 = str(string)
    for v in var:
        tmp2 = str(string).split(v)[-1]
        if len(tmp2) < len(tmp1):
            tmp1 = tmp2
    return tmp1

def do_recursively_in_class(func):
    """Enables function to take either element or iterable as input.
    
    Returns
    -------
    function
    """
    @wraps(func)
    def call_recursively(my_class_instance, x, **kwargs):
        if isinstance(x, (list,set,tuple)):
            f = lambda x: func(my_class_instance, x, **kwargs)
            out = {}
            pbar = lambda x: x
            if hasattr(my_class_instance, 'verbose'):
                if my_class_instance.verbose:
                    pbar = lambda x: tqdm(x)
            return dict(zip(x,map(f,pbar(x))))
        else:
            return func(my_class_instance, x, **kwargs)
        
    return call_recursively


def graph_to_dict(graph):
    """
    Map entities in graph to a dict.
    
    Parameters 
    ----------
    graph : rdflib.Graph
    
    Returns
    -------
    dict 
        On the form {entity : list of literals connected to entity}
    """
    entities = graph.subjects()
    d = defaultdict(list)
    
    for e in entities:
        d[e] = [str(o) for o in graph.objects(subject=e) if isinstance(o,Literal)]
    return d


