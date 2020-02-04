#utils 
from SPARQLWrapper import SPARQLWrapper, JSON
from functools import wraps
from rdflib import Literal
from collections import defaultdict

nan_values = ['nan', float('nan'),'--','-X','NA','NC',-1,'','sp.', -1,'sp,','var.','variant','NR']

def tanimoto(fp1, fp2):
    fp1_count = fp1.count('1')
    fp2_count = fp2.count('1')
    both_count = (fp1 & fp2).count('1')
    return float(both_count) / (fp1_count + fp2_count - both_count)


def remove_missing(item):
    if item in nan_values:
        item = 'NaN'
    return item

def test_endpoint(endpoint):
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
        print(e.message)
        return set()


def test_endpoint(endpoint):
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
    
def query_graph(graph, q):
    try:
        return list(graph.query(q))
    except Exception as e:
        print(e.message)
        return []
    

def prefixes(initNs):
    q = ''
    for k,i in initNs.items():
        q += "PREFIX\t"+k+':\t' + '<'+str(i)+'>\n'
    return q

def strip_namespace(string, var = ['/']):
    tmp1 = str(string)
    for v in var:
        tmp2 = str(string).split(v)[-1]
        if len(tmp2) < len(tmp1):
            tmp1 = tmp2
    return tmp1


def do_recursively_in_class(func):
    @wraps(func)
    def call_recursively(my_class_instance, x, **kwargs):
        if isinstance(x, (list,set,tuple)):
            return {j: func(my_class_instance, j, **kwargs) for j in x}
        else:
            return func(my_class_instance, x, **kwargs)
        
    return call_recursively


def graph_to_dict(graph):
    entities = graph.subjects()
    d = defaultdict(list)
    
    for e in entities:
        d[e] = [str(o) for o in graph.objects(subject=e) if isinstance(o,Literal)]
    return d


