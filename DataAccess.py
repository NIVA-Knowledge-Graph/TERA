
from rdflib import Graph, Namespace
from rdflib.namespace import RDF, RDFS, OWL
UNIT = Namespace('http://qudt.org/vocab/unit#')
from utils import query_endpoint, test_endpoint, query_graph, prefixes

class API:
    def __init__(self, namespace=None, endpoint=None, dataobject=None, name='API'):
        """
        endpoint :: str
            sparql endpoint url
        dataobject :: DataObject 
            see DataAggregation
        """
        if not endpoint and not dataobject:
            raise NotImplementedError
        
        if endpoint:
            if test_endpoint(endpoint):
                self.endpoint = endpoint
                self.use_endpoint = True
                self.namespace = Namespace(namespace)
            else:
                print('Did not reach endpoint at: ', endpoint)
                
        if dataobject:
            self.dataobject = dataobject
            self.use_endpoint = False
            self.namespace = dataobject.namespace
        
        self.name = name
            
        initNs = {'rdf':RDF, 'ns':self.namespace, 'owl':OWL, 'rdfs':RDFS, 'unit':UNIT}
        self.base_query = prefixes(initNs)
            
    def query(self, q, var):
        q = self.base_query + q
        if self.use_endpoint:
            return query_endpoint(self.endpoint, q, var)
        else:
            return query_graph(self.dataobject.graph, q)
        
    def query_type(self, t):
        """
        Return entities of type t.
        """
        q = """
            select ?s where {
                ?s rdf:type <%s>
            }
        """ % t
        return self.query(q, 's')
    
    def query_child(self, t):
        """
        Return children of t.
        """
        q = """
            select ?s where {
                ?s rdfs:subClassOf <%s> .
            }
        """ % t
        return self.query(q, 's')
    
    def query_label(self, t):
        """
        Return entities with label t.
        """
        q = """
            select ?s where {
                ?s rdfs:label "%s" . 
            }
        """ % t
        return self.query(q, 's')
    
    def query_parent(self, t):
        """
        Return parent of t.
        """
        q = """
            select ?s where {
                <%s> rdfs:subClassOf ?s .
            }
        """ % t
        return self.query(q, 's')
    
    def query_siblings(self, t, depth=1):
        if depth == -1: depth = '1,'
        q = """
            select ?s where {
                <%s> rdfs:subClassOf{%s} ?s .
            }
        """ % (t,str(depth))
        parents = self.query(q,'s')
        out = set()
        while parents:
            p = parents.pop(0)
            q = """
            select ?s where {
                ?s rdfs:subClassOf{%s} <%s> .
            }
            """ % (str(depth),t)
            out |= self.query(q,'s')
        return s
    
    def query_labels(self, t):
        q = """
            select ?s where {
                <%s> rdfs:label ?s .
            }
        """ % t
        return self.query(q, 's')
    
    def query_alt_labels(self, t):
        q = """
            select ?p ?s where {
                <%s> ?p ?s .
                ?p rdfs:subPropertyOf rdfs:label .
            } filter (isLiteral(?s))
        """ % t
        return self.query(q, ['p','s'])
    
    def construct_subgraph(self, t):
        """
        Return all triples *-connected to t.
        """
        out = set()
        tmp = set([t])
        visited = set()
        while tmp:
            curr = tmp.pop()
            visited.add(curr)
            q = """
                select ?s ?p ?o {
                    values ?s { <%s> }
                    ?s ?p ?o
                }
            """ % curr
            res = self.query(q, ['s','p','o'])
            out |= res
            tmp |= set([o for _,_,o in res])
            tmp -= visited
        
        return out
    
    
class TaxonomyAPI(API):
    def __init__(self, 
                 namespace='https://www.ncbi.nlm.nih.gov/taxonomy', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'NCBI API'):
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
    def query_species(self):
        return self.query_type(self.namespace['Taxon'])
    
    def query_division(self, t):
        if isinstance(t, (list,set)):
            return {a:self.query_division(a) for a in t}
        return self.query_subclassof(a)
    
    def query_ssd(self, t):
        if isinstance(t, (list,set)):
            return {a:self.query_subclassof(t) for a in t}
        return self.query_subclassof(t)
    
    def query_ranks(self):
        return self.query_type(self.namespace['Rank'])
    
    def query_rank(self, t):
        if isinstance(t, (list,set)):
            return {a:self.query_rank(t) for a in t}
        return self.query_subclassof(t)
    
class TraitsAPI(API):
    def __init__(self, namespace='https://eol.org/schema/terms/', endpoint = None, dataobject = None, name = 'EOL API'):
        super(TraitsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    def query_concervation_status(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_concervation_status(a) for a in t}
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/ontology/voc/SPMInfoItems#ConservationStatus> ?h .
            }
        """ % t
        return self.query(q,'h')

    def query_extinct_status(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_extinct_status(a) for a in t}
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/schema/terms/ExtinctionStatus> ?h .
            }
        """ % t
        return self.query(q,'h')
        
    def query_endemic_to(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_endemic_to(a) for a in t}
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/terms/endemic> ?h .
            }
        """ % t
        return self.query(q,'h')

    def query_ecoregion(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_ecoregion(a) for a in t}
        q =  """
            SELECT ?h WHERE {
                <%s> <https://www.wikidata.org/entity/Q295469> ?h .
            }
        """ % t
        return self.query(q,'h')
        
    def query_habitat(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_habitat(a) for a in t}
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/dwc/terms/habitat> ?h .
            }
        """ % t
        return self.query(q,'h')
    
class EcotoxChemicalAPI(API):
    def __init__(self, namespace='https://cfpub.epa.gov/ecotox/', endpoint = None, dataobject = None, name = 'ECOTOX Chemical API'):
        super(EcotoxChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        
    def query_chemical_names(self,t):
        if isinstance(t,(list,set)):
            return {a:self.query_chemical_names(a) for a in t}
        return self.query_labels(t)
    
    def query_chemicals(self):
        return self.query_type(self.namespace['Chemical'])
    
class EcotoxTaxonomyAPI(TaxonomyAPI):
    def __init__(self, namespace='https://cfpub.epa.gov/ecotox/', endpoint = None, dataobject = None, name = 'ECOTOX Taxonomy API'):
        super(EcotoxTaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    
class EffectsAPI(API):
    def __init__(self, namespace='https://cfpub.epa.gov/ecotox/', endpoint = None, dataobject = None, name = 'ECOTOX Effects API'):
        super(EffectsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    def query_chemicals_from_species(self,t):
        """
        Return chemical involved in experiment with species t.
        """
        if isinstance(t,(list,set)):
            return {a:self.query_chemicals_from_species(a) for a in t}
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species <%s> .
            ?t ns:chemical ?c .
            } 
        """ % t
        return self.query(q,'c')
    
    def query_species_from_chemicals(self,t):
        """
        Return species involved in experiment with chemical t.
        """
        if isinstance(t,(list,set)):
            return {a:self.query_species_from_chemicals(a) for a in t}
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            ?t ns:chemical <%s> .
            }
        """ % t
        return self.query(q,'c')
    
    def query_chemical(self):
        """
        Return chemicals used in at least one experiment.
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:chemical ?c .
            }
        """ % t
        return self.query(q,'c')
    
    def query_species(self):
        """
        Return species used in at least one experiment.
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            }
        """ % t
        return self.query(q,'c')
    
    def query_endpoint(self, c, s):
        if isinstance(c,(list,set)):
            return {a:self.query_endpoint(a,s) for a in c}
        if isinstance(s,(list,set)):
            return {a:self.query_endpoint(c,a) for a in s}
        q = """
            SELECT ?cc ?cu ?ep ?ef ?sd ?sdu {
                ?test rdf:type ns:Test ;
                  ns:chemical <%s> ;
                   ns:species <%s> ;
                   ns:hasResult [ 
                   ns:endpoint ?ep ;
                   ns:effect ?ef ;
                   ns:concentration [rdf:value ?cc ; 
                                        unit:units ?cu] ] .
               
                OPTIONAL {
                    ?test ns:studyDuration [rdf:value ?sd ;
                                            unit:units ?sdu] .
                }
            }""" % (str(c), str(s))
        
        return self.query(q, ['cc','cu','ep','ef','sd','sdu'])

        
class ChemicalAPI(API):
    def __init__(self, namespace='', endpoint = None, dataobject = None, name = 'Chemical API'):
        super(ChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        #TODO
    
    
    
    
    
    
    
    
    
    
