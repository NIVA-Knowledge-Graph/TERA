
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS, OWL
UNIT = Namespace('http://qudt.org/vocab/unit#')
from .utils import query_endpoint, test_endpoint, query_graph, prefixes, do_recursively_in_class, strip_namespace, tanimoto
from typing import Union

import pubchempy

from .DataIntegration import InchikeyToCas, InchikeyToChEBI, InchikeyToChEMBL, InchikeyToMeSH, InchikeyToPubChem, NCBIToEcotox, NCBIToEOL

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
            
        self.initNs = {'rdf':RDF, 
                       'ns':self.namespace, 
                       'owl':OWL, 
                       'rdfs':RDFS, 
                       'unit':UNIT,
                       'mesh':Namespace('http://id.nlm.nih.gov/mesh/'),
                       'obo':Namespace('http://purl.obolibrary.org/obo/'),
                       'pubchem':Namespace('http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#'),
                       'compound':Namespace('http://rdf.ncbi.nlm.nih.gov/pubchem/compound/')}
        self.base_query = prefixes(self.initNs)
            
    def query(self, q, var):
        """Pass SPARQL to graph or endpoint"""
        q = self.base_query + q
        if self.use_endpoint:
            return query_endpoint(self.endpoint, q, var)
        else:
            return query_graph(self.dataobject.graph, q)
        
    def query_type(self, t):
        """Return entities of type t."""
        q = """
            select ?s where {
                ?s rdf:type <%s>
            }
        """ % t
        return self.query(q, 's')
    
    def query_child(self, t):
        """Return children of t."""
        q = """
            select ?s where {
                ?s rdfs:subClassOf <%s> .
            }
        """ % t
        return self.query(q, 's')
    
    def query_label(self, t):
        """Return entities with label t."""
        q = """
            select ?s where {
                ?s rdfs:label "%s" . 
            }
        """ % t
        return self.query(q, 's')
    
    def query_parent(self, t):
        """Return parent of t."""
        q = """
            select ?s where {
                <%s> rdfs:subClassOf ?s .
            }
        """ % t
        return self.query(q, 's')
    
    def query_siblings(self, t, depth=1):
        """Return (depth-1)-cusins of t."""
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
        """Get rdfs:labels of t"""
        q = """
            select ?s where {
                <%s> rdfs:label ?s .
            }
        """ % t
        return self.query(q, 's')
    
    def query_alt_labels(self, t):
        """Get literals of t where p =< rdfs:label."""
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
    
    @do_recursively_in_class
    def convert_id(self, id_: Union[URIRef, str, list, set],f,t, strip=False):
        """
        Convert between types of ids used in data.
        
        Args:
            from_ :: str 
                input id type, inchikey, cas, chemical_id or cid.
            to_ :: str 
                output id type, inchikey, cas, chemical_id or cid.
            id_ :: list or element
                list of ids, 'no mapping' if no mapping between from_ and to_ exists.
            strip :: bool
                remove namespace from inputs
        Return:
            converted ids.
        Raises:
            NotImplementedError: if from_ or to_ not in (inchikey, cas, cid, mesh, chebi, chemble)
        """

        if f == t: return id_
    
        if not hasattr(self, 'mappings'):
            raise AttributeError(self.name + ' has not attribute mappings.')
    
        if isinstance(id_, (list, set)):
            return {i:self.convert_id(i) for i in id_}
        if isinstance(id_, URIRef):
            id_ = str(id_)
        if strip:
            id_ = strip_namespace(id_, ['/','#','CID'])
        
        if f == self.base_identifier and t in self.mappings:
            return self.mapping[t].convert(id_)
        
        if f in self.mapping:
            return self.convert_id(self.mappings[f].convert(id_,reverse=True),
                                   f=self.base_identifier,t=t)
        
        raise NotImplementedError('From %s to %s is not supported. \n Supported from/to values are %s', (f,t,','.join(self.mappings.keys())))
    
class TaxonomyAPI(API):
    def __init__(self, 
                 namespace=None, 
                 endpoint = None,
                 dataobject = None,
                 mapping_dataobject = None,
                 name = 'Taxonomy API'):
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
        self.mappings = {'eol':NCBIToEOL()}
        if mapping_dataobject:
            self.mappings['ecotox'] = NCBIToEcotox(dataobject, mapping_dataobject)
        elif test_endpoint:
            self.mappings['ecotox'] = EndpointMapping(endpoint)
        
        self.base_identifier = 'ncbi'
        
    def get_species(self):
        """All species"""
        return self.query_type(self.namespace['Taxon'])
    
    @do_recursively_in_class
    def get_division(self, t: Union[URIRef, str, list, set]):
        """Return all species in division t"""
        return self.query_subclassof(a)
    
    @do_recursively_in_class
    def get_ssd(self, t: Union[URIRef, str, list, set]):
        """Return all species in ssd t"""
        return self.query_subclassof(t)
    
    def get_ranks(self):
        """Get all taxonomic ranks"""
        return self.query_type(self.namespace['Rank'])
    
    @do_recursively_in_class
    def get_rank(self, t: Union[URIRef, str, list, set]):
        """Return all taxons with rank t"""
        return self.query_subclassof(t)

class ChemicalAPI(API):
    def __init__(self, namespace='', endpoint = None, dataobject = None, name = 'Chemical API'):
        super(ChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        
        self.mappings = {'cas':InchikeyToCas(),
                         'cid':InchikeyToPubChem(),
                         'chebi':InchikeyToChEBI(),
                         'chemble':InchikeyToChEMBL(),
                         'mesh':InchikeyToMeSH()}
        self.base_identifier = 'inchikey'
   
    @do_recursively_in_class
    def get_fingerprint(self, id_: Union[URIRef, str, list, set], f='inchikey', strip = False):
        """Get bin fingerprints for id_."""
        c = self.convert_id(id_, f, 'cid', strip)
        
        fp = None
        try:
            fp = Compound.from_cid(c).fingerprint
            fp = bin(int(fp, 16))
        except pubchempy.BadRequestError as e:
            print(c,e)
        except pubchempy.NotFoundError as e:
            print(c,e)

        return fp
    
    @do_recursively_in_class
    def get_names(self, id_: Union[URIRef, str, list, set], f='inchikey', strip = False):
        """Get synonyms for id_."""
        c = self.convert_id(id_, f, 'cid', strip)
        out = []
        try:
            out = Compound.from_cid(c).synonyms
        except pubchempy.BadRequestError as e:
            print(c,e)
        except pubchempy.NotFoundError as e:
            print(c,e)

        return out
    
    @do_recursively_in_class
    def class_hierarchy(self, id_: Union[URIRef, str, list, set], f='inchikey', strip=False):
        """Return all triples connceted to id_."""
        a = self.convert_id(id_, f, 'cid', strip = strip)
        b = self.convert_id(id_, f, 'mesh', strip = strip)
        a = self.initNs['compound'][a]
        b = self.initNs['mesh'][b]
        return self.construct_subgraph(a) | self.construct_subgraph(b)
        
    @do_recursively_in_class
    def get_features(self, id_: Union[URIRef, str, list, set], params=None, f='inchikey', strip=False):
        """Return chemical features.
            params :: list 
                properties to return
                ex: params = ['charge','molecular_weight','xlogp']
                To see all avalible features use which_features() . 
        """
        id_ = self.convert_id(id_, f, 'cid', strip=strip)
        id_ = self.initNs['compound'][id_]
        
        out = dict()
        try:
            if params:
                out = Compound.from_cid(c).to_dict(properties = params)
            else:
                out = Compound.from_cid(c).to_dict()
        except pubchempy.NotFoundError as e:
            print(c,e)
        except pubchempy.BadRequestError as e:
            print(c,e)

        return out
    
    @do_recursively_in_class
    def which_features(self, id_: Union[URIRef, str, list, set], f='inchikey', strip=False):
        """Chemical features avalible."""
        return self.get_features(id_, f, strip).keys()
    
    @do_recursively_in_class
    def simiarity(self, id_: Union[URIRef, str, list, set], ids, f='inchikey',strip=False):
        """Returns chemical simiarity between id and ids"""
        fp = self.get_fingerprint(id_, f, strip)
        fps = self.get_fingerprint(ids, f, strip)
        return {i:tanimoto(fp,f) for i,f in zip(ids,fps) if f and fp}
        
    def compounds(self):
        """Return all compounds."""
        q = """
            SELECT ?s {
            ?s  ?o  ?z
            FILTER (isURI(?s) && STRSTARTS(str(?s), str(compound:) ) )
            }
            """
        return self.query(q, var = 's')
    
class TraitsAPI(TaxonomyAPI):
    def __init__(self, namespace='https://eol.org/schema/terms/', endpoint = None, dataobject = None, name = 'EOL API'):
        super(TraitsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    @do_recursively_in_class
    def get_concervation_status(self,t: Union[URIRef, str, list, set]):
        """Return concervation status of t."""
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/ontology/voc/SPMInfoItems#ConservationStatus> ?h .
            }
        """ % t
        return self.query(q,'h')

    @do_recursively_in_class
    def get_extinct_status(self,t: Union[URIRef, str, list, set]):
        """Return extinct status (true/false) of t. """
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/schema/terms/ExtinctionStatus> ?h .
            }
        """ % t
        return self.query(q,'h')
        
    @do_recursively_in_class
    def get_endemic_to(self,t: Union[URIRef, str, list, set]):
        """Return endemic region of t."""
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/terms/endemic> ?h .
            }
        """ % t
        return self.query(q,'h')
    
    @do_recursively_in_class
    def get_ecoregion(self,t: Union[URIRef, str, list, set]):
        """Return ecoregion of t. """
        q =  """
            SELECT ?h WHERE {
                <%s> <https://www.wikidata.org/entity/Q295469> ?h .
            }
        """ % t
        return self.query(q,'h')
        
    @do_recursively_in_class
    def get_habitat(self,t: Union[URIRef, str, list, set]):
        """Return habiat of t."""
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/dwc/terms/habitat> ?h .
            }
        """ % t
        return self.query(q,'h')
    
class EcotoxChemicalAPI(ChemicalAPI):
    def __init__(self, namespace='https://cfpub.epa.gov/ecotox/', endpoint = None, dataobject = None, name = 'ECOTOX Chemical API'):
        super(EcotoxChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        
    @do_recursively_in_class
    def query_chemical_names(self,t: Union[URIRef, str, list, set]):
        return self.query_labels(t)
    
    def query_chemicals(self):
        return self.query_type(self.namespace['Chemical'])
    
class EcotoxTaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Taxonomy API'):
        super(EcotoxTaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
class NCBITaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://www.ncbi.nlm.nih.gov/taxonomy', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'NCBI API'):
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
    
class EffectsAPI(API):
    def __init__(self, namespace='https://cfpub.epa.gov/ecotox/', endpoint = None, dataobject = None, name = 'ECOTOX Effects API'):
        super(EffectsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    @do_recursively_in_class
    def query_chemicals_from_species(self,t: Union[URIRef, str, list, set]):
        """Return chemical involved in experiment with species t."""
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species <%s> .
            ?t ns:chemical ?c .
            } 
        """ % t
        return self.query(q,'c')
    
    @do_recursively_in_class
    def query_species_from_chemicals(self,t: Union[URIRef, str, list, set]):
        """Return species involved in experiment with chemical t."""
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            ?t ns:chemical <%s> .
            }
        """ % t
        return self.query(q,'c')
    
    def query_chemicals(self):
        """Return chemicals used in at least one experiment."""
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:chemical ?c .
            }
        """ % t
        return self.query(q,'c')
    
    def query_species(self):
        """Return species used in at least one experiment."""
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            }
        """ % t
        return self.query(q,'c')
    
    def query_endpoint(self, 
                       c: Union[URIRef, str, list, set], 
                       s: Union[URIRef, str, list, set]):
        """
        Return endpoints that use chemical c and species s.
        If list is none: c = query_chemicals and s = query_species
        """
        if not c:
            c = self.query_chemicals()
        if not s:
            s = self.query_species()
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


