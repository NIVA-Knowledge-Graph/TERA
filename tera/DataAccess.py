"""
A set of APIs to access data created with DataAggregation and DataIntegration moduels.
"""

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS, OWL
UNIT = Namespace('http://qudt.org/vocab/unit#')
from .utils import query_endpoint, test_endpoint, query_graph, prefixes, do_recursively_in_class, strip_namespace, tanimoto
from typing import Union

import pubchempy

from .DataIntegration import InchikeyToCas, InchikeyToChEBI, InchikeyToChEMBL, InchikeyToMeSH, InchikeyToPubChem, NCBIToEcotox, NCBIToEOL

class API:
    def __init__(self, 
                 namespace=None, 
                 endpoint=None, 
                 dataobject=None, 
                 name='API'):
        """API for accessing data sets. 
        
        Parameters
        ----------
        namespace : str, default None 
            Base URI for API
        
        endpoint : str, default None
            SPARQL endpoint URL

        dataobject : tera.DataObject, default None
            see DataAggregation
        
        Raises
        ------
        AssertionError 
            * If both endpoint and dataobject is None.
        
        """
        assert endpoint or dataobject
        
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
        """Pass SPARQL to graph or endpoint.
        
        Parameters
        ----------
        q : str
            sparql query
                
        var : str or list
            Bindings to return from query.
        
        Returns
        -------
        set
        """
        q = self.base_query + q
        if self.use_endpoint:
            return query_endpoint(self.endpoint, q, var)
        else:
            return query_graph(self.dataobject.graph, q)
        
    def query_type(self, t):
        """Return entities of type.
            
        Parameters
        ----------
        t : str or rdflib.URIRef
            Type URI.
    
        Returns 
        -------
        set 
        """
        q = """
            select ?s where {
                ?s rdf:type <%s>
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_child(self, t):
        """Return children.
        
        Parameters
        ----------
        t : str or rdflib.URIRef 
            Parent URI.
        
        Returns 
        -------
        set
        """
        q = """
            select ?s where {
                ?s rdfs:subClassOf <%s> .
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_label(self, t):
        """Return entities with label t.
        
        Parameters
        ----------
        t : str or rdflib.URIRef
        
        Returns 
        -------
        set 
        """
        q = """
            select ?s where {
                ?s rdfs:label "%s" . 
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_parent(self, t):
        """Return parent of t.
        
        Parameters
        ----------
            t : str or rdflib.URIRef
            
        Returns 
        -------
        set
        """
        q = """
            select ?s where {
                <%s> rdfs:subClassOf ?s .
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_siblings(self, t, depth=1):
        """Return (depth-1)-cousins.
        
        Parameters
        ----------
        t : str or rdflib.URIRef
        
        depth : int, default 1 
            Number of generation to search. 1 -> siblings, 1 -> 1st cousins, etc. 
        
        Returns 
        -------
        set
        """
        if depth == -1: depth = '1,'
        q = """
            select ?s where {
                <%s> rdfs:subClassOf{%s} ?s .
            }
        """ % (str(t),str(depth))
        parents = self.query(q,'s')
        out = set()
        while parents:
            p = parents.pop(0)
            q = """
            select ?s where {
                ?s rdfs:subClassOf{%s} <%s> .
            }
            """ % (str(depth),str(t))
            out |= self.query(q,'s')
        return s
    
    def query_alt_labels(self, t):
        """Get literals where prop =< rdfs:label.
        
        Parameters
        ----------
        t : str or rdflib.URIRef
        
        Returns 
        -------
        set
        """
        q = """
            select ?p ?s where {
                <%s> ?p ?s .
                ?p rdfs:subPropertyOf rdfs:label .
            } filter (isLiteral(?s))
        """ % str(t)
        return self.query(q, ['p','s'])
    
    def construct_subgraph(self, t):
        """Return all triples connected to input.
        
        Parameters
        ----------
        t : str or rdflib.URIRef
        
        Returns 
        -------
        set
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
            """ % str(curr)
            res = self.query(q, ['s','p','o'])
            out |= res
            tmp |= set([o for _,_,o in res])
            tmp -= visited
        
        return out
    
    @do_recursively_in_class
    def convert_id(self, id_: Union[URIRef, str, list, set],f,t, strip=False):
        """
        Convert between types of ids used in data.
        
        Parameters
        ----------
        f : str 
            input id type.
            
        t : str
            output id type.
            
        id_ : element or list 
            list of ids, 'no mapping' if no mapping between f and t exists.
            
        strip : bool  
            remove namespace from inputs
        
        Returns 
        -------
        str
            
        Raises
        ------
        NotImplementedError
            * If cannot convert between f and t.
        """

        if f == t: return id_
    
        if not hasattr(self, 'mappings'):
            raise AttributeError(self.name + ' has not attribute mappings.')
    
        if isinstance(id_, URIRef):
            id_ = str(id_)
        if strip:
            id_ = strip_namespace(id_, ['/','#','CID'])
        
        if f == self.base_identifier and t in self.mappings:
            return self.mappings[t].convert(id_)
        
        if f in self.mappings:
            return self.convert_id(self.mappings[f].convert(id_,reverse=True),
                                   f=self.base_identifier,t=t)
        
        raise NotImplementedError('From %s to %s is not supported.  Supported from/to values are %s', (f,t,','.join(self.mappings.keys())))
    
    def avalible_convertions(self):
        """Returns id types that can be converted between.
        
        Returns
        -------
        set
        """
        return set([self.base_identifier]) | set(self.mappings.keys())
    
class TaxonomyAPI(API):
    def __init__(self, 
                 namespace=None, 
                 endpoint=None,
                 dataobject = None,
                 mappings = {'eol',NCBIToEOL()},
                 name = 'Taxonomy API'):
        """Base class for accessing taxonomic data. 
        
        Parameters
        ----------
        dataobject : tera.Taxonomy 
            Data set to access using API. 
            
        mappings : dict 
            Mappings (tera.Alignment) from base_identifier (eg. ncbi) to other datasets. 
        """
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
        self.mappings = mappings
        self.base_identifier = 'ncbi'
        
    def get_taxa(self):
        """Return all taxa in taxonomy.
        
        Returns
        -------
        set
        """
        return self.query_type(self.namespace['Taxon'])
    
    @do_recursively_in_class
    def get_division(self, t: Union[URIRef, str, list, set]):
        """Return all taxa in division.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            Division URI 
        
        Returns 
        -------
        set
        """
        return self.query_subclassof(a)
    
    @do_recursively_in_class
    def get_ssd(self, t: Union[URIRef, str, list, set]):
        """Return all taxa in SSD.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set 
            SSD URI 
        
        Returns 
        -------
        set 
        """
        return self.query_subclassof(t)
    
    def get_ranks(self):
        """Return all ranks (taxonomic level).
        
        Returns 
        -------
        set
        """
        return self.query_type(self.namespace['Rank'])
    
    @do_recursively_in_class
    def get_rank(self, t: Union[URIRef, str, list, set]):
        """Return all taxa with rank.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            Rank URI
        
        Returns 
        -------
        set 
        """
        return self.query_subclassof(t)

class ChemicalAPI(API):
    def __init__(self, 
                 namespace=None, 
                 endpoint = None, 
                 dataobject = None, 
                 mappings = {'cas':InchikeyToCas(),
                            'cid':InchikeyToPubChem(),
                            'chebi':InchikeyToChEBI(),
                            'chemble':InchikeyToChEMBL(),
                            'mesh':InchikeyToMeSH()},
                 name = 'Chemical API'):
        """
        Base class for accessing chemical data. 
        
        Parameters
        ----------
        dataobject : tera.DataObject 
            Data set to access using API. 
        
        mappings : dict 
            Mappings from base_identifier (eg. ncbi) to other datasets. 
        """
        super(ChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        self.mappings = mappings
        self.base_identifier = 'inchikey'
   
    @do_recursively_in_class
    def get_fingerprint(self, id_: Union[URIRef, str, list, set], f='inchikey', strip = False):
        """Get binary fingerprints.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        f : str 
            Input identifier type. 
        
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        str
        """
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
        """Get synonyms.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        f : str 
            Input identifier type. 
            
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        list
        """
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
        """Return all triples connceted to input.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        f : str 
            Input identifier type. 
            
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        set 
        """
        a = self.convert_id(id_, f, 'cid', strip = strip)
        b = self.convert_id(id_, f, 'mesh', strip = strip)
        a = self.initNs['compound'][a]
        b = self.initNs['mesh'][b]
        return self.construct_subgraph(a) | self.construct_subgraph(b)
        
    @do_recursively_in_class
    def get_features(self, id_: Union[URIRef, str, list, set], params=None, f='inchikey', strip=False):
        """Return chemical features.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        params : list 
            Properties to return.
            eg. params = ['charge','molecular_weight','xlogp'] 
            To see all avalible features use which_features(). 
        
        f : str 
            Input identifier type. 
            
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        dict 
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
        """Chemical features avalible.
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        f : str 
            Input identifier type. 
            
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        list 
        """
        return [p for p in dir(Compound) if isinstance(getattr(Compound, p), property)]
    
    @do_recursively_in_class
    def simiarity(self, id_: Union[URIRef, str, list, set], ids, f='inchikey',strip=False):
        """Returns chemical simiarity between id and ids
        
        Parameters
        ----------
        id_ : rdflib.URIRef, str, list, or set 
            URI or identifier. 
        
        ids : list or set 
            URI or identifiers to compare against. 
        
        f : str 
            Input identifier type. 
            
        strip : bool 
            Remove namespace. Should be true if URI is passed.
        
        Returns 
        -------
        dict 
        """
        fp = self.get_fingerprint(id_, f, strip)
        fps = self.get_fingerprint(ids, f, strip)
        return {i:tanimoto(fp,f) for i,f in fps.items() if f and fp}
        
    def compounds(self):
        """Return all compounds.
        
        Returns 
        -------
        set 
        """
        q = """
            SELECT ?s {
            ?s  ?o  ?z
            FILTER (isURI(?s) && STRSTARTS(str(?s), str(compound:) ) )
            }
            """
        return self.query(q, var = 's')
    
class TraitsAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://eol.org/schema/terms/', 
                 endpoint = None, 
                 dataobject = None,
                 mappings = {'eol',NCBIToEOL()},
                 name = 'EOL API'):
        """
        Class for accessing EOL traits data.
        
        Parameters
        ----------
        namespace : str 
        
        endpoint : str 
        
        dataobject : tera.DataObject
        
        mapping : dict 
        """
        super(TraitsAPI, self).__init__(namespace, 
                                        endpoint, 
                                        dataobject, 
                                        mappings,
                                        name)
    
    @do_recursively_in_class
    def get_concervation_status(self,t: Union[URIRef, str, list, set]):
        """Return concervation status of t.
        
        Parameters
        ---------- 
        t : rdflib.URIRef, str, list, or set 
            URI 
        
        Returns 
        -------
        str 
        """
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/ontology/voc/SPMInfoItems#ConservationStatus> ?h .
            }
        """ % str(t)
        return self.query(q,'h')

    @do_recursively_in_class
    def get_extinct_status(self,t: Union[URIRef, str, list, set]):
        """Return extinct status (true/false).
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            URI
        
        Returns 
        -------
        str 
        """
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/schema/terms/ExtinctionStatus> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
        
    @do_recursively_in_class
    def get_endemic_to(self,t: Union[URIRef, str, list, set]):
        """Return endemic region.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            URI
        
        Returns 
        -------
        str 
        """
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/terms/endemic> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
    
    @do_recursively_in_class
    def get_ecoregion(self,t: Union[URIRef, str, list, set]):
        """Return ecoregion.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            URI
        
        Returns 
        -------
        str 
        """
        q =  """
            SELECT ?h WHERE {
                <%s> <https://www.wikidata.org/entity/Q295469> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
        
    @do_recursively_in_class
    def get_habitat(self,t: Union[URIRef, str, list, set]):
        """Return habiat.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, or set
            URI
        
        Returns 
        -------
        str 
        """
        q = """
            SELECT ?h WHERE {
                <%s> <http://rs.tdwg.org/dwc/terms/habitat> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
    
class EcotoxChemicalAPI(ChemicalAPI):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Chemical API'):
        """
        Class for accessing chemical data in Ecotox. 
        
        Parameters
        ----------
        namespace : str 
        
        endpoint : str 
        
        dataobject : tera.DataObject
        """
        super(EcotoxChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        
    @do_recursively_in_class
    def query_chemical_names(self,t: Union[URIRef, str, list, set]):
        """
        Return chemical names.
        
        Parameters
        ---------- 
        t : rdflib.URIRef, str, list, set 
            URI 
        
        Returns 
        -------
        str
        """
        return self.query_labels(t)
    
    def query_chemicals(self):
        """Return set of all chemicals.
        
        Returns
        -------
        set
        """
        return self.query_type(self.namespace['Chemical'])
    
class EcotoxTaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Taxonomy API'):
        """Class for accessing Ecotox taxonomic data.
        
        Parameters
        ----------
        namespace : str 
        
        endpoint : str 
        
        dataobject : tera.DataObject
        """
        super(EcotoxTaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
class NCBITaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://www.ncbi.nlm.nih.gov/taxonomy', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'NCBI API'):
        """Class for accessing NCBI taxonomic data.
        
        Parameters
        ----------
        namespace : str 
        
        endpoint : str 
        
        dataobject : tera.DataObject
        """
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
    
class EffectsAPI(API):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Effects API'):
        """Class for accessing Ecotox effect data.
        
        Parameters
        ----------
        namespace : str 
        
        endpoint : str 
        
        dataobject : tera.DataObject
        """
        super(EffectsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    @do_recursively_in_class
    def get_chemicals_from_species(self,t: Union[URIRef, str, list, set]):
        """Return chemical involved in experiment with certain species.
  
        Parameters
        ----------
        t : rdflib.URIRef, str, list, set 
            Species URI 
        
        Returns 
        -------
        set 
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species <%s> .
            ?t ns:chemical ?c .
            } 
        """ % str(t)
        return self.query(q,'c')
    
    @do_recursively_in_class
    def get_species_from_chemicals(self, t: Union[URIRef, str, list, set]):
        """Return species involved in experiment using chemical.
        
        Parameters
        ----------
        t : rdflib.URIRef, str, list, set 
            Chemical URI 
        
        Returns 
        -------
        set 
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            ?t ns:chemical <%s> .
            }
        """ % str(t)
        return self.query(q,'c')
    
    def get_chemicals(self):
        """Return chemicals used in at least one experiment.
        
        Returns 
        -------
        set 
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:chemical ?c .
            }
        """
        return self.query(q,'c')
    
    def get_species(self):
        """Return species used in at least one experiment.
        
        Returns 
        -------
        set 
        """
        q = """
        select ?c where {
            ?t rdf:type ns:Test .
            ?t ns:species ?c .
            }
        """
        return self.query(q,'c')
    
    def get_endpoint(self, 
                       c: Union[URIRef, str, list, set], 
                       s: Union[URIRef, str, list, set]):
        """
        Return endpoints that use chemical c and species s.
        
        Parameters
        ---------- 
        c : rdflib.URIRef, str, list, set 
            Chemical URIs. If None, c <- query_chemicals
        s : rdflib.URIRef, str, list, set 
            Species URIs. If None, s <- query_species
        
        Returns 
        -------
        dict 
            On the form {chemical:{species:endpoint values}}
        """
        if not c:
            c = self.query_chemicals()
        if not s:
            s = self.query_species()
            
        if isinstance(c,(list,set)):
            return {a:self.get_endpoint(a,s) for a in c}
        if isinstance(s,(list,set)):
            return {a:self.get_endpoint(c,a) for a in s}
        q = """
            SELECT ?cc ?cu ?ep ?ef ?sd ?sdu WHERE {
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
