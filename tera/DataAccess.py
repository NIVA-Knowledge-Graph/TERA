
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
        """
        Args:
        namespace : str \n
            Base URI for API. \n
        endpoint : str \n 
            SPARQL endpoint url \n 
        dataobject : tera.DataObject \n  
            see DataAggregation \n 
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
        """Pass SPARQL to graph or endpoint.,
        
        Args:
            q : str \n
                sparql query \n
            var : str or list \n
                Bindings to return from query.
        Return:
            Query result.
        """
        q = self.base_query + q
        if self.use_endpoint:
            return query_endpoint(self.endpoint, q, var)
        else:
            return query_graph(self.dataobject.graph, q)
        
    def query_type(self, t):
        """Return entities of type.
        Args:
            t : str or rdflib.URIRef \n
            Type URI.
        Return:
            set \n
            Entities of type t. 
        """
        q = """
            select ?s where {
                ?s rdf:type <%s>
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_child(self, t):
        """Return children.
        Args:
            t : str or rdflib.URIRef \n
            Parent URI.
        Return:
            set \n
            Entities in t. 
        """
        q = """
            select ?s where {
                ?s rdfs:subClassOf <%s> .
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_label(self, t):
        """Return entities with label t.
        Args:
            t : str or rdflib.URIRef \n
            URI.
        Return:
            set \n
            Labels of t. 
        """
        q = """
            select ?s where {
                ?s rdfs:label "%s" . 
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_parent(self, t):
        """Return parent of t.
        Args:
            t : str or rdflib.URIRef \n
            Child URI.
        Return:
            set \n
            Parents of t.
        """
        q = """
            select ?s where {
                <%s> rdfs:subClassOf ?s .
            }
        """ % str(t)
        return self.query(q, 's')
    
    def query_siblings(self, t, depth=1):
        """Return (depth-1)-cousins.
        Args:
            t : str or rdflib.URIRef \n
                URI. \n 
            depth : int \n 
                Number of generation to search. 1 -> siblings, 1 -> 1st cousins, etc. 
        Return:
            set \n
            Cousins of t.
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
        Args:
            t : str or rdflib.URIRef \n
                URI. \n 
        Return:
            set \n
            Alt labels of t.
        """
        q = """
            select ?p ?s where {
                <%s> ?p ?s .
                ?p rdfs:subPropertyOf rdfs:label .
            } filter (isLiteral(?s))
        """ % str(t)
        return self.query(q, ['p','s'])
    
    def construct_subgraph(self, t):
        """
        Return all triples connected to input.
        Args:
            t : str or rdflib.URIRef \n
                URI
        Return:
            set \n
            Set of triples connected by arbitrary links to input. 
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
        Args:
            f : str \n 
                input id type. \n
            f : str  \n
                output id type.\n
            id_ : element or list \n 
                list of ids, 'no mapping' if no mapping between f and t exists. \n
            strip : bool \n  
                remove namespace from inputs
        Return:
            str or dict \n
            If single input, returns string. Else, dict on the form {from:to}. 
        Raises:
            NotImplementedError: if cannot convert between f and t.
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
        
        raise NotImplementedError('From %s to %s is not supported. \n Supported from/to values are %s', (f,t,','.join(self.mappings.keys())))
    
    def avalible_convertions(self):
        """Returns id types that can be converted between."""
        return set([self.base_identifier]) | set(self.mappings.keys())
    
class TaxonomyAPI(API):
    def __init__(self, 
                 namespace=None, 
                 endpoint=None,
                 dataobject = None,
                 mappings = {'eol',NCBIToEOL()},
                 name = 'Taxonomy API'):
        """
        Base class for accessing taxonomic data. 
        Args:
            dataobject : tera.Taxonomy \n
                Data set to access using API. \n
            mappings : dict \n
                Mappings from base_identifier (eg. ncbi) to other datasets. 
        """
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
        self.mappings = mappings
        self.base_identifier = 'ncbi'
        
    def get_taxa(self):
        """Return all taxa in taxonomy."""
        return self.query_type(self.namespace['Taxon'])
    
    @do_recursively_in_class
    def get_division(self, t: Union[URIRef, str, list, set]):
        """Return all taxa in division.
        Args:
            t : rdflib.URIRef, str, list, or set \n
            Division URI \n
        Return:
            set \n
            Set of taxa in devision.
        """
        return self.query_subclassof(a)
    
    @do_recursively_in_class
    def get_ssd(self, t: Union[URIRef, str, list, set]):
        """Return all taxa in SSD.
        Args:
            t : rdflib.URIRef, str, list, or set \n
            SSD URI \n
        Return:
            set \n
            Set of taxa in SSD.
        """
        return self.query_subclassof(t)
    
    def get_ranks(self):
        """Return all ranks (taxonomic level).
        Return:
            set \n
            Set of ranks.
        """
        return self.query_type(self.namespace['Rank'])
    
    @do_recursively_in_class
    def get_rank(self, t: Union[URIRef, str, list, set]):
        """Return all taxa with rank.
        Args:
            t : rdflib.URIRef, str, list, or set \n
            Rank URI \n
        Return:
            set \n
            Set of taxa at rank.
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
        Args:
            dataobject : tera.DataObject \n
                Data set to access using API. \n
            mappings : dict \n
                Mappings from base_identifier (eg. ncbi) to other datasets. 
        """
        super(ChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        self.mappings = mappings
        self.base_identifier = 'inchikey'
   
    @do_recursively_in_class
    def get_fingerprint(self, id_: Union[URIRef, str, list, set], f='inchikey', strip = False):
        """Get binary fingerprints.
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            fp : str \n
            Binary fingerprint.
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
        """Get synonyms for id_.
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            list \n
            Synonyms.
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
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            set \n
            Triples.
        """
        a = self.convert_id(id_, f, 'cid', strip = strip)
        b = self.convert_id(id_, f, 'mesh', strip = strip)
        a = self.initNs['compound'][a]
        b = self.initNs['mesh'][b]
        return self.construct_subgraph(a) | self.construct_subgraph(b)
        
    @do_recursively_in_class
    def get_features(self, id_: Union[URIRef, str, list, set], params=None, f='inchikey', strip=False):
        """Return chemical features.
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            params : list \n
                Properties to return \n 
                eg. params = ['charge','molecular_weight','xlogp'] \n
                To see all avalible features use which_features(). \n
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            dict \n
            Chemical features.
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
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            list \n
            Chemical features.
        """
        return [p for p in dir(Compound) if isinstance(getattr(Compound, p), property)]
    
    @do_recursively_in_class
    def simiarity(self, id_: Union[URIRef, str, list, set], ids, f='inchikey',strip=False):
        """Returns chemical simiarity between id and ids
        Args:
            id_ : rdflib.URIRef, str, list, or set \n
                URI or identifier. \n
            ids : list or set \n
                URI or identifiers to compare against. 
            f : str \n
                Input identifier type. \n
            strip : bool \n
                Remove namespace. Should be true if URI is passed.
        Return:
            dict \n
            Chemical simiarity.
        """
        fp = self.get_fingerprint(id_, f, strip)
        fps = self.get_fingerprint(ids, f, strip)
        return {i:tanimoto(fp,f) for i,f in fps.items() if f and fp}
        
    def compounds(self):
        """Return all compounds.
        Return:
            set \n
            All compounds.
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
        """
        super(TraitsAPI, self).__init__(namespace, 
                                        endpoint, 
                                        dataobject, 
                                        mappings,
                                        name)
    
    @do_recursively_in_class
    def get_concervation_status(self,t: Union[URIRef, str, list, set]):
        """Return concervation status of t.
        Args: 
            t : rdflib.URIRef, str, list, or set \n
            URI 
        Return:
            concervation status
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
        Args:
            t : rdflib.URIRef, str, list, or set\n
                URI
        Return:
            str \n
                extinct status
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
        Args:
            t : rdflib.URIRef, str, list, or set\n
                URI
        Return:
            str \n
                endemic to"""
        q = """
            SELECT ?h WHERE {
                <%s> <http://eol.org/terms/endemic> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
    
    @do_recursively_in_class
    def get_ecoregion(self,t: Union[URIRef, str, list, set]):
        """Return ecoregion.
        Args:
            t : rdflib.URIRef, str, list, or set\n
                URI
        Return:
            str \n
                ecoregion"""
        q =  """
            SELECT ?h WHERE {
                <%s> <https://www.wikidata.org/entity/Q295469> ?h .
            }
        """ % str(t)
        return self.query(q,'h')
        
    @do_recursively_in_class
    def get_habitat(self,t: Union[URIRef, str, list, set]):
        """Return habiat.
        Args:
            t : rdflib.URIRef, str, list, or set\n
                URI
        Return:
            str \n
                habiat"""
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
        """
        super(EcotoxChemicalAPI, self).__init__(namespace, endpoint, dataobject, name)
        
    @do_recursively_in_class
    def query_chemical_names(self,t: Union[URIRef, str, list, set]):
        """
        Return chemical names.
        Args: 
            t : rdflib.URIRef, str, list, set \n
            URI 
        Return:
            str\n
            Chemical label.
        """
        return self.query_labels(t)
    
    def query_chemicals(self):
        """Return set of all chemicals."""
        return self.query_type(self.namespace['Chemical'])
    
class EcotoxTaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Taxonomy API'):
        """Class for accessing Ecotox taxonomic data."""
        super(EcotoxTaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
        
class NCBITaxonomyAPI(TaxonomyAPI):
    def __init__(self, 
                 namespace='https://www.ncbi.nlm.nih.gov/taxonomy', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'NCBI API'):
        """Class for accessing NCBI taxonomic data."""
        super(TaxonomyAPI, self).__init__(namespace, endpoint, dataobject, name)
    
class EffectsAPI(API):
    def __init__(self, 
                 namespace='https://cfpub.epa.gov/ecotox/', 
                 endpoint = None, 
                 dataobject = None, 
                 name = 'ECOTOX Effects API'):
        """Class for accessing Ecotox effect data."""
        super(EffectsAPI, self).__init__(namespace, endpoint, dataobject, name)
    
    @do_recursively_in_class
    def get_chemicals_from_species(self,t: Union[URIRef, str, list, set]):
        """Return chemical involved in experiment with certain species.
        Args:
            t : rdflib.URIRef, str, list, set \n
                Species URI 
        Return:
            set \n
                Chemical URIs
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
        Args:
            t : rdflib.URIRef, str, list, set \n
                Chemical URI 
        Return:
            set \n
                Species URIs
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
        Return:
            set \n
            Chemical URIs
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
        Return:
            set \n
            Species URIs"""
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
        Args: 
            c : rdflib.URIRef, str, list, set \n
                Chemical URIs. If None, c <- query_chemicals\n
            s : rdflib.URIRef, str, list, set \n
                Species URIs. If None, s <- query_species\n
        Return:
            dict \n
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
