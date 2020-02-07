
### tests

from tera.DataAggregation import Taxonomy, PubChem, Effects, Traits, EcotoxTaxonomy, EcotoxChemicals
from tera.DataAccess import TaxonomyAPI, TraitsAPI, EcotoxChemicalAPI, EcotoxTaxonomyAPI, EffectsAPI, ChemicalAPI
from tera.DataIntegration import NCBIToEOL, NCBIToEcotox

directory = 'test_data/'

t = Taxonomy(directory=directory)
tr = Traits(directory=directory)
p = PubChem(directory=directory)
e = Effects(directory=directory)
et = EcotoxTaxonomy(directory=directory)
ec = EcotoxChemicals(directory=directory)

for g in [t,p,tr,e,et,ec]:
    print(g.name, '#triples:', len(g.graph))

eapi = EffectsAPI(dataobject=e)
tapi = TaxonomyAPI(dataobject=t, mappings={'eol':NCBIToEOL,
                                           'ecotox':NCBIToEcotox(t,et)})
capi = ChemicalAPI(dataobject=p)

species = list(eapi.get_species())
chemicals = list(eapi.get_chemicals())

endpoints = eapi.get_endpoint(chemicals[:10],species[:10])

for k1 in endpoints:
    for k2 in endpoints[k1]:
        tmp = endpoints[k1][k2]
        if tmp: 
            print(capi.convert_id(k1,f='cas',t='inchikey',strip=True),
              tapi.convert_id(k2,f='ecotox',t='ncbi',strip=False))
            print(tmp)
