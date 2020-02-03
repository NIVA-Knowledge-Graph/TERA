
### tests

from DataAggregation import Taxonomy, PubChem, Effects, Traits, EcotoxTaxonomy, EcotoxChemicals
from DataIntegration import CasToInchikey, NCBIToEOL, InchikeyToPubChem
from DataAccess import TaxonomyAPI, TraitsAPI, EcotoxChemicalAPI, EcotoxTaxonomyAPI, EffectsAPI

directory = 'test_data/'

t = Taxonomy(directory=directory)
tr = Traits(directory=directory)
p = PubChem(directory=directory)
e = Effects(directory=directory)
et = EcotoxTaxonomy(directory=directory)
ec = EcotoxChemicals(directory=directory)

for g in [t,p,tr,e,et,ec]:
    print(g.name, '#triples:', len(g.graph))

tapi = TaxonomyAPI(dataobject=t)
trapi = TraitsAPI(dataobject=tr)
ecapi = EcotoxChemicalAPI(dataobject=ec)
etapi = EcotoxChemicalAPI(dataobject=et)
eapi = EffectsAPI(dataobject=e)

sp = tapi.query_species()
eol_sp = NCBIToEOL().convert(sp, strip=True)
eol_sp = [tr.namespace[e] for e in eol_sp if e != 'no mapping']

print(len(trapi.query_concervation_status(eol_sp)))





