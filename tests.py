
### tests

from tera import Taxonomy, PubChem, Effects, Traits, EcotoxTaxonomy, EcotoxChemicals
from tera import CasToInchikey, NCBIToEOL, InchikeyToPubChem
from tera import TaxonomyAPI, TraitsAPI, EcotoxChemicalAPI, EcotoxTaxonomyAPI, EffectsAPI

directory = 'test_data/'

t = Taxonomy(directory=directory)
tr = Traits(directory=directory)
p = PubChem(directory=directory)
e = Effects(directory=directory)
et = EcotoxTaxonomy(directory=directory)
ec = EcotoxChemicals(directory=directory)

for g in [t,p,tr,e,et,ec]:
    print(g.name, '#triples:', len(g.graph))

tapi = TaxonomyAPI(dataobject=t, mapping_dataobject=et)

help(tapi)

