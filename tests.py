
### tests

from DataAggregation import Taxonomy, PubChem, Effects, Traits, EcotoxTaxonomy, EcotoxChemicals
from DataIntegration import CasToInchikey, NCBIToEOL, InchikeyToPubChem

directory = 'test_data/'

t = Taxonomy(directory=directory)
tr = Traits(directory=directory)
p = PubChem(directory=directory)
e = Effects(directory=directory)
et = EcotoxTaxonomy(directory=directory)
ec = EcotoxChemicals(directory=directory)

for g in [t,p,tr,e,et,ec]:
    print(g.name, '#triples:', len(g.graph))
    


