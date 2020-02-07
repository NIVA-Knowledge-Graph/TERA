# TERA
Code base for the publication [TERA: the Toxicological Effect and Risk Assessment Knowledge Graph](https://arxiv.org/abs/1908.10128). 

Documentation can be found [here](https://erik-bm.github.io/TERA/).

## Data sources
* [ECOTOX](https://cfpub.epa.gov/ecotox/)
* [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
* [PubChem](https://pubchem.ncbi.nlm.nih.gov)
* [ChEMBL](https://www.ebi.ac.uk/chembl/)
* Encyclopedia of Life ([EOL](https://www.eol.org))
* [Wikidata](https://www.wikidata.org/wiki/Wikidata:Main_Page)
* Other data sources can be easily added.

The raw data can be downloaded by running
```
bash download_data.sh
```

## Install
```
pip3 install -r requirements.txt
python3 setup.py install
```

## Examples
See `tests.py`. More to come...


## Related publications

-  Erik B. Myklebust, Ernesto Jimenez-Ruiz, Zofia C. Rudjord, Raoul Wolf, Knut Erik Tollefsen: [**Integrating  semantic  technologies  in  environmental  risk  assessment:  A  vision**](https://s3.amazonaws.com/setac.mms.uploads/m_48/extended_abstracts/49766_Myklebust/EBMyklebust_et_al_Semantics_and_risk_assessment.pdf).  29th Annual Meeting of the Society of Environmental Toxicology and Chemistry ([SETAC](https://helsinki.setac.org/)) (2019)

- Erik B. Myklebust, Ernesto Jimenez-Ruiz, Jiaoyan Chen, Raoul Wolf, Knut Erik Tollefsen. [**Knowledge Graph Embedding for Ecotoxicological Effect Prediction**](https://arxiv.org/abs/1907.01328). International Semantic Web Conference 2019. *Best Student Paper in the In-Use track of ISWC 2019*

- Erik B. Myklebust, Ernesto Jimenez-Ruiz, Jiaoyan Chen, Raoul Wolf, Knut Erik Tollefsen. [**TERA: the Toxicological Effect and Risk Assessment Knowledge Graph**](https://arxiv.org/abs/1908.10128). Submitted to a conference, 2019. 

