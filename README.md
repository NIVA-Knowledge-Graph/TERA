# TERA
Code base for the publication [TERA: the Toxicological Effect and Risk Assessment Knowledge Graph](https://arxiv.org/abs/1908.10128). 

Documentation can be found [here](https://niva-knowledge-graph.github.io/TERA/).

TERA snapshot is avalible as a Zenodo dataset: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4244313.svg)](https://doi.org/10.5281/zenodo.4244313)


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

- Erik B. Myklebust, Ernesto Jimenez Ruiz, Jiaoyan Chen, Raoul Wolf, Knut Erik Tollefsen. **Prediction of Adverse Biological Effects of Chemicals Using Knowledge Graph Embeddings**. Accepted for publication in the Semantic Web Journal, 2021. ([arXiv](https://arxiv.org/abs/2112.04605)) ([Paper](http://semantic-web-journal.org/content/prediction-adverse-biological-effects-chemicals-using-knowledge-graph-embeddings-0)) ([REPOSITORY](https://github.com/NIVA-Knowledge-Graph/KGs_and_Effect_Prediction_2020))

-  Erik B. Myklebust, Ernesto Jimenez-Ruiz, Zofia C. Rudjord, Raoul Wolf, Knut Erik Tollefsen: [**Integrating  semantic  technologies  in  environmental  risk  assessment:  A  vision**](https://s3.amazonaws.com/setac.mms.uploads/m_48/extended_abstracts/49766_Myklebust/EBMyklebust_et_al_Semantics_and_risk_assessment.pdf).  29th Annual Meeting of the Society of Environmental Toxicology and Chemistry ([SETAC](https://helsinki.setac.org/)) (2019)

- Erik B. Myklebust, Ernesto Jimenez-Ruiz, Jiaoyan Chen, Raoul Wolf, Knut Erik Tollefsen. [**Knowledge Graph Embedding for Ecotoxicological Effect Prediction**](https://arxiv.org/abs/1907.01328). International Semantic Web Conference 2019. *Best Student Paper in the In-Use track of ISWC 2019*

- Erik B. Myklebust, Ernesto Jimenez-Ruiz, Jiaoyan Chen, Raoul Wolf, Knut Erik Tollefsen. [**TERA: the Toxicological Effect and Risk Assessment Knowledge Graph**](https://arxiv.org/abs/1908.10128). arXiv:1908.10128, 2019. 

