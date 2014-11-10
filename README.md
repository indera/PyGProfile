#PyGProfile

### Installation

Clone this repository and then:

```bash
$ cd pygprofile/
$ python setup.py install --user
```

This will install the scripts in the pygprofile/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Summary

- pygp_go.py - TopGO, GOstats and goProfiles wrapper for running GO enrichment analysis from DESEQ2 results or custom input. Also a GO term finder included
- pygp_path.py - Wrapper for SPIA pathway enrichment analysis and also a simple hypergeometric test for pathway gene enrichment.
