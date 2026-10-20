
## QC

Analyze clashes of protein structure given pdb
    - input pdb files are defined in SQL table chain_pdb.
    - create pipeline in pipelines/qc/chain_clashes.py. 
    - collect all methods from some scripts in src/structure_profiling.
      move them to the class Clashes defined structure/clashes.py.
    - clashes results are exported to SQL table chain_clashes.
      SQL statements are defined in src/sql/
    - add test cases in tests/

Analyze chain continuity and chain break of protein structure given pdb
    - input pdb files are defined in SQL table chain_pdb.
    - create pipeline in pipelines/qc/chain_breaks.py. 
    - collect all methods from some scripts in src/structure_profiling.
      move them to the class Breaks defined qc/breaks.py.
    - Analytic results are exported to SQL table chain_breaks.
      SQL statements are defined in src/sql/
    - add test cases in tests/

Analyze chain torsion and angles of protein structure given pdb
    - input pdb files are defined in SQL table chain_pdb.
    - create pipeline in pipelines/qc/chain_torsions.py. 
    - collect all methods from some scripts in src/structure_profiling.
      move them to the class Breaks defined qc/torsions.py.
    - Analytic results are exported to SQL table chain_torsions.
      SQL statements are defined in src/sql/
    - add test cases in tests/
    - add a command to bin/5_qc.sh


    Analyze pipelines/pdb_qc.py. identify any other methods that not covered by chain_clashes.py,
    or chain_breaks.py or chain_torsions.py. categorize them by rationals and 
    create one or multiple pipelines in pipelines/qc/.