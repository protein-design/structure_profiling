
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