'''
Chain-continuity / chain-break analysis for one Bio.PDB chain.

Collected here from src/structure_profiling (previously spread over
geometry/geometry_distance.py, geometry/geometry.py and qc/qc_screen.py):
- missing_residues:  residues lacking backbone atoms (N/CA/C/O)
- ca_distances:      consecutive CA(i)-CA(i+1) distances (ideal ~3.8 A)
- peptide_distances: consecutive peptide C(i)-N(i+1) distances
                     (ideal 1.329 A)
- chain_continuity:  localise real discontinuities ('breaks') and
                     stretched-but-intact peptides ('strained') per
                     consecutive residue pair

run() returns the summary exported by pipelines/qc/chain_breaks.py to
the SQL table chain_breaks; check() gates it on the Breaks cutoffs
(the BACKBONE_GATES values of the former qc_screen screens).

Cutoffs:
    cn_break = 2.0 A   C-N above this is a real break (ideal 1.329 A)
    ca_break = 4.5 A   CA-CA above this is a break (ideal ~3.8 A)
    cn_strain = 1.5 A  C-N above this but <= cn_break is a strained peptide
A pair is also a break when C-N is unmeasurable (missing atoms) or the
residue numbering skips (resnum_gap != 1).

references:
wwPDB OneDep validation reports (X-ray + EM guides, updated 2025).
'''
import numpy as np
import pandas as pd


class Breaks:
    '''chain-continuity analysis of one Bio.PDB chain (breaks / fragments).'''

    # continuity cutoffs (A)
    cn_break:float = 2.0       # C-N above this is a real break
    ca_break:float = 4.5       # CA-CA above this is a break
    cn_strain:float = 1.5      # C-N above this but <= cn_break is strained

    # gates (fail -> reason), the former BACKBONE_GATES of qc_screen.py
    GATES = {
        'num_missing': 0,      # no residue may lack N/CA/C/O atoms
        'num_breaks': 0,       # real discontinuities fail
        'num_fragments': 1,    # one contiguous chain
        'num_strained': 3,     # design tolerance (tune)
        'ca_max_dist': 5.0,    # ideal CA-CA ~3.8 A
        'pep_max_dist': 5.0,   # ideal peptide C-N 1.329 A
    }

    def __init__(self, chain):
        self.chain = chain

    # ------------------------------------------------------------------
    # residue selection (same filter as GeometryDistance / Clashes)
    # ------------------------------------------------------------------

    def standard_residues(self):
        '''standard amino-acid residues (skip hetero / water).

        Cutoff: none - filters on residue id[0] == ' ' only.
        '''
        return [r for r in self.chain if r.id[0] == " "]

    # ------------------------------------------------------------------
    # chain continuity
    # ------------------------------------------------------------------

    def missing_residues(self) -> list:
        '''check if every residues has expected backbone atoms: N, CA, C, O
        return: residues with missing atoms
        Cutoff: no tolerance - any residue lacking N/CA/C/O is listed;
        Breaks.check() fails when num_missing > 0.
        '''
        missing = []
        for r in self.standard_residues():
            atoms = {a.name for a in r}
            required = {"N", "CA", "C", "O"}
            if not required <= atoms:
                missing.append((r.id[1], required - atoms))
        return missing

    def chain_continuity(self, cn_break: float = None,
            ca_break: float = None, cn_strain: float = None):
        '''Localise backbone breaks / chain-continuity failures for one chain.

        One row per consecutive standard-residue pair: real discontinuities
        land in 'breaks', stretched-but-intact peptides in 'strained';
        num_fragments = 1 + num_breaks (a clean chain is one fragment).

        Cutoff: cn_break=2.0 A (Breaks.cn_break), C-N above this is a real
        break (ideal 1.329 A); ca_break=4.5 A (Breaks.ca_break), CA-CA
        above this is a break (ideal ~3.8 A); cn_strain=1.5 A
        (Breaks.cn_strain), C-N above this but <= cn_break is a stretched
        peptide. A pair is also a break when C-N is unmeasurable
        (missing atoms) or residue numbering skips (resnum_gap != 1).
        '''
        cn_break = self.cn_break if cn_break is None else cn_break
        ca_break = self.ca_break if ca_break is None else ca_break
        cn_strain = self.cn_strain if cn_strain is None else cn_strain
        res = self.standard_residues()
        # breaks:  one row per real discontinuity (chain break)
        # strained: one row per stretched but intact peptide
        # frag: fragment counter, starts at 0, += 1 per break
        breaks, strained, frag = [], [], 0

        for i in range(len(res) - 1):
            a, b = res[i], res[i + 1]
            d_cn = float(np.linalg.norm(a['C'].coord - b['N'].coord)) \
                if ('C' in a and 'N' in b) else None
            d_ca = float(np.linalg.norm(a['CA'].coord - b['CA'].coord)) \
                if ('CA' in a and 'CA' in b) else None
            num_gap = b.id[1] - a.id[1]           # residue-numbering gap

            row = {
                'prev_res': a.id[1],
                'prev_name': a.resname,
                'next_res': b.id[1],
                'next_name': b.resname,
                'resnum_gap': num_gap,
                'd_C_N': d_cn,
                'd_CA_CA': d_ca
            }

            # pair is broken: C-N unmeasurable or > cn_break, resnum_gap != 1,
            # or d_CA_CA > ca_break
            # incremented exactly once for every break, so num_fragments = 1 + frag
            broken = (d_cn is None or d_cn > cn_break or num_gap != 1
                or (d_ca is not None and d_ca > ca_break))
            if broken:
                frag += 1
                breaks.append({**row, 'fragment_after': frag})
            # pair is not broken and 1.5 Å < d_C_N ≤ 2.0 Å
            elif d_cn is not None and d_cn > cn_strain:
                strained.append(row)

        # column-complete frames also when no pair qualifies
        pair_cols = ['prev_res', 'prev_name', 'next_res', 'next_name',
                     'resnum_gap', 'd_C_N', 'd_CA_CA']
        breaks = pd.DataFrame(breaks, columns=pair_cols + ['fragment_after'])
        strained = pd.DataFrame(strained, columns=pair_cols)
        return {
            'num_residues': len(res),
            'num_fragments': 1 + frag,
            'num_breaks': len(breaks),
            'num_strained': len(strained),
            'breaks': breaks,
            'strained': strained,
        }

    def ca_distances(self):
        '''detect bad Ca connections. No large gap
        calculate distance between Ca-Ca atoms
        Cutoff: ideal CA-CA ~3.8 A; Breaks.check() fails when
        ca_max_dist > 5 A.
        '''
        ca = [r["CA"].coord for r in self.standard_residues() if "CA" in r]
        if len(ca) < 2:
            return np.array([])
        # distance between ca-ca
        ca2 = np.asarray(ca[1:]) - np.asarray(ca[:-1])
        return np.linalg.norm(ca2, axis=1)

    def peptide_distances(self):
        """peptide C-N distance
        https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111178&ver=0&trm=peptide&org=
        Cutoff: ideal C-N = 1.329 A (peptide bond); chain_continuity()
        flags C-N > 2.0 A as a break; Breaks.check() fails when
        pep_max_dist > 5 A.
        """
        residues = self.standard_residues()
        pep_dist = [
            np.linalg.norm(residues[i]["C"].coord - residues[i+1]["N"].coord)
            for i in range(len(residues)-1)
            if "C" in residues[i] and "N" in residues[i+1]
        ]
        return np.array(pep_dist)


    # ------------------------------------------------------------------
    # summary + gates (used by pipelines/qc/chain_breaks.py)
    # ------------------------------------------------------------------

    def run(self, keep_details:bool=True) -> dict:
        '''compute all chain-continuity metrics of one chain.

        Covers the continuity half of Geometry.run() and the former
        backbone_qc_screen(): missing backbone atoms (num_missing),
        numbering-gap / broken peptides (num_breaks), the resulting
        fragment count (num_fragments = 1 + num_breaks), stretched
        peptides (num_strained) plus the max consecutive CA-CA
        (ca_max_dist, ideal ~3.8 A) and peptide C-N distance
        (pep_max_dist, ideal 1.329 A).

        Returns a dict of summary statistics (n_residues, num_missing,
        num_breaks, num_fragments, num_strained, ca_max_dist,
        pep_max_dist); raw tables are nested under 'details' when
        keep_details.
        '''
        res = self.standard_residues()
        missing = self.missing_residues()
        cont = self.chain_continuity()
        ca = self.ca_distances()
        pep = self.peptide_distances()
        info = {
            'n_residues': len(res),
            'num_missing': len(missing),
            'num_breaks': cont['num_breaks'],
            'num_fragments': cont['num_fragments'],
            'num_strained': cont['num_strained'],
            'ca_max_dist': round(float(ca.max()), 2) if len(ca) else 0.0,
            'pep_max_dist': round(float(pep.max()), 2) if len(pep) else 0.0,
        }
        if keep_details:
            info['details'] = {
                'missing_residues': missing,
                'breaks': cont['breaks'],
                'strained': cont['strained'],
            }
        return info

    def check(self, info:dict=None) -> tuple:
        '''gate on the summary statistics returned by run().

        Cutoffs (fail -> reason), see Breaks.GATES:
          num_missing > 0      -> 'missing_residues' (no tolerance);
          num_breaks > 0       -> 'chain_breaks' (broken / numbering-gap
                                  peptides);
          num_fragments > 1    -> 'chain_fragments' (one contiguous chain);
          num_strained > 3     -> 'chain_strain' (design tolerance);
          ca_max_dist > 5 A    -> 'ca_distance' (ideal ~3.8 A);
          pep_max_dist > 5 A   -> 'peptide_distance' (ideal 1.329 A).
        '''
        if info is None:
            info = self.run(keep_details=False)
        if not info:
            return False, 'null_breaks'
        gates = self.GATES
        # no missing atoms in any residues
        if info.get('num_missing', 0) > gates['num_missing']:
            return False, 'missing_residues'
        # chain continuity: any real break (broken / numbering-gap peptide)
        if info.get('num_breaks', 0) > gates['num_breaks']:
            return False, 'chain_breaks'
        # a broken chain splits into fragments: num_fragments = num_breaks + 1
        if info.get('num_fragments', 1) > gates['num_fragments']:
            return False, 'chain_fragments'
        # stretched-but-intact peptides (design tolerance)
        if info.get('num_strained', 0) > gates['num_strained']:
            return False, 'chain_strain'
        # good cutoff 3.8A, allowed 5A
        if info.get('ca_max_dist', 0.0) > gates['ca_max_dist']:
            return False, 'ca_distance'
        # ideal C-N 1.329 A; fail on a break (max > 5 A)
        if info.get('pep_max_dist', 0.0) > gates['pep_max_dist']:
            return False, 'peptide_distance'
        return True, 'pass'
