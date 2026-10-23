'''Empirical (Top8000-style) Ramachandran grids (MolProbity convention).

Split out of geometry_angle.py: residue-type specific 20x20 phi/psi count
grids (GLY / PRO / pre-PRO / general). Favoured = top 98% most populated
bins of that type, allowed = top 99.95% (MolProbity convention).
'''
import numpy as np


# ===========================================================================
# Empirical (Top8000-style) Ramachandran grids (MolProbity convention)
# ===========================================================================


class RamaGrid:
    '''Residue-type specific 20x20 phi/psi count grids with MolProbity
    favoured (98%) / allowed (99.95%) / outlier classification.

    grids: dict of GRID_TYPES -> (20, 20) integer count array (counts of a
    high-resolution reference set, e.g. the Richardson Top8000).
    '''

    # MolProbity plots are residue-type specific: GLY, PRO, pre-PRO and general
    GRID_TYPES = ('general', 'GLY', 'PRO', 'pre-PRO')
    NBIN = 20                       # 20 x 20 over phi/psi in [-180, 180)
    BINW = 360.0 / NBIN             # 18 degrees

    def __init__(self, grids):
        for key in self.GRID_TYPES:
            grids.setdefault(key, np.zeros((self.NBIN, self.NBIN),
                                           dtype=np.int64))
        self.grids = {k: np.asarray(v, dtype=np.int64)
                      for k, v in grids.items()}
        self.cuts = {k: self._cuts(self.grids[k]) for k in self.GRID_TYPES}

    @classmethod
    def _bin_index(cls, deg):
        '''phi/psi in [-180, 180) -> bin index in [0, NBIN);
        None if deg is None / NaN / inf'''
        if deg is None or not np.isfinite(deg):
            return None
        i = int(np.floor((deg + 180.0) / cls.BINW))
        return min(max(i, 0), cls.NBIN - 1)

    @staticmethod
    def rama_grid_key(resname, prev_resname=None):
        '''grid key for one residue: GLY / PRO / pre-PRO / general'''
        if resname == 'GLY':
            return 'GLY'
        if resname == 'PRO':
            return 'PRO'
        if prev_resname == 'PRO':
            return 'pre-PRO'
        return 'general'

    # -- learning the 98% / 99.95% density cuts -----------------------------
    @staticmethod
    def _cuts(counts):
        '''density thresholds such that favoured=top 98%, allowed=top 99.95%'''
        if counts.sum() <= 0:
            return None, None
        order = np.argsort(counts.ravel())[::-1]          # densest bins first
        cum = np.cumsum(counts.ravel()[order])
        total = cum[-1]
        fav_idx = int(np.searchsorted(cum, 0.980 * total))
        alw_idx = int(np.searchsorted(cum, 0.9995 * total))
        return counts.ravel()[order[fav_idx]], counts.ravel()[order[alw_idx]]

    @classmethod
    def fit(cls, df):
        '''learn grids from a DataFrame(resname, phi, psi [, prev_resname])'''
        grids = {k: np.zeros((cls.NBIN, cls.NBIN), dtype=np.int64)
                 for k in cls.GRID_TYPES}
        for _, row in df.iterrows():
            ip, is_ = cls._bin_index(row['phi']), cls._bin_index(row['psi'])
            if ip is None or is_ is None:
                continue
            key = cls.rama_grid_key(row['resname'],
                                    row.get('prev_resname', None))
            grids[key][ip, is_] += 1
        return cls(grids)

    @classmethod
    def load_npz(cls, path):
        '''load grids from an .npz with one (20,20) array per GRID_TYPES'''
        data = np.load(path)
        return cls({k: data[k] for k in cls.GRID_TYPES if k in data})

    @classmethod
    def demo(cls):
        '''synthetic grids for demo only - NOT real Top8000 data.'''
        rng = np.random.default_rng(0)
        grids = {}
        # approximate alpha / beta / alpha-L Gaussians per residue type
        blobs = {
            'general': [(-62, -41, 18, 22), (-115, 140, 25, 28)],
            'GLY': [(-62, -41, 18, 22), (-115, 140, 25, 28), (60, 40, 20, 20)],
            'PRO': [(-60, 140, 12, 25), (-60, -30, 12, 25)],
            'pre-PRO': [(-115, 140, 25, 28), (-60, -30, 20, 20)],
        }
        for key, means in blobs.items():
            grid = np.zeros((cls.NBIN, cls.NBIN))
            for (mp, ms, sp, ss) in means:
                n = 25000
                phi = rng.normal(mp, sp, n)
                psi = rng.normal(ms, ss, n)
                for p, s in zip(phi, psi):
                    ip, is_ = cls._bin_index(p), cls._bin_index(s)
                    if ip is not None and is_ is not None:
                        grid[ip, is_] += 1
            grids[key] = grid.astype(np.int64)
        return cls(grids)

    def classify(self, phi, psi, resname, prev_resname=None):
        ''''favoured' / 'allowed' / 'outlier' / 'NA' for one (phi, psi).

        Cutoffs: favoured if bin count >= favoured cut (top 98% density),
        allowed if >= allowed cut (top 99.95%), else outlier; None/NaN
        phi/psi -> 'NA'; empty grids (cuts None) -> 'allowed'.
        '''
        ip, is_ = self._bin_index(phi), self._bin_index(psi)
        if ip is None or is_ is None:
            return 'NA'
        key = self.rama_grid_key(resname, prev_resname)
        count = int(self.grids[key][ip, is_])
        fav_cut, alw_cut = self.cuts[key]
        if count <= 0 or fav_cut is None:
            return 'outlier' if count <= 0 else 'allowed'
        if count >= fav_cut:
            return 'favoured'
        if count >= alw_cut:
            return 'allowed'
        return 'outlier'

