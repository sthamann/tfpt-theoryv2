# file: e8cascade/solver_patch.py
from e8cascade.solver import E8CascadeSolver as Base

class SolverPatched(Base):
    # -----------------------------------------------------------------
    #  ►► NEU ◄◄  zwei getrennte Skalierungsfaktoren für das
    #             gravit. Portal (g1/g2 gemeinsam, g3 separat)
    # -----------------------------------------------------------------
    cR_factor  : float = 1.0   # wirkt auf g1 u. g2
    cR3_factor : float = 1.0   # zusätzl. Faktor nur für g3

    def set_thresholds(self, d):
        self.thresholds = d
        self.Sigma_F = d.get('MSigma', 1e9)
        self.N_R     = d.get('MNR',    1e12)

    def _beta_coeffs(self, mu, y):
        beta = super()._beta_coeffs(mu, y)

        # ------------------  Grav-Portal anpassen  -------------------
        if 'g1_grav' in beta:
            beta['g1'] += self.cR_factor * beta.pop('g1_grav')
        if 'g2_grav' in beta:
            beta['g2'] += self.cR_factor * beta.pop('g2_grav')
        if 'g3_grav' in beta:
            beta['g3'] += self.cR_factor * self.cR3_factor * beta.pop('g3_grav')

        # Schwellen
        if mu < self.Sigma_F:
            beta.pop('g3_sigma', None)
        if mu < self.N_R:
            beta.pop('g1_nr', None)
        return beta