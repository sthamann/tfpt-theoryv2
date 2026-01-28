from __future__ import annotations

from tfpt_unconventional.modules.cosmo_history_sampler import CosmoHistorySamplerModule
from tfpt_unconventional.modules.flavor_holdout_search import FlavorHoldoutSearchModule
from tfpt_unconventional.modules.gravity_gaugefix_ga import GravityGaugeFixingGAModule
from tfpt_unconventional.modules.matching_metamorphic_audit import MatchingMetamorphicAuditModule
from tfpt_unconventional.modules.omega_b_aps_bridge import OmegaBApsBridgeModule
from tfpt_unconventional.modules.threshold_graph_audit import ThresholdGraphAuditModule
from tfpt_unconventional.modules.torsion_regime_designer import TorsionRegimeDesignerModule


def get_unconventional_module_registry():
    modules = [
        MatchingMetamorphicAuditModule(),
        CosmoHistorySamplerModule(),
        OmegaBApsBridgeModule(),
        ThresholdGraphAuditModule(),
        FlavorHoldoutSearchModule(),
        GravityGaugeFixingGAModule(),
        TorsionRegimeDesignerModule(),
    ]
    return {m.module_id: m for m in modules}

