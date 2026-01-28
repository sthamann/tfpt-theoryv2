from __future__ import annotations

from typing import Dict

from theoryv3_suite.modules.axion_dm_audit import AxionDmAuditModule
from theoryv3_suite.modules.alpha_backreaction_sensitivity_audit import AlphaBackreactionSensitivityAuditModule
from theoryv3_suite.modules.baryon_consistency_audit import BaryonConsistencyAuditModule
from theoryv3_suite.modules.constant_factory_audit import ConstantFactoryAuditModule
from theoryv3_suite.modules.dark_energy_exponential_audit import DarkEnergyExponentialAuditModule
from theoryv3_suite.modules.dark_energy_norm_half_origin_audit import DarkEnergyNormHalfOriginAuditModule
from theoryv3_suite.modules.defect_partition_g5_audit import DefectPartitionG5AuditModule
from theoryv3_suite.modules.flavor_pattern_audit import FlavorPatternAuditModule
from theoryv3_suite.modules.g5_crosslink_audit import G5CrosslinkAuditModule
from theoryv3_suite.modules.g5_origin_audit import G5OriginAuditModule
from theoryv3_suite.modules.pmns_tm1_audit import PmnsTm1AuditModule
from theoryv3_suite.modules.seed_invariants_audit import SeedInvariantsAuditModule
from theoryv3_suite.modules.yukawa_exponent_index_audit import YukawaExponentIndexAuditModule
from theoryv3_suite.modules.yukawa_index_mapping_audit import YukawaIndexMappingAuditModule


def get_module_registry() -> Dict[str, object]:
    modules = [
        SeedInvariantsAuditModule(),
        DefectPartitionG5AuditModule(),
        AlphaBackreactionSensitivityAuditModule(),
        G5OriginAuditModule(),
        DarkEnergyExponentialAuditModule(),
        DarkEnergyNormHalfOriginAuditModule(),
        FlavorPatternAuditModule(),
        PmnsTm1AuditModule(),
        YukawaExponentIndexAuditModule(),
        YukawaIndexMappingAuditModule(),
        BaryonConsistencyAuditModule(),
        AxionDmAuditModule(),
        G5CrosslinkAuditModule(),
        ConstantFactoryAuditModule(),
    ]
    return {m.module_id: m for m in modules}
