from __future__ import annotations

from tfpt_suite.modules.aps_eta_gluing import ApseEtaGluingModule
from tfpt_suite.modules.arrow_mechanism import ArrowMechanismModule
from tfpt_suite.modules.birefringence_tomography import BirefringenceTomographyModule
from tfpt_suite.modules.baryogenesis_mechanism import BaryogenesisMechanismModule
from tfpt_suite.modules.bbn_consistency import BbnConsistencyModule
from tfpt_suite.modules.boltzmann_transfer import BoltzmannTransferModule
from tfpt_suite.modules.brst_ghost_deriver import BrstGhostDeriverModule
from tfpt_suite.modules.core_invariants import CoreInvariantsModule
from tfpt_suite.modules.bbn_neff_sanity import BbnNeffSanityModule
from tfpt_suite.modules.dm_alternative_channels import DmAlternativeChannelsModule
from tfpt_suite.modules.flavor_topology_mapper import FlavorTopologyMapperModule
from tfpt_suite.modules.gw_background_predictor import GwBackgroundPredictorModule
from tfpt_suite.modules.mass_spectrum_minimal import MassSpectrumMinimalModule
from tfpt_suite.modules.mass_spectrum_deriver import MassSpectrumDeriverModule
from tfpt_suite.modules.koide_constraints import KoideConstraintsModule
from tfpt_suite.modules.qed_anomalies_audit import QedAnomaliesAuditModule
from tfpt_suite.modules.torsion_condensate import TorsionCondensateModule
from tfpt_suite.modules.ufe_gravity_normalization import UfeGravityNormalizationModule
from tfpt_suite.modules.chiral_index_three_cycles import ChiralIndexThreeCyclesModule
from tfpt_suite.modules.bounce_perturbations import BouncePerturbationsModule
from tfpt_suite.modules.primordial_spectrum_builder import PrimordialSpectrumBuilderModule
from tfpt_suite.modules.gw_background_bounds import GwBackgroundBoundsModule
from tfpt_suite.modules.discrete_consistency_uniqueness import DiscreteConsistencyUniquenessModule
from tfpt_suite.modules.defect_partition_derivation import DefectPartitionDerivationModule
from tfpt_suite.modules.alpha_on_shell_bridge import AlphaOnShellBridgeModule
from tfpt_suite.modules.discrete_complexity_minimizer import DiscreteComplexityMinimizerModule
from tfpt_suite.modules.effective_action_r2 import EffectiveActionR2Module
from tfpt_suite.modules.g2_and_lamb_shift_proxy import G2AndLambShiftProxyModule
from tfpt_suite.modules.mobius_cusp_classification import MobiusCuspClassificationModule
from tfpt_suite.modules.mobius_delta_calibration import MobiusDeltaCalibrationModule
from tfpt_suite.modules.topology_phase_map import TopologyPhaseMapModule
from tfpt_suite.modules.seesaw_block import SeesawBlockModule
from tfpt_suite.modules.omega_b_conjecture_scan import OmegaBConjectureScanModule
from tfpt_suite.modules.axion_fa_derivation import AxionFaDerivationModule
from tfpt_suite.modules.axion_dm_pipeline import AxionDmPipelineModule
from tfpt_suite.modules.axion_scenario_matrix import AxionScenarioMatrixModule
from tfpt_suite.modules.dark_energy_paths import DarkEnergyPathsModule
from tfpt_suite.modules.global_consistency_test import GlobalConsistencyTestModule
from tfpt_suite.modules.torsion_bounds_mapping import TorsionBoundsMappingModule
from tfpt_suite.modules.torsion_dm_pipeline import TorsionDmPipelineModule
from tfpt_suite.modules.torsion_observable_designer import TorsionObservableDesignerModule
from tfpt_suite.modules.torsion_falsifiability_snr import TorsionFalsifiabilitySnrModule
from tfpt_suite.modules.torsion_observable_spin_fluid import TorsionObservableSpinFluidModule
from tfpt_suite.modules.ckm_full_pipeline import CkmFullPipelineModule
from tfpt_suite.modules.mobius_z3_yukawa_generator import MobiusZ3YukawaGeneratorModule
from tfpt_suite.modules.pmns_full_pipeline import PmnsFullPipelineModule
from tfpt_suite.modules.flavor_joint_objective_scan import FlavorJointObjectiveScanModule
from tfpt_suite.modules.uncertainty_propagator import UncertaintyPropagatorModule
from tfpt_suite.modules.pmns_mechanism_bridge import PmnsMechanismBridgeModule
from tfpt_suite.modules.pmns_z3_breaking import PmnsZ3BreakingModule
from tfpt_suite.modules.k_calibration import KCalibrationModule
from tfpt_suite.modules.cosmo_threshold_history import CosmoThresholdHistoryModule
from tfpt_suite.modules.cosmo_reheating_policy_v106 import CosmoReheatingPolicyV106Module
from tfpt_suite.modules.msbar_matching_map import MsbarMatchingMapModule
from tfpt_suite.modules.matching_finite_pieces import MatchingFinitePiecesModule
from tfpt_suite.modules.below_mt_eft_cascade import BelowMtEftCascadeModule
from tfpt_suite.modules.stability_unitarity_audit import StabilityUnitarityAuditModule
from tfpt_suite.modules.likelihood_engine import LikelihoodEngineModule
from tfpt_suite.modules.alpha_precision_audit import AlphaPrecisionAuditModule
from tfpt_suite.modules.two_loop_rg_fingerprints import TwoLoopRgFingerprintsModule
from tfpt_suite.modules.unification_gate import UnificationGateModule
from tfpt_suite.modules.predictions_dashboard import PredictionsDashboardModule
from tfpt_suite.modules.qft_completeness_ledger import QftCompletenessLedgerModule
from tfpt_suite.modules.anomaly_cancellation_audit import AnomalyCancellationAuditModule
from tfpt_suite.modules.baryogenesis_placeholder import BaryogenesisPlaceholderModule
from tfpt_suite.modules.arrow_of_time_proxy import ArrowOfTimeProxyModule


def get_module_registry():
    modules = [
        CoreInvariantsModule(),
        UfeGravityNormalizationModule(),
        BrstGhostDeriverModule(),
        ChiralIndexThreeCyclesModule(),
        MassSpectrumMinimalModule(),
        MassSpectrumDeriverModule(),
        BbnNeffSanityModule(),
        BbnConsistencyModule(),
        KoideConstraintsModule(),
        TwoLoopRgFingerprintsModule(),
        UnificationGateModule(),
        AlphaPrecisionAuditModule(),
        PredictionsDashboardModule(),
        QftCompletenessLedgerModule(),
        AnomalyCancellationAuditModule(),
        EffectiveActionR2Module(),
        BouncePerturbationsModule(),
        GwBackgroundBoundsModule(),
        GwBackgroundPredictorModule(),
        CosmoThresholdHistoryModule(),
        CosmoReheatingPolicyV106Module(),
        KCalibrationModule(),
        PrimordialSpectrumBuilderModule(),
        BoltzmannTransferModule(),
        MatchingFinitePiecesModule(),
        MsbarMatchingMapModule(),
        BelowMtEftCascadeModule(),
        StabilityUnitarityAuditModule(),
        ApseEtaGluingModule(),
        DiscreteConsistencyUniquenessModule(),
        DefectPartitionDerivationModule(),
        AlphaOnShellBridgeModule(),
        DiscreteComplexityMinimizerModule(),
        MobiusCuspClassificationModule(),
        MobiusDeltaCalibrationModule(),
        TopologyPhaseMapModule(),
        SeesawBlockModule(),
        BirefringenceTomographyModule(),
        OmegaBConjectureScanModule(),
        AxionFaDerivationModule(),
        AxionDmPipelineModule(),
        AxionScenarioMatrixModule(),
        DarkEnergyPathsModule(),
        TorsionCondensateModule(),
        TorsionBoundsMappingModule(),
        TorsionObservableSpinFluidModule(),
        TorsionObservableDesignerModule(),
        TorsionFalsifiabilitySnrModule(),
        TorsionDmPipelineModule(),
        DmAlternativeChannelsModule(),
        CkmFullPipelineModule(),
        MobiusZ3YukawaGeneratorModule(),
        PmnsFullPipelineModule(),
        FlavorTopologyMapperModule(),
        FlavorJointObjectiveScanModule(),
        UncertaintyPropagatorModule(),
        PmnsMechanismBridgeModule(),
        PmnsZ3BreakingModule(),
        GlobalConsistencyTestModule(),
        LikelihoodEngineModule(),
        BaryogenesisPlaceholderModule(),
        BaryogenesisMechanismModule(),
        G2AndLambShiftProxyModule(),
        QedAnomaliesAuditModule(),
        ArrowOfTimeProxyModule(),
        ArrowMechanismModule(),
    ]
    return {m.module_id: m for m in modules}

