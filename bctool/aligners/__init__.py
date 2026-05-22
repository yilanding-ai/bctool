from .base import AlignerBase
from .bwameth import BwaMethAligner
from .bsbolt import BSBoltAligner
from .bsmap import BSMAPAligner
from .walt import WaltAligner
from .abismal import AbismalAligner
from .batmeth2 import Batmeth2Aligner
from .basal import BasalAligner
from .hisat3n import Hisat3nAligner, Hisat3nRepeatAligner
from .bismark import BismarkBwt2E2EAligner, BismarkHis2Aligner
from .bsseeker2 import (
    BSseeker2BwtAligner,
    BSseeker2Soap2Aligner,
    BSseeker2Bwt2E2EAligner,
    BSseeker2Bwt2LocalAligner,
)
from .conversion_aware import ConversionAwareAligner

ALIGNER_REGISTRY = {
    "bwameth": BwaMethAligner,
    "bsbolt": BSBoltAligner,
    "bsmap": BSMAPAligner,
    "walt": WaltAligner,
    "abismal": AbismalAligner,
    "batmeth2": Batmeth2Aligner,
    "basal": BasalAligner,
    "hisat3n": Hisat3nAligner,
    "hisat3n_repeat": Hisat3nRepeatAligner,
    "bismark_bwt2_e2e": BismarkBwt2E2EAligner,
    "bismark_his2": BismarkHis2Aligner,
    "bsseeker2_bwt": BSseeker2BwtAligner,
    "bsseeker2_soap2": BSseeker2Soap2Aligner,
    "bsseeker2_bwt2_e2e": BSseeker2Bwt2E2EAligner,
    "bsseeker2_bwt2_local": BSseeker2Bwt2LocalAligner,
    "conversion_aware": ConversionAwareAligner,
}


def get_aligner(name):
    cls = ALIGNER_REGISTRY.get(name)
    if cls is None:
        raise ValueError(f"Unknown aligner: {name}. Use --list-aligners to see available ones.")
    return cls
