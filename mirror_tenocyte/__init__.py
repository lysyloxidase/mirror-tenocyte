"""mirror_tenocyte — mixed-chirality (D/L) collagen matrices and tenocyte
mechanotransduction, in silico.

Despite the name, this is not really about mirror-image organisms or even about
pure mirror matrices. Integrin recognition is stereospecific, so a *pure* D
matrix presents no bindable epitope and the cell behaves as if on a
non-adhesive surface — that answer is immediate and uninteresting. Everything
worth modelling here concerns **mixtures**, which is a chiral biomaterials
question with a cell-biological readout.

Modules along the axis:

    chirality   chi -> bindable ligand density (stereospecific recognition)
    assembly    D/L composition -> mechanical competence (fibrillogenesis),
                plus the collagen assembly hierarchy
    model       six coupled ODEs: integrin -> FAK -> Rho -> tension -> YAP -> collagen
    reduced     what the model ACTUALLY is: a function of two scalars
    clutch      explicit motor-clutch in pN/nm/s — diagnostic only, deliberately
                NOT imported by `model`
    analysis    mode discrimination, observables, sensitivity

Important structural fact, verified in `reduced.verify_reduction`: the ODE
cascade reads the environment only through `effective_ligand` and
`effective_stiffness`, so every prediction is a function of those two numbers.
The cascade supplies the shape of that function and names measurable
intermediates, but adds no degrees of freedom.

Hypothesis-generating, not validated. No tenocyte-specific data is used.
See README.md for the current position and HISTORY.md for what was retracted.
"""

from .parameters import (
    get_default_parameters, get_default_environment, get_preset, PRESETS,
)
from .model import (
    derivatives, run_to_steady_state, run_timecourse, phenotype_index,
    classify_phenotype, STATE_NAMES,
)
from .chirality import effective_ligand, chirality_specific_gain
from .assembly import (
    assembly_competence, effective_stiffness, heterochiral_helix_probability,
)
from .reduced import ReducedModel, verify_reduction, env_coordinates

__all__ = [
    # parameters
    "get_default_parameters", "get_default_environment", "get_preset", "PRESETS",
    # model
    "derivatives", "run_to_steady_state", "run_timecourse",
    "phenotype_index", "classify_phenotype", "STATE_NAMES",
    # chirality
    "effective_ligand", "chirality_specific_gain",
    # assembly
    "assembly_competence", "effective_stiffness",
    "heterochiral_helix_probability",
    # the honest reduction
    "ReducedModel", "verify_reduction", "env_coordinates",
]

__version__ = "0.2.0"
