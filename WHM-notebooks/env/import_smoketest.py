#!/usr/bin/env python
"""WHM environment import smoke-test.

Mirrors the imports in the first code cell of WHMcode.ipynb and prints a
version + OK line for each. Exits non-zero if any required import fails.

classy (CLASS) is treated as OPTIONAL here: the CLASS build is handled in a
separate hand-off. If classy is not importable it is reported as PENDING
rather than failing the whole smoke-test. Set WHM_REQUIRE_CLASSY=1 to make a
missing classy a hard failure.

NOTE: this deliberately imports the LOCAL modified WHM-CAMB, not pip's camb.
"""
import importlib
import os
import sys

GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"


def _ver(mod):
    for attr in ("__version__", "version", "VERSION"):
        v = getattr(mod, attr, None)
        if isinstance(v, str) and v:
            return v
    # Fall back to installed-distribution metadata (e.g. euclidemu2 has no __version__).
    try:
        import importlib.metadata as _md
        return _md.version(getattr(mod, "__name__", "").split(".")[0])
    except Exception:  # noqa: BLE001
        return "(version unknown)"


failures = []


def check(label, fn, required=True):
    try:
        result = fn()
        print(f"{GREEN}[ OK ]{RESET} {label}: {result}")
    except Exception as e:  # noqa: BLE001
        tag = f"{RED}[FAIL]{RESET}" if required else f"{YELLOW}[PEND]{RESET}"
        print(f"{tag} {label}: {type(e).__name__}: {e}")
        if required:
            failures.append(label)


def imp(name):
    def _f():
        m = importlib.import_module(name)
        return _ver(m)
    return _f


print("=" * 60)
print(f"Python {sys.version.split()[0]}  ({sys.executable})")
import platform
print(f"Platform: {platform.platform()}  machine={platform.machine()}")
print("=" * 60)

check("euclidemu2", imp("euclidemu2"))
check("numpy", imp("numpy"))
check("scipy", imp("scipy"))
check("matplotlib", imp("matplotlib"))


def _camb():
    import camb
    # confirm this is the LOCAL modified WHM-CAMB with brieden modes
    from camb.nonlinear import halofit_brieden2025_tweaked  # noqa: F401
    return f"{camb.__version__}  (path={os.path.dirname(camb.__file__)}, brieden modes present)"


check("camb (local WHM-CAMB)", _camb)
check("baccoemu", imp("baccoemu"))


def _bacco_instantiate():
    import baccoemu
    baccoemu.Matter_powerspectrum()
    return "Matter_powerspectrum() instantiated"


check("baccoemu.Matter_powerspectrum()", _bacco_instantiate)


def _velo_lpt():
    from velocileptors.LPT.lpt_rsd_fftw import LPT_RSD  # noqa: F401
    import velocileptors
    return f"{_ver(velocileptors)}  LPT_RSD ok"


def _velo_moment():
    from velocileptors.LPT.moment_expansion_fftw import MomentExpansion  # noqa: F401
    return "MomentExpansion ok"


check("velocileptors.LPT.lpt_rsd_fftw.LPT_RSD", _velo_lpt)
check("velocileptors.LPT.moment_expansion_fftw.MomentExpansion", _velo_moment)


def _classy():
    from classy import Class  # noqa: F401
    import classy
    return _ver(classy)


check("classy (CLASS)", _classy, required=bool(os.environ.get("WHM_REQUIRE_CLASSY")))

# Not in WHMcode.ipynb cell 0, but a required WHM dependency per the project goal.
check("FlamingoBaryonResponseEmulator", imp("FlamingoBaryonResponseEmulator"))

print("=" * 60)
if failures:
    print(f"{RED}SMOKETEST FAILED{RESET}: {len(failures)} required import(s) failed: {failures}")
    sys.exit(1)
print(f"{GREEN}SMOKETEST PASSED{RESET}: all required imports succeeded.")
