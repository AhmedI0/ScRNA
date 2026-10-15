"""Run the CellBender 0.3.2 command line under recent PyTorch (used by run_cellbender.py).

Pyro's parameter store tags every parameter it hands out with a weak reference
to its unconstrained value (``param.unconstrained``). Recent PyTorch pickles a
parameter's Python attributes together with its data, and weak references
cannot be pickled, so CellBender's checkpoint (``torch.save`` of the whole
model) fails and the run stops before the posterior is computed. The tags are
removed only while ``torch.save`` runs and put back afterwards (pyro would set
them again on the next parameter access anyway); the model is not changed.

    python cellbender_compat.py remove-background --input ... --output ...
"""
from __future__ import annotations

import sys

import torch

_torch_save = torch.save


def _tagged_tensors(obj) -> list[torch.Tensor]:
    import pyro

    tensors = list(obj.parameters()) if isinstance(obj, torch.nn.Module) else []
    store = pyro.get_param_store()
    return tensors + [store[name] for name in store.keys()] + list(store._params.values())


def _save_without_weakrefs(obj, f, *args, **kwargs):
    held = []
    for t in _tagged_tensors(obj):
        ref = t.__dict__.pop("unconstrained", None)
        if ref is not None:
            held.append((t, ref))
    try:
        return _torch_save(obj, f, *args, **kwargs)
    finally:
        for t, ref in held:
            t.unconstrained = ref


torch.save = _save_without_weakrefs

if __name__ == "__main__":
    from cellbender.base_cli import main

    sys.argv[0] = "cellbender"
    sys.exit(main())
