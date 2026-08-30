"""TxNova: unannotated spliced residual loci from bulk RNA-seq BAMs."""

__version__ = "0.1.10"
USER_AGENT = f"txnova/{__version__}"

from txnova.errors import TxNovaError

__all__ = ["USER_AGENT", "TxNovaError", "__version__"]
