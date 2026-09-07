"""TxNova: unannotated spliced residual loci from bulk RNA-seq BAMs."""

from txnova.errors import TxNovaError

__version__ = "0.1.11"
USER_AGENT = f"txnova/{__version__}"


__all__ = ["USER_AGENT", "TxNovaError", "__version__"]
