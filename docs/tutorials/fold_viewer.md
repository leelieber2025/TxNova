# TxNova fold viewer

This is the 3Dmol widget TxNova writes as `fold/<locus>.html` when
`coding.fold: true`. Color is pLDDT (blue = confident). The model below
is a packaged example ORF so the widget still renders if ESMFold is
unreachable. It is not a GSE221720 residual.

The view is **inline** (same pattern as py3Dmol / omicverse): the
HTML/JavaScript lives on this page, not in an iframe. Drag to rotate,
scroll to zoom.

```{raw} html
:file: _fold_embed.html
```

```{raw} html
<p class="txnova-fold-note"><a href="../fold_viewer_example.html" target="_blank">Open as a standalone page</a></p>
```
