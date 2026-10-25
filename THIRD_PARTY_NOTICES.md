# Licence scope and third-party resources

The MIT grant applies to original project code distributed in this repository, including the runtime-binding sources, environment/input checkers, historical builder and explicitly identified path derivatives. Original documentation, project-generated aggregate tables and original metadata curation are offered under CC BY 4.0, as stated in README. Source accessions, source-derived expression values and factual annotations retain their provenance; this grant does not claim exclusive rights in those facts. No licence in this release overrides rights in any external resource or source dataset. The retained dependency lock describes externally licensed packages and does not relicense them.

Upstream resources include DoRothEA, PROGENy, Reactome, GO/QuickGO annotations and Matrisome definitions. This release does not distribute their source files, gene-set memberships or downloadable package archives. It distributes original code, metadata, aggregate family counts/reference descriptors and retained feature-level anchor result tables, retaining collection labels and feature identifiers to describe the analysis. The two selected fgsea tables explicitly omit their derived `leadingEdge` membership columns; their other field strings and row order are unchanged. Users reconstructing upstream analyses must obtain the exact frozen resources separately and observe their applicable version-specific terms.

In particular, a software-package licence is not a blanket grant covering every constituent network. Reactome website terms must not be substituted for the declaration attached to a particular archived release. No permission to redistribute an unlicensed third-party membership file is implied.

The project's non-commercial academic use does not restrict the MIT grant for its original code to non-commercial use. The grant neither claims ownership of external material nor changes external data-use conditions.

The selected JGW03 count/sample/annotation tables derive from the three linked public GEO submissions; their source rights are not replaced by the code licence. See `data/jgw03_inputs/README.md` and the [GEO data-use disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). The frozen allocation registries are project-generated computational inputs, not additional biological samples.

This release excludes article PDFs, private diagnostic errors, machine-specific custody logs, credentials and software dependency tarballs. It is not the entire study evidence archive.

GSE277064 adds three unchanged source count gzip files and project-curated sample/title metadata. The associated Sato et al. article identifies the accession and is CC BY-NC-ND 4.0; see `data/g277_inputs/README.md`. No article content is redistributed, and no source-derived count is relicensed by MIT or the original-curation CC BY grant. The input command preserves fractional values and source identities.

G277 saved-preprocessing projections contain project-selected numerical inputs and mapped gene annotations. Source-derived values retain upstream provenance and applicable terms. They include no fitted effect/significance table, activity score or gene-set/network membership. Local resource reconstruction requires the separately authenticated frozen resources under their terms.

The PC12 dose-comparison extension distributes project-generated retained statistical result columns and identifiers. It excludes gene-set memberships, network edges and unsupported approximate edgeR interval fields. The original source accession and export/projection hashes are in `data/pc12_dose_saved/PROVENANCE.json`. Existing upstream licence boundaries remain unchanged.

The preparation-identity extension contains original interpretive curation and links to source records. It does not redistribute full articles, chemical workbooks or private review documents. Its CC BY 4.0 grant covers the original curation only.

GSE298414 preparation results add project-generated statistical tables, expression-QC summaries, inferred score matrices, mapped annotation identifiers and original source-label curation. No gene-set memberships, network edges, full article or supplemental workbook is redistributed. These results do not relicense the public GEO inputs or external resources. The separate source ledger preserves conflicting source labels rather than replacing original metadata.

GSE298414 constituent comparisons add original statistical outputs, public-library QC and inferred score matrices, representative annotation identifiers and curated metadata. They share 72 input libraries with the earlier preparation extension. No article files, supplemental workbook, gene-set memberships, network edges or scientific model binaries are added. External resources retain their respective terms.
