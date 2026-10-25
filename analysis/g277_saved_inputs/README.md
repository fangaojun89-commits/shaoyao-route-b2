# Connect G277 saved preprocessing and original consumers

Reuse successful [count/design preparation](../g277_preparation/README.md) and [resource construction](../resource_build/README.md) outputs. This command joins saved preprocessing to the original annotation block, Reactome constructor and natural-course gene-method input validator. It never reruns a fit, score, filtering, TMM calculation, voom transformation or label allocation.

Use R4.6.1 and an explicit library containing data.table, limma, AnnotationDbi, org.Mm.eg.db **3.23.0**, jsonlite and digest. The actual tested versions/session are recorded. This uses an existing Windows library; no new Linux restore is claimed. Keep the exact verified frozen RDS under the resource root. Run from any working directory with a new output path:

```sh
python -I -B pipeline/connect_g277_saved_inputs.py --prepared-root /existing/g277-prepared --resource-root /existing/resource-build --rscript /path/to/Rscript --r-library /path/to/R-library --output /new/g277-native-inputs
```

The original annotation block must reproduce all bytes of the retained 55,416-row annotation table. Each of five domains must match unchanged raw counts, saved mask AND, original designs/contrasts, retained normalization metadata and saved voom E/weights. The original Reactome constructor must reproduce exact ordered index hashes; PROGENy/DoRothEA coverage uses minsize10 without ULM scoring. Three natural-course views must pass the unchanged original gene-method validator; invented negative-count, altered-design and zero-weight cases must be rejected.

Five `prepared/*_native_inputs.rds` outputs contain the saved preprocessing view and reconstructed Reactome indices, plus the separate frozen-resource identity. Annotation, commands, input/source hashes, actual exit, versions and session are retained. Locally reconstructed memberships remain under upstream terms and are never automatically published. These are inputs for separately requested future scientific execution. Portable full-model orchestration and whole-study cold restoration remain unfinished; this command does not repeat accepted production models.
