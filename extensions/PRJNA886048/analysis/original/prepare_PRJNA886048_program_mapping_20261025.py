"""Map frozen resources to the new reference without reading expression results."""
import collections
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = Path('R:/Shaoyao_storage/data/PRJNA886048/reference_preparation')
OUT = Path('R:/Shaoyao_storage/data/PRJNA886048/program_mapping_20261025')


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write(name, rows, fields):
    with (OUT / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    resources = ROOT / 'data/processed/gene_sets'
    paths = dict(genes=REF / 'genes.tsv', targets=REF / 'kallisto_tx2gene.tsv',
                 orthologs=resources / 'ncbi_human_mouse_rat_one_to_one.tsv',
                 Reactome=resources / 'reactome97_human_long.tsv',
                 PROGENy=resources / 'progeny_mouse_top500.tsv',
                 DoRothEA=resources / 'dorothea_mouse_AB.tsv')
    genes = read(paths['genes'])
    assert len({r['gene_id'] for r in genes}) == len(genes)
    represented = {r['gene_id'] for r in read(paths['targets'])}
    assert represented <= {r['gene_id'] for r in genes}
    # Ambiguity is determined in the complete annotation, not by expression rank.
    symbol_ids = collections.defaultdict(set)
    for r in genes:
        if r['gene_name']:
            symbol_ids[r['gene_name']].add(r['gene_id'])
    orth = [r for r in read(paths['orthologs']) if r['other_tax_id'] == '10090']
    assert all(r['human_tax_id'] == '9606' for r in orth)
    assert len({r['human_symbol'] for r in orth}) == len(orth)
    assert len({r['other_symbol'] for r in orth}) == len(orth)
    human_mouse = {r['human_symbol']: r['other_symbol'] for r in orth}
    edges, excluded = [], []
    membership = collections.defaultdict(set)
    names = {}
    for collection in ('Reactome', 'PROGENy', 'DoRothEA'):
        for r in read(paths[collection]):
            if collection == 'Reactome':
                assert r['tax_id'] == '9606' and r['release'] == '97'
                term, label, symbol = r['term_id'], r['term_name'], r['gene_symbol']
                mouse, weight, confidence = human_mouse.get(symbol, ''), '', ''
            elif collection == 'PROGENy':
                assert r['tax_id'] == '10090' and r['release'] == '1.34.0'
                term = label = r['pathway']
                symbol = mouse = r['gene_symbol']
                weight, confidence = r['weight'], ''
            else:
                assert r['tax_id'] == '10090' and r['release'] == '1.23.0'
                assert r['confidence'] in ('A', 'B')
                term = label = r['tf']
                symbol = mouse = r['target']
                weight, confidence = r['mor'], r['confidence']
            key = (collection, term)
            if key in names:
                assert names[key] == label
            names[key] = label
            membership[key].add(symbol)
            ids = symbol_ids.get(mouse, set())
            reason = ('NO_FROZEN_ONE_TO_ONE_ORTHOLOG' if not mouse else
                      'NO_EXACT_ENSEMBL_SYMBOL' if not ids else
                      'AMBIGUOUS_ENSEMBL_SYMBOL' if len(ids) != 1 else
                      'NO_REFERENCE_SEQUENCE' if not ids <= represented else '')
            base = dict(collection=collection, term_id=term, term_name=label,
                        source_symbol=symbol, mouse_symbol=mouse)
            if reason:
                excluded.append(dict(base, reason=reason))
            else:
                edges.append(dict(base, gene_id=next(iter(ids)), weight=weight,
                                  confidence=confidence))
    # Preserve signed footprints separately from unsigned set membership.
    assert len({(r['collection'], r['term_id'], r['gene_id']) for r in edges}) == len(edges)
    counts = collections.Counter((r['collection'], r['term_id']) for r in edges)
    coverage = [dict(collection=c, term_id=t, term_name=names[(c, t)],
                     source_members=len(membership[(c, t)]), reference_members=counts[(c, t)],
                     reference_size_10_to_500=10 <= counts[(c, t)] <= 500)
                for c, t in sorted(names)]
    OUT.mkdir(exist_ok=False)
    write('mapped_edges.tsv', edges, ['collection', 'term_id', 'term_name', 'source_symbol',
                                      'mouse_symbol', 'gene_id', 'weight', 'confidence'])
    write('excluded_edges.tsv', excluded, ['collection', 'term_id', 'term_name', 'source_symbol',
                                          'mouse_symbol', 'reason'])
    write('reference_coverage.tsv', coverage, ['collection', 'term_id', 'term_name',
                                              'source_members', 'reference_members',
                                              'reference_size_10_to_500'])
    summary = {}
    for c in ('Reactome', 'PROGENy', 'DoRothEA'):
        rows = [r for r in coverage if r['collection'] == c]
        summary[c] = dict(terms=len(rows), reference_edges=sum(r['reference_members'] for r in rows),
                          terms_10_to_500_at_reference=sum(r['reference_size_10_to_500'] for r in rows),
                          excluded_edges_by_reason=dict(collections.Counter(
                              r['reason'] for r in excluded if r['collection'] == c)))
    receipt = dict(status='REFERENCE_MAPPING_COMPLETE_EXPRESSION_ELIGIBILITY_PENDING',
                   sources={k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                            for k, p in paths.items()},
                   represented_reference_genes=len(represented), collections=summary,
                   rules=['Exact symbols; no outcome-based alias selection',
                          'Exclude ambiguous symbols using the full annotation',
                          'Reactome uses frozen strict human-mouse one-to-one orthologs',
                          'Native mouse PROGENy/DoRothEA weights and signs retained',
                          'Reapply 10-500 size bounds after expression filtering; retain all terms now'],
                   expression_data_read=False, biological_model_run=False)
    (OUT / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
