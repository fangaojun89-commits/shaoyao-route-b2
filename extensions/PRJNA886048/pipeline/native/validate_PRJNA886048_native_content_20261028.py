"""Read-only native content gate. Does not authenticate execution or publish acceptance."""
import argparse,csv,hashlib,json,math
from pathlib import Path

TX2GENE_SHA256='b4bd74bd2e0056d950940686e31bcec66eacdcea2873d4a70cb1df431ee0fda5'

def require(value,message):
    if not value:raise ValueError(message)
def integer(v):return type(v) is int and v>=0
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def validate(native,tx2gene):
    reference_sha=hashlib.sha256(tx2gene.read_bytes()).hexdigest()
    require(reference_sha==TX2GENE_SHA256,'Transcript reference identity differs')
    qc=load(native/'fastp.json'); info=load(native/'quant/run_info.json')
    require(info.get('kallisto_version')=='0.46.1' and info.get('index_version')==10,'Wrong kallisto/index version')
    require(type(info.get('n_bootstraps')) is int and info['n_bootstraps']==0,'Unexpected bootstraps')
    for key in ['n_targets','n_processed','n_pseudoaligned','n_unique']:
        require(integer(info.get(key)),'Invalid count: '+key)
    require(info['n_targets']==146190 and 0<info['n_processed'] and 0<=info['n_unique']<=info['n_pseudoaligned']<=info['n_processed'],'Invalid quantification counts')
    before=qc['summary']['before_filtering']; after=qc['summary']['after_filtering']
    for summary in [before,after]:
        for k in ['total_reads','total_bases','q20_bases','q30_bases']:require(integer(summary[k]),'Invalid FASTP count')
        require(0<=summary['q30_bases']<=summary['q20_bases']<=summary['total_bases'],'Invalid quality base ordering')
        for k in ['q20_rate','q30_rate','gc_content']:
            require(type(summary[k]) in (int,float) and math.isfinite(summary[k]) and 0<=summary[k]<=1,'Invalid quality fraction')
        require(summary['q30_rate']<=summary['q20_rate'],'Invalid quality rate ordering')
    raw=before['total_reads']; retained=after['total_reads']
    require(raw%2==0 and retained%2==0 and 0<retained<=raw,'Invalid pair accounting')
    filtering=qc['filtering_result']
    require(all(integer(v) for v in filtering.values()),'Invalid filtering count')
    require(sum(filtering.values())==raw and filtering['passed_filter_reads']==retained,'Filtering sum mismatch')
    require(retained==2*info['n_processed'],'FASTP/kallisto pair mismatch')
    with tx2gene.open(encoding='utf-8-sig') as f:
        reader=csv.reader(f,delimiter='\t');header=next(reader)
        require(header==['transcript_id_versioned','transcript_id','gene_id','gene_name','length'],'Wrong transcript reference columns')
        expected=[row[0] for row in reader]
    require(len(expected)==146190 and len(set(expected))==146190,'Invalid transcript reference')
    ids=[];counts=[];tpms=[]
    with (native/'quant/abundance.tsv').open(encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        require(reader.fieldnames==['target_id','length','eff_length','est_counts','tpm'],'Wrong abundance columns')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()),'Malformed abundance row')
            values=[float(row[k]) for k in ['length','eff_length','est_counts','tpm']]
            require(all(math.isfinite(v) and v>=0 for v in values) and values[0]>0,'Invalid abundance value')
            ids.append(row['target_id']); counts.append(values[2]); tpms.append(values[3])
    require(len(ids)==len(set(ids))==146190 and set(ids)==set(expected),'Transcript universe differs')
    c=math.fsum(counts);t=math.fsum(tpms)
    require(abs(c-info['n_pseudoaligned'])<max(1,info['n_pseudoaligned']*1e-6),'Estimated count sum mismatch')
    require(abs(t-1e6)<5,'TPM sum mismatch')
    return {'status':'PASS_NATIVE_CONTENT_ONLY','raw_pairs':raw//2,'retained_pairs':retained//2,'pseudoaligned_pairs':info['n_pseudoaligned'],'targets':len(ids),'estimated_count_sum':c,'tpm_sum':t,'tx2gene_sha256':hashlib.sha256(tx2gene.read_bytes()).hexdigest(),'acceptance_published':False,'execution_provenance_verified':False}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True);p.add_argument('--tx2gene',type=Path,required=True);a=p.parse_args()
    print(json.dumps(validate(a.native,a.tx2gene),indent=2))
