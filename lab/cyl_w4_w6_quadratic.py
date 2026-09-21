"""Joint degree-two scalar ladder: exact tangent tests, OpenAI Codex 2026-09-18.

Reuses the two promoted engines at da1e2df. No fitted source functions.
Default: through spin19, automatically through21 if consistent at19.
Exit2: exact violated identity; exit3: finite consistency only; exit1: failure.
--relation mode: exit0 for that identity alone, exit1 for a nonzero kernel.
"""
from pathlib import Path
import argparse,importlib.util,json,sys,time
from flint import fmpq as F,fmpq_mat

BASE={'U':(1,0,0),'XU':(1,1,0),'YU':(1,0,1),'XX':(1,2,0),'XY':(1,1,1),'YY':(1,0,2)}
KINDS=dict(BASE)
for q in (2,3,4):
    for prefix,(_,a,b) in BASE.items():
        KINDS[('U' if prefix=='U' else prefix if prefix.endswith('U') else prefix+'U')+str(2*q)]=(q,a,b)
TOP3={name:val for name,val in KINDS.items() if val[0]<=3}

def load(root):
    spec=importlib.util.spec_from_file_location('linear_core',root/'lab/cyl_w4_linear.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def dump(out,name,data):
    (out/name).write_text(json.dumps(data,indent=2)+'\n')

def locate(mo,spins,kernels,kinds,maxdrop):
    steps=[]
    for drop in range(maxdrop+1):
        b=mo.block(spins,kernels,kinds=kinds,maxdrop=maxdrop,lastdrop=drop,witness=True);steps.append(b)
        if b['rank']!=b['augmented_rank']:
            return {'spin':max(spins),'degree':(max(spins)+1)//2-drop,'drop':drop,'prefix_layers':steps,'relation':b['relations'][0]}
    return None

def compatible_relation(mo,spins,kernels,maxdrop):
    """Source-side null relation even when every target passes (not an obstruction)."""
    rows=[];entries=[];rhs=[];sm=max(spins);ref=min(spins);m=mo.m
    for s in spins:
        k=(s+1)//2
        align=m.moment_problem__beta_scale(s,ref)*mo.w**(3*(sm-s)//4)*(1-mo.w)**((sm-s)//2)
        for drop in range(min(k,maxdrop)+1):
            for a in range(k-drop+1):
                b=k-drop-a
                if (a,b)==(k,0):continue
                entry={}
                for kind in KINDS:
                    p=m.moment_problem__spatial_coefficient(kernels[kind,s],a,b)*align
                    for (j,_,__),v in p.to_dict().items():entry[kind,j]=v
                entries.append(entry);rows.append((s,a,b));rhs.append(mo.target(s)[a,b])
    cols=sorted(set().union(*(e.keys() for e in entries)));n=len(rows)
    rr,_=fmpq_mat([[e.get(c,F(0)) for c in cols]+[F(i==j) for j in range(n)] for i,e in enumerate(entries)]).rref()
    for i in range(n):
        if any(rr[i,j] for j in range(len(cols))):continue
        vec=[rr[i,len(cols)+j] for j in range(n)]
        if not any(vec):continue
        last=next(v for v in reversed(vec) if v)
        terms=[[str(v/last),*rows[j]] for j,v in enumerate(vec) if v]
        target=sum((v*b for v,b in zip(vec,rhs)),F(0))/last
        if target:raise ArithmeticError('unexpected incompatible relation after consistent ranks')
        return {'terms':terms,'target':'0','spins':spins,'maxdrop':maxdrop,'normalization':'last nonzero coefficient one; selected without target'}
    return None

def save_witness(mo,kernels,first,out,prefix):
    if not first:return None
    terms=first['relation']['terms'];verified={}
    for kind in KINDS:
        kernel=mo.relation(terms,lambda s,kind=kind:kernels[kind,s])
        if kernel:raise ArithmeticError('nonzero identity '+prefix+' '+kind)
        verified[kind]=True;dump(out,f'{prefix}_{kind}.json',{'kind':kind,'terms':terms})
    dump(out,prefix+'_natural.json',{'kind':'U','terms':first['relation']['natural_terms'],'target':first['relation']['natural_target'],'normalization':'last nonzero coefficient one'})
    bad=[row[:] for row in terms];bad[0][0]=str(2*F(bad[0][0]));tamper={}
    for kind in KINDS:
        p=mo.relation(bad,lambda s,kind=kind:kernels[kind,s]);tamper[kind]=mo.m.moment_problem__encode(p)
    if not any(tamper.values()):raise ArithmeticError('vacuous coefficient control')
    kind=next(k for k,v in tamper.items() if v)
    dump(out,prefix+'_tampered.json',{'kind':kind,'terms':bad})
    extensions={}
    for q in (2,3,4):
        for a in range(4):
            b=3-a;p=mo.relation(terms,lambda s,q=q,a=a,b=b:mo.response((s+1)//2,q,a,b))
            extensions[f'q{q}_X{a}_Y{b}']={'annihilated':not bool(p),'kernel':mo.m.moment_problem__encode(p)}
    return {'all_24_identities_zero':verified,'tampered_coefficient_kernels':tamper,'out_of_class_degree3_controls':extensions,'target_Sol1':'1','target_Sol2':'1/2'}

def run(root,out,stop19=False):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=True);lc=load(root);m=lc.load_core(root)
    lc.prepare(m,out);mo=lc.Moment(m,out);kernels={};blocks=[];topblocks=[];first=None;first3=None;zero_checks=[]
    last=19;s=3
    while s<=last:
        for kind,(q,a,b) in KINDS.items():kernels[kind,s]=mo.response((s+1)//2,q,a,b)
        spins=list(range(3 if s%4==3 else 5,s+1,4))
        b=mo.block(spins,kernels,kinds=KINDS,maxdrop=3);b3=mo.block(spins,kernels,kinds=TOP3,maxdrop=2)
        b3all=mo.block(spins,kernels,kinds=KINDS,maxdrop=2)
        if any(b3[x]!=b3all[x] for x in ('rank','augmented_rank','columns')):raise ArithmeticError('W8 top-three invisibility')
        zero_checks.append(s);blocks.append(b);topblocks.append(b3)
        print('EXACT',spins,'four layers',b['rank'],b['augmented_rank'],'top three',b3['rank'],b3['augmented_rank'],'seconds',round(time.monotonic()-start,3),flush=True)
        if first is None and b['rank']!=b['augmented_rank']:first=locate(mo,spins,kernels,KINDS,3)
        if first3 is None and b3['rank']!=b3['augmented_rank']:first3=locate(mo,spins,kernels,TOP3,2)
        if s==19 and first is None and first3 is None and not stop19:last=21
        s+=2
    # Old class relation, generated with its own complete old kind list.
    old=mo.block([3,7,11,15,19],kernels,kinds=lc.KINDS,maxdrop=2,witness=True)
    if (old['rank'],old['augmented_rank'])!=(78,79):raise ArithmeticError('old nine-kind baseline')
    oldterms=old['relations'][0]['terms'];escape={}
    for kind in lc.KINDS:
        if mo.relation(oldterms,lambda s,kind=kind:kernels[kind,s]):raise ArithmeticError('old identity')
    for kind in ('XXU4','XXU6'):
        p=mo.relation(oldterms,lambda s,kind=kind:kernels[kind,s]);escape[kind]=m.moment_problem__encode(p)
        if not p:raise ArithmeticError('new quadratic kind did not break old identity')
    dump(out,'linear_baseline.json',{'status':'REPLAY','block':old,'new_quadratic_kinds_escape':escape})
    omitted={}
    for q,g in [(2,2),(3,1)]:
        k=5;p=mo.Z
        for gg,ex,v,pols in mo.expanded:
            if gg!=g:continue
            for u,pol in enumerate(pols):
                ell=k-gg-u-q-1
                if ell>=0:p+=v*m.moment_problem__binom(ex,ell)*pol*mo.V**ell
        p=m.moment_problem__grade(p*mo.X**2,max(0,k-3),k)/m.moment_problem__normalizer(k)
        if not p:raise ArithmeticError('omitted-grade control vacuous')
        omitted[f'q{q}_omit_g{g}_spin9']=m.moment_problem__encode(p)
    result={'status':'EXACT','fit':'none','momentum_bound':2,'kinds':KINDS,'spins_through':last,'blocks':blocks,'top_three_prefix_blocks':topblocks,'first_failure':first,'first_top_three_failure':first3,'W8_top_three_invisibility_checked_at':zero_checks,'omitted_grade_controls':omitted,'scope':'analytic tangent zero at regular pillow exact points; fixed A1 through kappa^0; same continued-Beta/IBP functional; degree<=2 at every spectral order'}
    result['witnesses']={}
    for label,failure in [('relation_quadratic_four',first),('relation_quadratic_top3',first3)]:
        if failure:result['witnesses'][label]=save_witness(mo,kernels,failure,out,label)
    result['top_three_layer_blocks']=[mo.block(list(range(start_s,last+1,4)),kernels,kinds=TOP3,maxdrop=2,witness=True) for start_s in (3,5)]
    result['lower_layer_freedom']=[]
    for start_s in (3,5):
        spins=list(range(start_s,last+1,4))
        upper=mo.block(spins,kernels,kinds=KINDS,maxdrop=1)
        third=mo.block(spins,kernels,kinds=KINDS,maxdrop=2)
        fourth=mo.block(spins,kernels,kinds=KINDS,maxdrop=3)
        result['lower_layer_freedom'].append({'spins':spins,'upper_two_rows':upper['rows'],'upper_two_rank':upper['rank'],'kminus2_rows':third['rows']-upper['rows'],'kminus2_rank_gain':third['rank']-upper['rank'],'kminus3_rows':fourth['rows']-third['rows'],'kminus3_rank_gain':fourth['rank']-third['rank'],'all_lower_targets_formally_free':fourth['rank']-upper['rank']==fourth['rows']-upper['rows']})
    result['compatible_relations']=[]
    if first is None and first3 is None:
        for start_s in (3,5):
            rel=compatible_relation(mo,list(range(start_s,last+1,4)),kernels,3)
            if not rel:continue
            label=f'relation_compatible_mod{start_s%4}'
            for kind in KINDS:
                if mo.relation(rel['terms'],lambda s,kind=kind:kernels[kind,s]):raise ArithmeticError('compatible identity '+kind)
                dump(out,f'{label}_{kind}.json',{'kind':kind,'terms':rel['terms'],'target':'0','normalization':rel['normalization']})
            bad=[row[:] for row in rel['terms']];bad[0][0]=str(2*F(bad[0][0]));badkernels={}
            for kind in KINDS:
                badkernels[kind]=m.moment_problem__encode(mo.relation(bad,lambda s,kind=kind:kernels[kind,s]))
            if not any(badkernels.values()):raise ArithmeticError('compatible relation control vacuous')
            badkind=next(k for k,v in badkernels.items() if v)
            dump(out,f'{label}_tampered.json',{'kind':badkind,'terms':bad})
            rel['label']=label;rel['bad_kind']=badkind;rel['tampered_kernels']=badkernels
            rel['out_of_class_degree3_controls']={}
            # These sources require at most grade3 in the four tested layers.
            # q=1,degree3 would need grade4 and is deliberately not claimed.
            for q in (2,3,4):
                for a in range(4):
                    bb=3-a
                    value=mo.relation(rel['terms'],lambda s,q=q,a=a,bb=bb:mo.response((s+1)//2,q,a,bb))
                    rel['out_of_class_degree3_controls'][f'q{q}_X{a}_Y{bb}']={'annihilated':not bool(value),'kernel':m.moment_problem__encode(value)}
            result['compatible_relations'].append(rel)
    # Same independently implemented full Riccati engine, expanded kind inventory.
    lc.KINDS=KINDS.copy()
    independent=lc.independent_check(m,mo,kernels,result,out,last)
    if 'nine_identities_zero' in independent:independent['all_24_identities_zero']=independent.pop('nine_identities_zero')
    for b in independent.get('top_three_layer_relations',[]):b['all_24_zero']=b.pop('all_nine_zero')
    # These target-zero relations follow independently because all 24 kernels
    # and all A1 tangents have just matched the full Riccati engine.
    independent['compatible_relations_via_verified_kernels']=[{'label':r['label'],'terms':len(r['terms']),'target':r['target']} for r in result['compatible_relations']]
    # Re-evaluate every lower-layer target without reading Delta_A1: A1
    # derivatives in this subtraction were independently verified above.
    checked_targets={};target_inventory=[];tv=F(-1,3)
    for s in range(3,last+1,2):
        k=(s+1)//2
        for drop in range(min(k,3)+1):
            for a in range(k-drop+1):
                b=k-drop-a;value=F(0)
                if drop>=2:
                    for field,sign in [('Ward',1),('A1',-1)]:
                        nn,dd=map(m.check_a1__univariate,m.check_a1__coefficient(out,1,k,a,drop,field))
                        value+=sign*(nn.derivative()(tv)*dd(tv)-nn(tv)*dd.derivative()(tv))/dd(tv)**2
                    target_inventory.append([s,a,b])
                if value!=mo.target(s)[a,b]:raise ArithmeticError('Ward-minus-A1 target differs from Delta_A1')
                checked_targets[s,a,b]=value
    for rel in result['compatible_relations']:
        if sum((F(c)*checked_targets[s,a,b] for c,s,a,b in rel['terms']),F(0)):raise ArithmeticError('compatible target independently nonzero')
    independent['separate_Ward_minus_A1_target_checks']=len(target_inventory)
    independent['separate_target_inventory']=target_inventory
    dump(out,'lower_layer_freedom.json',{'status':'EXACT','meaning':'rank of lower rows modulo upper rows; full row gain means arbitrary lower-layer target values are finitely consistent, not predicted','blocks':result['lower_layer_freedom']})
    dump(out,'quadratic_independent.json',independent)
    enc=m.moment_problem__encode
    dump(out,'quadratic_kernels.json',{f'{kind}:{s}':enc(p) for (kind,s),p in kernels.items()})
    dump(out,'quadratic_targets.json',{str(s):{f'{a},{b}':str(v) for (a,b),v in mo.target(s).items()} for s in range(3,last+1,2)})
    result['exit']=2 if first or first3 else 3
    result['independent_Riccati_verified']=True
    dump(out,'quadratic_results.json',result)
    compact=lambda f:None if f is None else {**{key:f[key] for key in ('spin','degree','drop')},'terms':len(f['relation']['terms']),'natural_target':f['relation']['natural_target']}
    summary={'status':'EXACT','exit':result['exit'],'spins_through':last,'first_failure':compact(first),'first_top_three_failure':compact(first3),'independent_A1_values':independent['A1_value_checks'],'independent_A1_tangents':independent['A1_tangent_checks'],'independent_kernels':independent['kernel_checks'],'compatible_relations':[{'spins':r['spins'],'terms':len(r['terms']),'target':r['target'],'label':r['label'],'bad_kind':r['bad_kind']} for r in result['compatible_relations']],'scope':result['scope'],'higher_orders':('P=2 scalar ladder closed by identities through spin '+str(min(f['spin'] for f in (first,first3) if f))+': q>=5 cannot affect k-3; q>=4 cannot affect k-2. Unrestricted degree remains open.' if result['exit']==2 else 'NOT CLOSED for P=2; finite consistency only. Previously closed P=1 class is unchanged.'),'operator':'NOT DONE; no all-spin solution or matrix computation','elapsed_seconds':round(time.monotonic()-start,3)}
    dump(out,'quadratic_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    return result['exit']

def verify(root,out,path):
    lc=load(root);lc.KINDS=KINDS.copy()
    return lc.verify_relation(root,out,path)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--out',type=Path,required=True);p.add_argument('--stop-at19',action='store_true');p.add_argument('--relation',type=Path)
    a=p.parse_args();sys.exit(verify(a.root.resolve(),a.out.resolve(),a.relation.resolve()) if a.relation else run(a.root.resolve(),a.out.resolve(),a.stop_at19))
