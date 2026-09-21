"""Nine-function scalar tangent test through spin 19, with exact adjoint identities.

OpenAI Codex, 2026-09-17. Sealed at 2d568a2. No fitted source functions.
Imports the previously promoted all-spin A1 certificate from --root/lab.
All generated files go under --out. The live research tree must be read-only.
An exit 2 excludes the stated analytic W2 + momentum-linear W4 truncation.
The default spin-19 run also bounds added higher orders of momentum degree <=1.
An exit 3 means finite consistency only, never an all-spin construction.
Argparse usage errors exit 64. Embedded provenance-only fields, including
those loaded through the A1 core, are consumed by no check; tampering with
those metadata fields is invisible by design.
"""
from pathlib import Path
import json,time,sys,importlib.util
from flint import fmpq as F,fmpq_mat

KINDS={'U':(1,0,0),'XU':(1,1,0),'YU':(1,0,1),'XX':(1,2,0),'XY':(1,1,1),'YY':(1,0,2),
       'U4':(2,0,0),'XU4':(2,1,0),'YU4':(2,0,1)}

def load_core(root):
    p=root/'lab/cyl_all_spin_a1_deficit.py'
    spec=importlib.util.spec_from_file_location('a1_core',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def dump(out,name,data):
    (out/name).write_text(json.dumps(data,indent=2)+'\n')

def prepare(core,out):
    source=out/'source';source.mkdir(exist_ok=True)
    for name,data in core.EMBEDDED_INPUTS.items():dump(source,name,data)
    core.derive_a1__derive(out);core.project_a1__run(out,source);core.check_a1__run(out);core.upper_layers_check__run(out)

class Moment:
    def __init__(self,m,out):
        self.m=m;self.out=out;self.C=m.moment_problem__C
        self.w,self.X,self.Y=self.C.gens();self.Z=self.C.constant(0)
        self.V,self.expanded=m.moment_problem__setup(out)
        self.cache={};self.base={};self.targets={}
    def response(self,k,q,a,b):
        m=self.m;key=k,q
        if key not in self.cache:
            p=self.Z
            for g,ex,v,pols in self.expanded:
                for u,pol in enumerate(pols):
                    ell=k-g-u-q-1
                    if ell>=0:p+=v*m.moment_problem__binom(ex,ell)*pol*self.V**ell
            self.cache[key]=p
        top=m.moment_problem__normalizer(k)
        if k not in self.base:
            base=self.Z
            for d in range(min(3,k)+1):
                for aa in range(k-d+1):
                    bb=k-d-aa
                    val=m.check_a1__value(m.check_a1__coefficient(self.out,1,k,aa,d),F(-1,3))
                    base+=val*self.X**aa*self.Y**bb
            self.base[k]=base
        raw=m.moment_problem__grade(self.cache[key]*self.X**a*self.Y**b,max(0,k-3),k)/top
        return raw-m.moment_problem__spatial_coefficient(raw,k,0)*self.base[k]
    def target(self,s):
        if s not in self.targets:self.targets[s]=self.m.moment_problem__targets((s+1)//2,self.out)
        return self.targets[s]
    def block(self,spins,kernels,kinds=KINDS,maxdrop=3,lastdrop=None,witness=False):
        m=self.m;sm=max(spins);ref=min(spins);rows=[];entries=[];rhs=[]
        for s in spins:
            k=(s+1)//2;limit=lastdrop if lastdrop is not None and s==sm else maxdrop
            for d in range(min(k,limit)+1):
                for a in range(k-d+1):
                    b=k-d-a
                    if (a,b)==(k,0):continue
                    entry={}
                    align=m.moment_problem__beta_scale(s,ref)*self.w**(3*(sm-s)//4)*(1-self.w)**((sm-s)//2)
                    for kind in kinds:
                        p=m.moment_problem__spatial_coefficient(kernels[kind,s],a,b)*align
                        for (j,_,__),v in p.to_dict().items():entry[kind,j]=v
                    rows.append((s,a,b));entries.append(entry);rhs.append(self.target(s)[a,b])
        cols=sorted(set().union(*(x.keys() for x in entries)));array=[[x.get(c,F(0)) for c in cols] for x in entries]
        rank=fmpq_mat(array).rank();aug=fmpq_mat([a+[b] for a,b in zip(array,rhs)]).rank()
        rec={'spins':spins,'rows':len(rows),'columns':len(cols),'rank':rank,'augmented_rank':aug,'maxdrop_previous':maxdrop,'lastdrop':lastdrop,'relations':[]}
        if witness and rank!=aug:
            nr=len(rows);rr,_=fmpq_mat([a+[F(i==j) for j in range(nr)] for i,a in enumerate(array)]).rref()
            for i in range(nr):
                if any(rr[i,j] for j in range(len(cols))):continue
                vec=[rr[i,len(cols)+j] for j in range(nr)]
                residual=sum((v*b for v,b in zip(vec,rhs)),F(0))
                if not residual:continue
                # First normalize without the target: the last nonzero coefficient is 1.
                last=next(v for v in reversed(vec) if v)
                natural=[[str(v/last),*rows[j]] for j,v in enumerate(vec) if v]
                terms=[[str(v/residual),*rows[j]] for j,v in enumerate(vec) if v]
                rec['relations'].append({'terms':terms,'target':'1','natural_terms':natural,'natural_target':str(residual/last)})
                break
            if not rec['relations']:raise ArithmeticError('rank defect without extracted witness')
        return rec
    def relation(self,terms,source):
        m=self.m;sm=max(row[1] for row in terms);ref=min(row[1] for row in terms)
        ans=self.Z
        for c,s,a,b in terms:
            p=m.moment_problem__spatial_coefficient(source(s),a,b)
            ans+=F(c)*m.moment_problem__beta_scale(s,ref)*self.w**(3*(sm-s)//4)*(1-self.w)**((sm-s)//2)*p
        return ans

def independent_check(m,mo,kernels,result,out,smax):
    """Full spectral-order Riccati, not the derivative-graded construction."""
    start=time.monotonic();C=mo.C;w,X,Y=C.gens();Z=mo.Z;D=m.moment_problem__D;poch=m.moment_problem__poch
    # Direct component formulas n=-(t+7)/(t+3), eps=-4/(t+3), cX=1,
    # d2=-1/((t+1)(t+3)), evaluated/differentiated at t=-1/3.
    nn,ee,cc,dd=F(-5,2),F(-3,2),F(1),F(-9,16)
    dn,de,dc,dd2=F(9,16),F(9,16),F(0),F(135,128)
    L=(-ee+nn*w)/2;dL=(-de+dn*w)/2
    V=-cc*w*X/2+((nn+2)*Y/2+F(1,4)+dd)*w*(1-w)-ee*(1-w)/24
    dV=-dc*w*X/2+(dn*Y/2+dd2)*w*(1-w)-de*(1-w)/24
    if V!=mo.V:raise ArithmeticError('independent base operator mismatch')
    R=[-L/2];Rt=[-dL/2]
    R.append((V-R[0]**2-D(R[0]))/2);Rt.append((dV-2*R[0]*Rt[0]-D(Rt[0]))/2)
    for j in range(1,smax):
        R.append(-(D(R[j])-j*L*R[j]+sum((R[i]*R[j-i] for i in range(j+1)),Z))/2)
        Rt.append(-(D(Rt[j])-j*(dL*R[j]+L*Rt[j])+2*sum((R[i]*Rt[j-i] for i in range(j+1)),Z))/2)
    norm={};tops={};jets={};checks=[]
    for s in range(3,smax+1,2):
        aa,bb=s*ee/2,s*nn/2;da,db=s*de/2,s*dn/2;vals={};dvals={}
        for (i,a,b),v in R[s].to_dict().items():
            bf=poch(aa,i)/poch(bb,i);df=bf*sum((da/(aa+j)-db/(bb+j) for j in range(i)),F(0))
            key=0,a,b;vals[key]=vals.get(key,F(0))+v*bf;dvals[key]=dvals.get(key,F(0))+v*df
        for (i,a,b),v in Rt[s].to_dict().items():
            key=0,a,b;dvals[key]=dvals.get(key,F(0))+v*poch(aa,i)/poch(bb,i)
        k=(s+1)//2;tops[s]=vals[0,k,0];norm[s]=C.from_dict(vals)/tops[s]
        jets[s]=(C.from_dict(dvals)-norm[s]*dvals[0,k,0])/tops[s]
        for drop in range(min(k,3)+1):
            for a in range(k-drop+1):
                b=k-drop-a;pair=m.check_a1__coefficient(out,1,k,a,drop)
                if m.check_a1__value(pair,F(-1,3))!=norm[s].to_dict().get((0,a,b),F(0)):raise ArithmeticError('direct A1 value')
                N,DD=map(m.check_a1__univariate,pair);tv=F(-1,3)
                dv=(N.derivative()(tv)*DD(tv)-N(tv)*DD.derivative()(tv))/DD(tv)**2
                if dv!=jets[s].to_dict().get((0,a,b),F(0)):raise ArithmeticError('direct A1 tangent')
                checks.append([s,a,b])
    def adjoints(first):
        js=[{} for _ in range(smax+1)];js[first]={0:C.constant(F(1,2))}
        for j in range(first,smax):
            acc={r:D(v)-j*L*v for r,v in js[j].items()}
            for r,v in js[j].items():acc[r+1]=acc.get(r+1,Z)+v
            for mm in range(j+1):
                for r,v in js[j-mm].items():acc[r]=acc.get(r,Z)+2*R[mm]*v
            js[j+1]={r:-v/2 for r,v in acc.items() if v}
        answer={}
        for s in norm:
            p=Z
            for r,v in js[s].items():
                for _ in range(r):v=-D(v)+s*L*v
                p+=v
            answer[s]=p
        return answer
    adj={q:adjoints(2*q+1) for q in (1,2,3,4)};ind={};kernel_checks=[]
    def independent_response(s,q,a,b):
        k=(s+1)//2;raw=adj[q][s]*X**a*Y**b
        answer=(raw-norm[s]*m.moment_problem__spatial_coefficient(raw,k,0))/tops[s]
        return m.moment_problem__grade(answer,max(0,k-3),k)
    for kind,(q,a,b) in KINDS.items():
        for s in norm:
            p=independent_response(s,q,a,b)
            if p!=kernels[kind,s]:raise ArithmeticError(f'independent kernel {kind,s}')
            ind[kind,s]=p;kernel_checks.append([kind,s])
    record={'status':'EXACT','method':'full spectral-order Riccati plus its t derivative and potential linearization; no derivative-grade truncation in this check',
        'A1_value_and_tangent_inventory':checks,'A1_value_checks':len(checks),'A1_tangent_checks':len(checks),'kernel_inventory':kernel_checks,'kernel_checks':len(kernel_checks)}
    if result['first_failure']:
        terms=result['first_failure']['relation']['terms'];ind_target=F(0)
        for kind in KINDS:
            if mo.relation(terms,lambda s,kind=kind:ind[kind,s]):raise ArithmeticError('independent identity '+kind)
        for c,s,a,b in terms:
            k=(s+1)//2;drop=k-a-b;field='Ward' if drop>=2 else 'A1'
            N,DD=map(m.check_a1__univariate,m.check_a1__coefficient(out,1,k,a,drop,field));tv=F(-1,3)
            dw=(N.derivative()(tv)*DD(tv)-N(tv)*DD.derivative()(tv))/DD(tv)**2
            ind_target+=F(c)*(dw-jets[s].to_dict().get((0,a,b),F(0)))
        if ind_target!=1:raise ArithmeticError('independent target')
        record['nine_identities_zero']=True;record['target_from_Ward_minus_direct_Riccati']=str(ind_target)
        record['higher_order_kernel_checks']=[]
        for q,a,b in [(3,1,0),(3,0,1),(3,2,0),(4,2,0)]:
            for s in norm:
                p=independent_response(s,q,a,b)
                if p!=mo.response((s+1)//2,q,a,b):raise ArithmeticError('independent higher-order response')
                record['higher_order_kernel_checks'].append([s,q,a,b])
        record['top_three_layer_relations']=[]
        for block in result['top_three_layer_blocks']:
            for relation in block['relations']:
                terms2=relation['terms'];target2=F(0)
                for kind in KINDS:
                    if mo.relation(terms2,lambda s,kind=kind:ind[kind,s]):raise ArithmeticError('independent top-three relation')
                for c,s,a,b in terms2:
                    k=(s+1)//2;drop=k-a-b;field='Ward' if drop>=2 else 'A1'
                    N,DD=map(m.check_a1__univariate,m.check_a1__coefficient(out,1,k,a,drop,field));tv=F(-1,3)
                    dw=(N.derivative()(tv)*DD(tv)-N(tv)*DD.derivative()(tv))/DD(tv)**2
                    target2+=F(c)*(dw-jets[s].to_dict().get((0,a,b),F(0)))
                if target2!=1:raise ArithmeticError('independent top-three target')
                record['top_three_layer_relations'].append({'spins':block['spins'],'terms':len(terms2),'all_nine_zero':True,'independent_target':str(target2)})
    # A tempting but wrong implementation omits the grade-one W4 Euler term.
    s=7;k=4
    naive=F(1,2)*m.moment_problem__binom(F(-1,2),k-3)*V**(k-3)*X/m.moment_problem__normalizer(k)
    omitted=ind['XU4',s]-m.moment_problem__grade(naive,k-3,k)
    if not omitted:raise ArithmeticError('missing-W4-derivative control is vacuous')
    record['missing_grade_one_W4_control']=m.moment_problem__encode(omitted)
    record['elapsed_seconds']=round(time.monotonic()-start,3);dump(out,'w4_linear_independent.json',record)
    print('EXACT independent Riccati:',len(checks),'values and derivatives;',len(kernel_checks),'kernels; all relations and controls pass; seconds',record['elapsed_seconds'],flush=True)
    return record

def run(root,out,smax=19):
    if smax not in (15,17,19):raise ValueError('this certificate supports smax 15, 17 or 19; default 19')
    start=time.monotonic();out.mkdir(parents=True,exist_ok=True);m=load_core(root)
    prepare(m,out);mo=Moment(m,out);kernels={};blocks=[];first=None
    for s in range(3,smax+1,2):
        for kind,(q,a,b) in KINDS.items():kernels[kind,s]=mo.response((s+1)//2,q,a,b)
        spins=list(range(3 if s%4==3 else 5,s+1,4))
        block=mo.block(spins,kernels);blocks.append(block)
        print('EXACT nine functions',spins,'rank',block['rank'],'augmented',block['augmented_rank'],'seconds',round(time.monotonic()-start,3),flush=True)
        if first is None and block['rank']!=block['augmented_rank']:
            layer=[]
            for d in range(4):
                r=mo.block(spins,kernels,lastdrop=d,witness=True);layer.append(r)
                if r['rank']!=r['augmented_rank']:
                    first={'spin':s,'degree':(s+1)//2-d,'drop':d,'prefix_layers':layer,'relation':r['relations'][0]};break
            print('EXACT FIRST FAILURE',first['spin'],'degree',first['degree'],'terms',len(first['relation']['terms']),flush=True)
    enc=m.moment_problem__encode
    baseline_matches=[]
    for s in range(3,min(smax,15)+1,2):
        old=m.moment_problem__kernel((s+1)//2,out,mo.V,mo.expanded)
        for kind,p in old.items():
            if p!=kernels[kind,s]:raise ArithmeticError('old seven-function kernel mismatch')
            baseline_matches.append([kind,s])
    oldkinds={key:val for key,val in KINDS.items() if key not in ('XU4','YU4')}
    baseline=mo.block([3,7,11,15],kernels,kinds=oldkinds,maxdrop=2,witness=True)
    if (baseline['rank'],baseline['augmented_rank'])!=(54,55):raise ArithmeticError('seven-function baseline rank')
    oldterms=baseline['relations'][0]['terms']
    if len(oldterms)!=42:raise ArithmeticError('seven-function relation inventory')
    for kind in oldkinds:
        if mo.relation(oldterms,lambda s,kind=kind:kernels[kind,s]):raise ArithmeticError('old relation no longer holds')
    if not mo.relation(oldterms,lambda s:kernels['XU4',s]):raise ArithmeticError('new-class baseline control is vacuous')
    dump(out,'seven_function_baseline.json',{'status':'EXACT','kernel_matches':baseline_matches,'block':baseline,'new_XU4_breaks_old_relation':True})
    dump(out,'w4_linear_kernels.json',{f'{kind}:{s}':enc(p) for (kind,s),p in kernels.items()})
    dump(out,'w4_linear_targets.json',{str(s):{f'{a},{b}':str(v) for (a,b),v in mo.target(s).items()} for s in range(3,smax+1,2)})
    result={'status':'EXACT','fit':'none','class':KINDS,'spins_through':smax,'blocks':blocks,'first_failure':first,
        'layers':'H0 starts at k-3; HX and HY start at k-2 and include derivative-grade-one terms at k-3',
        'exit':2 if first else 3,'scope':'necessary analytic tangent at regular pillow exact point; fixed A1 through kappa^0; same Beta/IBP functional; nine arbitrary functions; not all scalar spectral orders'}
    if first:
        terms=first['relation']['terms'];controls={};identities={}
        for kind in KINDS:
            p=mo.relation(terms,lambda s,kind=kind:kernels[kind,s])
            if p:raise ArithmeticError('nonzero nine-function kernel '+kind)
            identities[kind]=True;dump(out,f'relation_w4_linear_{kind}.json',{'kind':kind,'terms':terms})
        dump(out,'relation_w4_linear_natural.json',{'kind':'U','terms':first['relation']['natural_terms'],'target':first['relation']['natural_target'],'normalization':'last nonzero coefficient one, no target in defining normalization'})
        tamper=[row[:] for row in terms];tamper[0][0]=str(2*F(tamper[0][0]))
        for kind in KINDS:controls[kind]=enc(mo.relation(tamper,lambda s,kind=kind:kernels[kind,s]))
        if not any(controls.values()):raise ArithmeticError('coefficient tamper did not fire')
        extensions={}
        for q in (2,3,4):
            for degree in range(3):
                for a in range(degree+1):
                    b=degree-a
                    p=mo.relation(terms,lambda s,q=q,a=a,b=b:mo.response((s+1)//2,q,a,b))
                    extensions[f'q{q}_X{a}_Y{b}']={'annihilated':not bool(p),'kernel':enc(p)}
        result['nine_identities_zero']=identities;result['wrong_coefficient_control']=controls;result['higher_order_controls']=extensions
    # Through spin19, separately ask whether the top three layers alone fail.
    result['top_three_layer_blocks']=[mo.block(list(range(start_spin,smax+1,4)),kernels,maxdrop=2,witness=True) for start_spin in (3,5) if start_spin<=smax]
    for block in result['top_three_layer_blocks']:
        for rel in block['relations']:
            controls={}
            for kind in KINDS:
                if mo.relation(rel['terms'],lambda s,kind=kind:kernels[kind,s]):raise ArithmeticError('top-three relation')
                dump(out,f'relation_w4_linear_top3_{kind}.json',{'kind':kind,'terms':rel['terms']})
            dump(out,'relation_w4_linear_top3_natural.json',{'kind':'U','terms':rel['natural_terms'],'target':rel['natural_target'],'normalization':'last nonzero coefficient one'})
            for q,a,b in [(3,0,0),(3,1,0),(3,0,1),(3,2,0),(3,1,1),(3,0,2),(4,2,0),(2,2,0)]:
                p=mo.relation(rel['terms'],lambda s,q=q,a=a,b=b:mo.response((s+1)//2,q,a,b))
                controls[f'q{q}_X{a}_Y{b}']={'annihilated':not bool(p),'kernel':enc(p)}
            result['top_three_higher_order_controls']=controls
            if not all(controls[key]['annihilated'] for key in ['q3_X0_Y0','q3_X1_Y0','q3_X0_Y1','q4_X2_Y0']):raise ArithmeticError('degree-bound control')
            if all(controls[key]['annihilated'] for key in ['q3_X2_Y0','q3_X1_Y1','q3_X0_Y2']):raise ArithmeticError('quadratic W6 did not escape the top-three relation')
    independent=independent_check(m,mo,kernels,result,out,smax)
    result['independent_Riccati_verified']=True
    dump(out,'w4_linear_results.json',result)
    summary={'status':'EXACT','fit':'none','spins_through':smax,'first_failure':{key:first[key] for key in ('spin','degree','drop')} if first else None,
        'first_relation_terms':len(first['relation']['terms']) if first else None,'first_target_Sol1':'1' if first else None,'first_target_Sol2':'1/2' if first else None,
        'independent_A1_values':independent['A1_value_checks'],'independent_A1_tangents':independent['A1_tangent_checks'],'independent_kernels':independent['kernel_checks'],
        'top_three_layer_obstructions':[{'spins':b['spins'],'rank':b['rank'],'augmented_rank':b['augmented_rank'],'terms':len(b['relations'][0]['terms'])} for b in result['top_three_layer_blocks'] if b['relations']],
        'scope':result['scope'],'higher_orders':('At spin19 the top-three relation also excludes all added orders q>=3 with momentum degree <=1; quadratic W6 is not excluded. No bound on an unrestricted momentum-degree ladder is proved.' if any(b['relations'] for b in result['top_three_layer_blocks']) else 'NOT DONE in this shorter run; the spin19 top-three relation is required for the higher-linear-order exclusion.'),
        'operator':'NOT DONE; no scalar or matrix operator constructed','exit':result['exit']}
    dump(out,'w4_linear_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    return result['exit']

def verify_relation(root,out,path):
    """Adapter for the lab's existing {kind,terms} relation-file schema."""
    out.mkdir(parents=True,exist_ok=True);data=json.loads(path.read_text());kind=data.get('kind');terms=data.get('terms')
    if kind not in KINDS or not isinstance(terms,list) or not terms:raise ValueError('unknown kind or empty relation')
    for row in terms:
        if not isinstance(row,list) or len(row)!=4:raise ValueError('relation term must have four fields')
        c,s,a,b=row;F(c)
        if not all(type(v) is int for v in (s,a,b)) or s<3 or s%2!=1 or min(a,b)<0 or not 0<=(s+1)//2-a-b<=3:raise ValueError('invalid spin or layer')
    if len({row[1]%4 for row in terms})!=1:raise ValueError('relation mixes Gamma residue classes')
    m=load_core(root);prepare(m,out);mo=Moment(m,out);q,a,b=KINDS[kind]
    kernel=mo.relation(terms,lambda s:mo.response((s+1)//2,q,a,b))
    dump(out,'relation_verification.json',{'status':'EXACT','source':str(path),'kind':kind,'identity':not bool(kernel),'kernel':m.moment_problem__encode(kernel)})
    print('EXACT relation',kind,'PASS' if not kernel else 'FAIL',flush=True)
    return 0 if not kernel else 1

if __name__=='__main__':
    import argparse
    class UsageParser(argparse.ArgumentParser):
        def error(self, message):
            self.print_usage(sys.stderr)
            self.exit(64, f'{self.prog}: error: {message}\n')
    p=UsageParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--out',type=Path,required=True);p.add_argument('--smax',type=int,default=19);p.add_argument('--relation',type=Path)
    a=p.parse_args();sys.exit(verify_relation(a.root.resolve(),a.out.resolve(),a.relation.resolve()) if a.relation else run(a.root.resolve(),a.out.resolve(),a.smax))
