
#!/usr/bin/env python3
"""Build Master SNP Matrix V2 from Snippy *_snps.tab files."""
import argparse, glob, os
from pathlib import Path
import pandas as pd

ANN=["FTYPE","STRAND","NT_POS","AA_POS","LOCUS_TAG","GENE","PRODUCT","EFFECT"]
REQ=["CHROM","POS","TYPE","REF","ALT"]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("-i","--input",default=".")
    ap.add_argument("-o","--output",default=".")
    args=ap.parse_args()
    files=sorted(glob.glob(os.path.join(args.input,"*_snps.tab")))
    if not files: raise SystemExit("No *_snps.tab files found.")
    master={};samples=[];mismatch=[];per={};types={}
    for i,f in enumerate(files,1):
        s=os.path.basename(f).replace("_snps.tab","");samples.append(s);per[s]=0
        print(f"[{i}/{len(files)}] {s}")
        df=pd.read_csv(f,sep="\t",dtype=str).fillna("")
        miss=[c for c in REQ if c not in df.columns]
        if miss: raise ValueError(f"{f} missing {miss}")
        for _,r in df.iterrows():
            key=(r["CHROM"],r["POS"],r["TYPE"],r["REF"],r["ALT"])
            if key not in master:
                d={"Variant_ID": f"{r['CHROM']}_{r['POS']}_{r['TYPE']}_{r['REF']}_{r['ALT']}"}
                for c in REQ:d[c]=r[c]
                for a in ANN:d[a]=r[a] if a in df.columns else ""
                master[key]=d
                types[r["TYPE"]]=types.get(r["TYPE"],0)+1
            else:
                for a in ANN:
                    new=r[a] if a in df.columns else ""
                    old=master[key].get(a,"")
                    if old=="" and new!="":
                        master[key][a]=new
                    elif old!="" and new!="" and old!=new:
                        mismatch.append(f"{s}\t{master[key]['Variant_ID']}\t{a}\tOLD={old}\tNEW={new}")
            master[key][s]=f"{r['REF']}>{r['ALT']}"
            per[s]+=1
    rows=[];brows=[]
    meta=["Variant_ID"]+REQ+ANN
    for rec in master.values():
        row={k:rec.get(k,"") for k in meta}; brow=row.copy()
        for s in samples:
            if s in rec: row[s]=rec[s];brow[s]=1
            else: row[s]="0";brow[s]=0
        rows.append(row);brows.append(brow)
    m=pd.DataFrame(rows);b=pd.DataFrame(brows)
    m["POS"]=m["POS"].astype(int);b["POS"]=b["POS"].astype(int)
    m=m.sort_values(["CHROM","POS","TYPE","REF","ALT"])
    b=b.sort_values(["CHROM","POS","TYPE","REF","ALT"])
    out=Path(args.output);out.mkdir(parents=True, exist_ok=True)
    m.to_csv(out/"Master_SNP_Matrix.csv",index=False)
    b.to_csv(out/"Master_SNP_Binary.csv",index=False)
    (out/"Annotation_Mismatch_Report.txt").write_text("\n".join(mismatch) if mismatch else "No annotation mismatch detected.\n")
    summ=["Master SNP Matrix V2","",f"Samples: {len(samples)}",f"Unique variants: {len(m)}",""]
    for t,n in sorted(types.items()): summ.append(f"{t}: {n}")
    summ.append(f"\nAnnotation conflicts: {len(mismatch)}\n")
    summ.append("Per sample:")
    for s,n in per.items(): summ.append(f"{s}\t{n}")
    (out/"Variant_Summary.txt").write_text("\n".join(summ))
    print("Done.")
if __name__=="__main__":
    main()
