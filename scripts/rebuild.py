"""Rebuild empirical figures/tables from saved results; no experiments or bootstrap."""
from pathlib import Path
import argparse, os, shutil, subprocess, sys
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'build')
    parser.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
    parser.add_argument('--gs',default=os.environ.get('GS','gs'))
    args=parser.parse_args();out=args.output.resolve()
    # A rebuild may replace its own output, never the source results or paper figures.
    if out==ROOT or any(out==ROOT/x or ROOT/x in out.parents for x in ['results','figures','src','data','scripts','experiments','docs']):
        parser.error('Output must be build/ or a separate directory.')
    out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['REBUILD_OUTPUT']=str(out)
    for script in ['build_tables.py','plot_scaling.py','plot_network.py','plot_sachs.py']:
        subprocess.run([sys.executable,str(ROOT/'scripts'/script)],check=True,env=env)
    subprocess.run([args.rscript,str(ROOT/'scripts/plot_components.R')],check=True,env=env)
    figure=out/'figures/bma_confirmatory.pdf';embedded=figure.with_name('component_embedded.pdf')
    subprocess.run([args.gs,'-q','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite','-dCompatibilityLevel=1.7',
        '-dEmbedAllFonts=true','-dSubsetFonts=true','-dPDFSETTINGS=/prepress',
        '-sOutputFile='+str(embedded),'-c','<</NeverEmbed []>> setdistillerparams','-f',str(figure)],check=True)
    embedded.replace(figure)
    for path in (out/'figures').glob('*.pdf'):
        for page in PdfReader(path).pages:
            for font in page['/Resources'].get('/Font',{}).values():
                font=font.get_object()
                for child in font.get('/DescendantFonts',[font]):
                    child=child.get_object();descriptor=child.get('/FontDescriptor')
                    assert descriptor is not None,(path,child)
                    assert any(k in descriptor.get_object() for k in ['/FontFile','/FontFile2','/FontFile3']),path
    print('Rebuilt four empirical figures with embedded fonts in '+str(out))

if __name__=='__main__':main()
