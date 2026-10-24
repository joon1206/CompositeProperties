"""Fetch verifiable Crossref metadata for a manually curated DOI set.

Run only to refresh metadata; review notes are maintained separately. No
paper-level scientific claims are inferred from bibliographic metadata.
"""
from pathlib import Path
import csv,html,json,re,time,urllib.parse,urllib.request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research'/'literature'
DOIS=[
 '10.1016/S0266-3538(02)00017-9',
 '10.1177/002199830003400605',
 '10.1016/j.compositesa.2017.02.015',
 '10.1016/S1359-6454(98)00089-5',
 '10.1016/j.matpr.2017.06.113',
 '10.1016/j.compositesb.2013.10.049',
 '10.1016/S0266-3538(00)00073-7',
 '10.1016/j.jmrt.2026.01.011',
 '10.1016/j.tca.2016.08.014',
 '10.1177/002199836800200308',
 '10.1016/0001-6160(73)90064-3',
 '10.1016/0022-5096(63)90060-7',
 '10.1088/0370-1301/69/8/305',
 '10.1002/pen.760160512',
 '10.1016/j.compositesb.2013.02.002',
 '10.1016/0015-0568(80)90012-3',
 '10.1016/S0142-9418(01)00066-6',
 '10.1016/j.msea.2006.10.060',
 '10.1016/j.commatsci.2012.05.020',
 '10.1002/pc.750050303',
 '10.1137/120888971',
 '10.1016/j.apm.2020.04.003',
 '10.1016/0010-4361(69)90005-6',
 '10.1098/rspa.1957.0133',
 '10.1016/0022-5096(65)90015-3',
 '10.1061/(ASCE)AS.1943-5525.0000450',
 '10.1039/C3SM52959B',
 '10.1016/j.compscitech.2007.07.015',
 '10.1016/j.ijadhadh.2017.06.006',
 '10.1177/09673911231171479',
 '10.3795/KSME-A.2015.39.2.129',
 '10.1039/C1JM12405F',
]

# Brief notes only where the primary publisher abstract was inspected.
NOTES={
 DOIS[0]:('layered media','exact flat-interface thermoelastic continuity','Closest exact layered-media method; prevents methodology novelty claim'),
 DOIS[1]:('laminates','3D sublaminate homogenization; FE verification','Closest elastic homogenization precedent'),
 DOIS[2]:('glass fiber / polyimide / silica nanoparticles','3D unit-cell micromechanics with interphase','Closest hierarchical material precedent; different resin and architecture'),
 DOIS[3]:('silica-filled epoxy encapsulant','axisymmetric unit-cell and viscoelastic matrix','Particle arrangement affects local stress more than CTE in reported setup'),
 DOIS[4]:('silica-filled epoxy','experiment and empirical comparison','Measured 70 wt% case; composition and matrix differ from our assumed inputs'),
 DOIS[5]:('silica nanofiber / epoxy','experimental morphology comparison','Challenges spherical filler assumption'),
 DOIS[6]:('glass fiber / polypropylene angle-ply','experimental CTE and mean-field micromechanics','Creep can change CTE and challenges elastic-only modeling'),
 DOIS[7]:('woven glass / FR4','temperature-dependent measurement and multiscale model','Application-specific prior work; not a direct input dataset'),
 DOIS[8]:('glass fiber / silica bead / epoxy prepreg','cure kinetics','Shows relevant material architecture and unmodeled cure'),
 DOIS[9]:('general composites','thermoelastic energy bounds','Classical CTE baseline'),
 DOIS[11]:('multiphase solids','variational elastic bounds','Source of stiffness bounds'),
 DOIS[12]:('spherical particles','elastic/thermoelastic effective media','Classical particulate comparison'),
 DOIS[14]:('carbon / epoxy cross-ply','FE and measured CTE','Relevant validation method; different materials'),
 DOIS[17]:('fiber-reinforced polymer','viscoelastic unit-cell residual stress','Unmodeled cure and viscoelastic effects'),
 DOIS[18]:('random glass fiber / epoxy','FE and experiment','Different architecture; potential validation methodology'),
 DOIS[19]:('glass / epoxy laminate','experimental laminate theory test','Historical experimental comparison'),
 DOIS[20]:('rough-interface laminates','asymptotic homogenization','Flat perfect interface assumption limitation'),
 DOIS[21]:('periodic composites','thermoelastic Fourier eigenstrain','Established numerical homogenization alternative'),
 DOIS[26]:('monodisperse and polydisperse hard spheres','packing simulation','The 0.70 nested filler fraction requires a packing explanation'),
 DOIS[27]:('silica mineral / epoxy','experimental particle-shape comparison','Challenges spherical particle closure'),
 DOIS[28]:('silica nanoparticle / epoxy','thermomechanical and moisture experiments','Reports high-loading deviation of Kerner and Halpin-Tsai'),
 DOIS[29]:('colloidal silica / epoxy','dispersion experiment','High loading and dispersion context'),
 DOIS[30]:('silica microparticle / epoxy','experimental modulus and CTE','Potential validation source; only abstract inspected'),
 DOIS[31]:('bimodal silica / epoxy','experimental filler morphology','Challenges monodisperse packing assumption'),
}

FINDINGS={
 DOIS[0]:'Derives exact effective thermoelastic layer formulas under interface continuity',
 DOIS[1]:'Reports 3D finite-element agreement with exact sublaminate formulas',
 DOIS[3]:'Reported effective CTE less sensitive than local thermal stress to studied viscoelasticity and particle arrangement',
 DOIS[4]:'At 70 wt% silica reported E increase about 51-55% and CTE reduction about 25% versus its baseline',
 DOIS[5]:'At 8.74 vol% silica nanofiber reported about 40% CTE reduction versus neat epoxy',
 DOIS[6]:'Reported measured angle-ply CTE mismatch with simple elastic prediction and analysis including creep',
 DOIS[7]:'Measured temperature-dependent FR4 thermal properties and used multiscale modeling',
 DOIS[9]:'Derived energy-principle bounds and special exact composite CTE cases',
 DOIS[11]:'Derived variational elastic bounds for multiphase materials',
 DOIS[12]:'Derived elastic and thermoelastic effective-medium relations for spherical grains',
 DOIS[26]:'Computed characteristic monodisperse random close packing near volume fraction 0.64',
 DOIS[27]:'Reported particle-shape effects on silica/epoxy mechanical properties',
 DOIS[28]:'Reported slight Kerner/Halpin-Tsai deviation at high filler weight fraction',
 DOIS[31]:'Reported bimodal silica/epoxy CTE reduction with increasing nanoparticle content',
}
LIMITATIONS={
 DOIS[0]:'Flat perfectly bonded linear layers; no particle-level morphology',
 DOIS[1]:'Layer homogenization; does not identify constituent material properties',
 DOIS[2]:'Polyimide and silica nanoparticles rather than the present illustrative resin and loading',
 DOIS[3]:'Epoxy encapsulant particle cell, not a cross-ply laminate',
 DOIS[4]:'Reported filler load by weight; exact resin and packing not matched to present input',
 DOIS[5]:'Nanofiber morphology differs from assumed silica spheres',
 DOIS[6]:'Polypropylene and angle-ply, not the present epoxy cross-ply',
 DOIS[7]:'Woven FR4 and temperature dependence differ from constant-property unidirectional model',
 DOIS[26]:'Ideal hard-sphere packing cannot establish manufacturability of this resin mixture',
 DOIS[27]:'Mineral particle shapes differ from ideal spherical silica',
 DOIS[28]:'Tested nanoparticle epoxy systems, not this three-phase laminate',
}


def fetch(doi):
    url='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
    req=urllib.request.Request(url,headers={'User-Agent':'composite-research/0.1 (bibliographic verification)'})
    with urllib.request.urlopen(req,timeout=20) as f:return json.load(f)['message']


def clean(s):return re.sub(r'<[^>]+>','',html.unescape(s or '')).replace('\n',' ').strip()


def tex(s):
    return str(s).replace('&',r'\&').replace('%',r'\%').replace('#',r'\#').replace('_',r'\_')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[];bib=[];failed=[]
    for i,doi in enumerate(DOIS):
        try:m=fetch(doi)
        except Exception as e:failed.append((doi,str(e)));continue
        date=m.get('published',m.get('issued',{})).get('date-parts',[[None]])[0]
        authors=[(' '.join(filter(None,[a.get('given',''),a.get('family','')]))).strip() for a in m.get('author',[])]
        title=clean(m.get('title',[''])[0]);venue=clean(m.get('container-title',[''])[0] if m.get('container-title') else '')
        access='publisher abstract inspected' if doi in NOTES else 'metadata only'
        mat,method,relation=NOTES.get(doi,('','','Bibliographic context; scientific content not assessed'))
        rows.append({'doi':m.get('DOI',doi),'title':title,'authors':'; '.join(authors),
                     'year':date[0],'venue':venue,'url':m.get('URL','https://doi.org/'+doi),
                     'access':access,'materials':mat,'assumptions_methods':method,
                     'validated_findings':FINDINGS.get(doi,''),'limitations':LIMITATIONS.get(doi,''),'relationship':relation,
                     'novelty_effect':'challenges method novelty' if i in (0,1,2) else ''})
        key='ref'+str(i+1)
        fields={'author':' and '.join(authors),'title':title,'year':date[0],'journal':venue,'doi':m.get('DOI',doi),'url':m.get('URL','https://doi.org/'+doi)}
        bib.append('@article{'+key+',\n'+''.join(f'  {k} = {{{tex(v)}}},\n' for k,v in fields.items())+'}\n')
        time.sleep(.08)
    with (OUT/'papers.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'references.bib').write_text('\n'.join(bib),encoding='utf-8')
    (OUT/'metadata_fetch_log.json').write_text(json.dumps({'verified_count':len(rows),'failed':failed,'source':'Crossref REST API','date':'2026-10-24'},indent=2),encoding='utf-8')
    print('verified',len(rows),'failed',failed)


if __name__=='__main__':main()
