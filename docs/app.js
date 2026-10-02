'use strict';
const $=id=>document.getElementById(id), AA={ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',MSE:'M'};
const palettes={red:[[.63,.63,.63],[.73,.55,.55],[.75,.47,.47],[.77,.38,.38],[.79,.29,.29],[.82,.21,.21],[.84,.13,.13],[.88,0,0],[1,0,.55]],yellow:[[.63,.63,.63],[.66,.67,.51],[.68,.71,.43],[.70,.73,.36],[.72,.76,.28],[.74,.80,.20],[.76,.83,.12],[.79,.87,0],[1,.8,.16]],blue:[[.63,.63,.63],[.50,.58,.68],[.42,.55,.71],[.35,.52,.73],[.28,.49,.76],[.20,.46,.80],[.12,.43,.83],[0,.38,.87],[.6,0,1]],green:[[.63,.63,.63],[.50,.68,.56],[.42,.71,.53],[.35,.74,.49],[.26,.77,.44],[.19,.80,.41],[.12,.83,.37],[.01,.87,.31],[0,1,.83]]};
const bins=[0,50,60,70,80,90,95,98,100];
const branding=window.APP_BRANDING||{name:'CladeScope',logoUrl:'icon.svg'};
document.title=branding.name;$('brandName').textContent=branding.name;$('brandLogo').src=branding.logoUrl;document.querySelector('link[rel=icon]').href=branding.logoUrl;
let state={seqs:[],tree:null,reference:'',structureSeq:'',chain:'',selectedNode:null,selectedSite:null,clade1:null,clade2:null,scores:[],mapping:new Map(),reverse:new Map()},viewer,model,chainResidues=new Map();
function message(s){$('message').textContent=s}
function parseFasta(s){const seqs=[];for(const chunk of s.replace(/\r/g,'').split('>')){if(!chunk.trim())continue;const lines=chunk.trim().split('\n'),header=lines.shift(),id=header.trim().split(/\s+/)[0],seq=lines.join('').replace(/\s/g,'').toUpperCase().replace(/\./g,'-');if(!/^[A-Z-]+$/.test(seq))throw Error('Alignment contains unsupported characters. Use aligned protein FASTA.');seqs.push({id,seq})}if(!seqs.length)throw Error('No FASTA sequences found.');if(new Set(seqs.map(x=>x.id)).size!==seqs.length)throw Error('FASTA IDs must be unique.');if(seqs.some(x=>x.seq.length!==seqs[0].seq.length))throw Error('Sequences must have the same aligned length.');return seqs}
function parseNewick(text) {
  let i=0,n=0; const s=text.trim();
  function skip(){
    while(i<s.length){
      if(/\s/.test(s[i])){i++;continue}
      if(s[i]==='['){let depth=1;i++;while(i<s.length&&depth){if(s[i]==='[')depth++;if(s[i]===']')depth--;i++}if(depth)throw Error('Unclosed Newick comment.');continue}
      break;
    }
  }
  function label(){skip();if(s[i]==="'"){i++;let value='';while(i<s.length){if(s[i]==="'"){i++;if(s[i]==="'"){value+="'";i++;continue}return value}value+=s[i++]}throw Error('Unclosed quoted Newick name.')}
    let value='';while(i<s.length&&!/[,:();]/.test(s[i])){if(s[i]==='['){skip();continue}value+=s[i++]}return value.trim();
  }
  function node(){skip();const o={id:'node'+n++,name:'',children:[],length:null};if(s[i]==='('){i++;o.children.push(node());skip();while(s[i]===','){i++;o.children.push(node());skip()}if(s[i++]!==')')throw Error('Invalid Newick parentheses.')}
    o.name=label();skip();if(s[i]===':'){i++;skip();let value='';while(i<s.length&&!/[,();]/.test(s[i])){if(s[i]==='['){skip();continue}value+=s[i++]}o.length=Number(value);if(!value.trim()||!Number.isFinite(o.length)||o.length<0)throw Error('Invalid branch length.')}
    o.ids=o.children.length?o.children.flatMap(x=>x.ids):[o.name];if(!o.children.length&&!o.name)throw Error('Unnamed tree leaf.');return o;
  }
  const root=node();skip();if(s[i]===';')i++;skip();if(i!==s.length)throw Error('Unexpected text in Newick.');
  let order=0;function annotate(o){o.children.forEach(annotate);o.order=o.children.length?Math.min(...o.children.map(x=>x.order)):order++}annotate(root);return root;
}
function refreshTree(node){node.children.forEach(refreshTree);node.ids=node.children.length?node.children.flatMap(x=>x.ids):[node.name];node.order=node.children.length?Math.min(...node.children.map(x=>x.order)):node.order;return node}
function ladderizeTree(node,mode){node.children.forEach(x=>ladderizeTree(x,mode));node.children.sort((a,b)=>(mode==='ascending'?a.ids.length-b.ids.length:mode==='descending'?b.ids.length-a.ids.length:0)||a.order-b.order);return refreshTree(node)}
let rootSerial=0;
function rerootTree(root,selectedId){
  // A unary root stem has no effect on leaf relationships.
  while(root.children.length===1)root=root.children[0];
  const graph=new Map();let selected,parent;
  function collect(o,p){graph.set(o,[]);if(o.id===selectedId){selected=o;parent=p}o.children.forEach(x=>{collect(x,o);graph.get(o).push({node:x,length:x.length});graph.get(x).push({node:o,length:x.length})})}collect(root,null);
  if(!selected||!parent)throw Error('Select a non-root node to choose a root branch.');
  const edge=selected.length,half=edge===null?null:edge/2;
  function orient(o,from,length){const children=graph.get(o).filter(e=>e.node!==from).map(e=>orient(e.node,o,e.length));
    if(o===root&&children.length===1){const child=children[0];child.length=length===null||child.length===null?null:length+child.length;return child}
    return {...o,length,children};
  }
  return refreshTree({id:'root'+(++rootSerial),name:'New root',length:null,children:[orient(selected,parent,half),orient(parent,selected,half)]});
}
function validateTree(tree,seqs){const ids=new Set(seqs.map(x=>x.id));if(new Set(tree.ids).size!==tree.ids.length)throw Error('Tree leaf IDs must be unique.');const absent=tree.ids.filter(x=>!ids.has(x));if(absent.length)throw Error('Tree leaves missing from FASTA: '+absent.slice(0,6).join(', '));return seqs.filter(x=>!tree.ids.includes(x.id)).length}
function calculateSites(seqs,reference,ids,cutoff){const ref=seqs.find(x=>x.id===reference);const group=seqs.filter(x=>ids.includes(x.id));return [...ref.seq].map((aa,c)=>{const counts={};group.forEach(x=>counts[x.seq[c]]=(counts[x.seq[c]]||0)+1);const valid=group.length-(counts['-']||0)-(counts.X||0);const count=aa==='X'||aa==='-'?0:counts[aa]||0;const denom=valid/group.length<cutoff?group.length:valid;return {aa,counts,total:group.length,valid,count,denom,score:aa==='-'||!valid||!group.length?null:100*count/denom}})}
function color(score){if(score===null)return '#bfbf94';let b=0;for(let i=0;i<bins.length;i++)if(score>=bins[i])b=i;return '#'+palettes[$('palette').value][b].map(v=>Math.round(v*255).toString(16).padStart(2,'0')).join('')}
function cladeSites(seqs,ids,cutoff){
  const group=seqs.filter(x=>ids.includes(x.id)),length=seqs[0]?.seq.length||0;
  return Array.from({length},(_,column)=>{
    const counts={};for(const seq of group){const aa=seq.seq[column];counts[aa]=(counts[aa]||0)+1}
    const valid=group.length-(counts['-']||0)-(counts.X||0);
    const denom=group.length&&valid/group.length<cutoff?group.length:valid;
    const entries=Object.entries(counts).filter(([aa])=>aa!=='-'&&aa!=='X');
    const count=entries.length?Math.max(...entries.map(([,n])=>n)):0;
    const dominant=entries.filter(([,n])=>n===count).map(([aa])=>aa).sort();
    return {counts,total:group.length,valid,denom,count,dominant,score:!valid||!group.length?null:100*count/denom};
  });
}
function sharedConservation(site1,site2){
  if(site1.score===null||site2.score===null)return {score:null,dominant:[]};
  const aas=[...new Set([...Object.keys(site1.counts),...Object.keys(site2.counts)])].filter(aa=>aa!=='-'&&aa!=='X');
  const frequencies=aas.map(aa=>[aa,Math.min((site1.counts[aa]||0)/site1.denom,(site2.counts[aa]||0)/site2.denom)]);
  const fraction=Math.max(0,...frequencies.map(([,value])=>value));
  return {score:100*fraction,dominant:fraction?frequencies.filter(([,value])=>value===fraction).map(([aa])=>aa).sort():[]};
}
function comparisonColor(site1,site2,shared,threshold=50){
  if(site1.score===null||site2.score===null)return '#bfbf94';
  if(site1.score<threshold&&site2.score<threshold)return '#555555';
  return '#'+[site2.score/100,shared.score/100,site1.score/100].map(v=>Math.round(Math.max(0,Math.min(1,v))*255).toString(16).padStart(2,'0')).join('');
}
function calculateComparison(seqs,reference,ids1,ids2,cutoff,threshold=50){
  const ref=seqs.find(x=>x.id===reference),one=cladeSites(seqs,ids1,cutoff),two=cladeSites(seqs,ids2,cutoff);
  return one.map((site,i)=>{const shared=sharedConservation(site,two[i]);return {aa:ref.seq[i],clade1:site,clade2:two[i],shared:shared.score,sharedDominant:shared.dominant,comparisonColor:comparisonColor(site,two[i],shared,threshold)}});
}
function isComparison(){return $('mode').value==='comparison'}
function siteColor(site){return isComparison()?site.comparisonColor||'#bfbf94':color(site.score)}
function sameMembers(a,b){return !!a&&a.ids.length===b.ids.length&&a.ids.every(id=>b.ids.includes(id))}
function makeClade(node){return {ids:[...node.ids],label:node.name||'Clade ('+node.ids.length+' sequences)'}}
function comparisonGroups(ids1,ids2,includeOverlap=false){
  const overlap=ids1.filter(id=>ids2.includes(id));
  let one=[...ids1],two=[...ids2];
  if(!includeOverlap&&overlap.length){
    if(overlap.length===ids1.length&&ids1.length<ids2.length)two=ids2.filter(id=>!ids1.includes(id));
    else if(overlap.length===ids2.length&&ids2.length<ids1.length)one=ids1.filter(id=>!ids2.includes(id));
    else{one=ids1.filter(id=>!overlap.includes(id));two=ids2.filter(id=>!overlap.includes(id))}
  }
  return {one,two,overlap,excluded1:ids1.length-one.length,excluded2:ids2.length-two.length};
}
function effectiveGroups(){return comparisonGroups(state.clade1?.ids||[],state.clade2?.ids||[],$('includeOverlap').checked)}
function selectedIds(){if(!isComparison())return state.selectedNode.ids;const groups=effectiveGroups();return [...new Set([...groups.one,...groups.two])]}
function assignClade(which){state['clade'+which]=makeClade(state.selectedNode);message('');update()}
function options(id,values,value){$(id).replaceChildren(...values.map(v=>{const o=document.createElement('option');o.value=v;o.textContent=v||'(blank chain)';return o}));$(id).value=value}
function loadStructure(data,format){const scratch=$3Dmol.createViewer(document.createElement('div'));const parsed=scratch.addModel(data,format);const atoms=parsed.selectedAtoms({});if(!atoms.length)throw Error('No atoms could be read from this structure.');const chains=new Map();for(const a of atoms){if(a.atom!=='CA'||!AA[a.resn])continue;const ch=a.chain||'',key=String(a.resi)+':'+(a.icode||'');if(!chains.has(ch))chains.set(ch,new Map());if(!chains.get(ch).has(key))chains.get(ch).set(key,{key,resi:a.resi,icode:a.icode||'',aa:AA[a.resn]})}if(!chains.size)throw Error('No protein C-alpha atoms found.');chainResidues=new Map([...chains].map(([k,v])=>[k,[...v.values()]]));viewer.removeAllModels();model=viewer.addModel(data,format);state.chain=[...chains.keys()][0];options('chain',[...chains.keys()],state.chain);viewer.zoomTo();}
// Global alignment to observed residues; gaps account for unresolved structure residues.
function alignMap(full,observed){const n=full.length,m=observed.length;if(n*m>8000000)throw Error('Chain mapping is limited to 8 million sequence-pair cells in this prototype.');const cols=m+1,trace=new Uint8Array((n+1)*cols);let prev=new Int32Array(cols),cur=new Int32Array(cols);for(let j=1;j<=m;j++){prev[j]=-3*j;trace[j]=2}for(let i=1;i<=n;i++){cur[0]=-3*i;trace[i*cols]=1;for(let j=1;j<=m;j++){const d=prev[j-1]+(full[i-1]===observed[j-1]?3:-2),u=prev[j]-3,l=cur[j-1]-3;cur[j]=Math.max(d,u,l);trace[i*cols+j]=cur[j]===d?0:cur[j]===u?1:2}const t=prev;prev=cur;cur=t}let i=n,j=m;const pairs=[];while(i||j){const t=trace[i*cols+j];if(i&&j&&t===0){pairs.push([--i,--j])}else if(i&&(t===1||!j))i--;else j--}return pairs.reverse()}
function mapStructure(){const seq=state.seqs.find(x=>x.id===state.structureSeq),res=chainResidues.get(state.chain)||[];const columns=[];[...seq.seq].forEach((aa,i)=>{if(aa!=='-')columns.push(i)});const full=columns.map(i=>seq.seq[i]).join('');const pairs=alignMap(full,res.map(r=>r.aa).join(''));state.mapping=new Map();state.reverse=new Map();let matches=0;for(const [i,j] of pairs){state.mapping.set(columns[i],res[j]);state.reverse.set(res[j].key,columns[i]);if(full[i]===res[j].aa)matches++}const pct=pairs.length?100*matches/pairs.length:0;$('mapping').textContent=`Chain ${state.chain||'(blank)'} · ${pairs.length}/${res.length} observed residues mapped · ${pct.toFixed(1)}% sequence agreement`;if(pct<90)message('Mapping has less than 90% sequence agreement. Verify the structure-associated sequence and chain before interpreting residue colors.');return pct}
function renderTree(){const root=state.tree;let row=0;const nodes=[],edges=[];const maxDepth=(o,d=0)=>Math.max(d,...o.children.map(x=>maxDepth(x,d+1)));const depth=maxDepth(root)||1;function visit(o,d){o.x=24+d*(180/depth);o.y=o.children.length?0:30+row++*36;o.children.forEach(ch=>{visit(ch,d+1);edges.push([o,ch])});if(o.children.length)o.y=o.children.reduce((a,b)=>a+b.y,0)/o.children.length;nodes.push(o)}visit(root,0);const treeWidth=Math.max(560,...nodes.map(o=>o.x+28+(o.children.length?(o.name||'Clade')+' · '+o.ids.length:o.name).length*8+50));const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('width',treeWidth);svg.setAttribute('height',Math.max(180,row*36+35));svg.setAttribute('viewBox',`0 0 ${treeWidth} ${Math.max(180,row*36+35)}`);svg.style.minWidth=treeWidth+'px';function el(tag,attrs){const e=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);return e}for(const [a,b] of edges)svg.append(el('path',{d:`M ${a.x} ${a.y} V ${b.y} H ${b.x}`,class:'branch'}));for(const o of nodes){const g=el('g',{class:'node'+(o===state.selectedNode?' chosen':'')+(sameMembers(state.clade1,o)?' clade1':'')+(sameMembers(state.clade2,o)?' clade2':''),tabindex:0,role:'button','aria-label':'Compare '+(o.name||'clade')+' '+o.ids.length+' sequences'});g.append(el('circle',{cx:o.x,cy:o.y,r:o.children.length?6:4}));const t=el('text',{x:o.x+12,y:o.y+5});t.textContent=o.children.length?(o.name||'Clade')+' · '+o.ids.length:o.name;if(sameMembers(state.clade1,o)||sameMembers(state.clade2,o))t.textContent+=' ['+[sameMembers(state.clade1,o)?'1':null,sameMembers(state.clade2,o)?'2':null].filter(Boolean).join(', ')+']';g.append(t);const click=()=>{state.selectedNode=o;update()};g.onclick=click;g.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();click()}};svg.append(g)}$('tree').replaceChildren(svg);$('treecount').textContent=root.ids.length+' leaves · rooted cladogram';$('reroot').disabled=state.selectedNode===root;}
function renderMsa(){const memberships=effectiveGroups(),frag=document.createDocumentFragment(),group=new Set(selectedIds()),ref=state.seqs.find(x=>x.id===state.reference);const rows=[ref,...state.seqs.filter(x=>x!==ref)];const ruler=document.createElement('div');ruler.className='msa-row ruler';const title=document.createElement('div');title.className='msa-name';title.textContent='Alignment column';ruler.append(title);for(let c=0;c<ref.seq.length;c++){const cell=document.createElement('span');cell.className='cell';cell.textContent=(c+1)%10===0?c+1:'';ruler.append(cell)}frag.append(ruler);for(const seq of rows){if($('filter').checked&&seq!==ref&&!group.has(seq.id))continue;const row=document.createElement('div');row.className='msa-row'+(group.has(seq.id)?' active':'')+(seq===ref?' reference':'');const name=document.createElement('span');name.className='msa-name'+(memberships.one.includes(seq.id)?' clade1':'')+(memberships.two.includes(seq.id)?' clade2':'');name.textContent=(seq===ref?'★ ':'')+seq.id;name.title=seq.id;row.append(name);for(let c=0;c<seq.seq.length;c++){const b=document.createElement('button');b.className='cell'+(state.selectedSite===c?' selected':'');b.textContent=seq.seq[c];b.style.background=siteColor(state.scores[c]);const rgb=b.style.background.match(/\d+/g);if(rgb&&rgb.length>=3&&(.299*Number(rgb[0])+.587*Number(rgb[1])+.114*Number(rgb[2]))>180){b.style.color='#142c40';b.style.textShadow='none'}b.title=`${seq.id} · column ${c+1} · ${seq.seq[c]}`;b.dataset.column=c;b.onclick=()=>selectSite(c);row.append(b)}frag.append(row)}$('msa').replaceChildren(frag)}
function residueKey(a){return (a.chain||'')+'|'+String(a.resi)+':'+(a.icode||'')}
function isWater(a){return ['HOH','WAT','H2O','DOD','SOL','TIP3','TIP3P'].includes((a.resn||'').trim().toUpperCase())}
function structureSelections(){
  const polymerKeys=new Set();for(const [chain,residues] of chainResidues)for(const r of residues)polymerKeys.add(chain+'|'+r.key);
  const chainVisible=a=>$('showOtherChains').checked||(a.chain||'')===state.chain;
  return {protein:{predicate:a=>polymerKeys.has(residueKey(a))&&chainVisible(a)},ligand:{predicate:a=>!polymerKeys.has(residueKey(a))&&!isWater(a)&&chainVisible(a)},water:{predicate:a=>isWater(a)&&chainVisible(a)}};
}
function renderStructure(){
  if(!model)return;
  viewer.removeAllShapes();viewer.setStyle({},{});model.setClickable({},false);
  const selections=structureSelections(),scores=new Map();for(const [c,r] of state.mapping)scores.set(r.key,siteColor(state.scores[c]));
  const style={colorfunc:a=>(a.chain||'')===state.chain?scores.get(String(a.resi)+':'+(a.icode||''))||'#bfbf94':'#bfbf94'};
  viewer.setStyle(selections.protein,$('representation').value==='stick'?{stick:style}:{cartoon:style});
  if($('showLigands').checked)viewer.setStyle(selections.ligand,{stick:{colorscheme:'Jmol',radius:.16},sphere:{colorscheme:'Jmol',scale:.25}});
  if($('showWater').checked)viewer.setStyle(selections.water,{sphere:{colorscheme:'Jmol',scale:.25}});
  model.setClickable(selections.protein,true,a=>{if((a.chain||'')!==state.chain)return;const c=state.reverse.get(String(a.resi)+':'+(a.icode||''));if(c!==undefined)selectSite(c,true)});
  if(state.selectedSite!==null){const r=state.mapping.get(state.selectedSite);if(r)viewer.addStyle({chain:state.chain,resi:r.resi,...(r.icode?{icode:r.icode}:{})},{stick:{color:'#ffbe3e',radius:.25},sphere:{color:'#ffbe3e',scale:.25}})}viewer.render();
}
function selectSite(c,scroll=false){state.selectedSite=c;document.querySelectorAll('.cell[data-column]').forEach(e=>e.classList.toggle('selected',Number(e.dataset.column)===c));if(scroll)$('msa').scrollLeft=Math.max(0,c*23-200);renderStructure();renderDetail()}
function renderDetail(){
  const c=state.selectedSite;if(c===null){$('detail').textContent='Click an alignment cell or a structure residue.';return}
  const s=state.scores[c],r=state.mapping.get(c),ref=state.seqs.find(x=>x.id===state.reference);
  const fmt=x=>x===null?'No data':x.toFixed(1)+'%';
  const items=[['Alignment column',c+1],['Reference residue',s.aa==='-'?'Gap':ref.seq.slice(0,c+1).replace(/-/g,'').length+' · '+s.aa],['Structure',r?(state.chain||'(blank)')+':'+r.resi+r.icode:'Not mapped']];
  if(isComparison())items.push(['Clade 1 conservation',(s.clade1.dominant.join('/')||'—')+' · '+fmt(s.clade1.score)],['Clade 2 conservation',(s.clade2.dominant.join('/')||'—')+' · '+fmt(s.clade2.score)],['Shared conservation',(s.sharedDominant.join('/')||'—')+' · '+fmt(s.shared)]);else items.push(['Reference identity',fmt(s.score)],['Score denominator',s.count+' / '+s.denom]);
  const box=document.createElement('div');box.className='sitegrid';items.forEach(([name,value])=>{const d=document.createElement('div'),l=document.createElement('label'),v=document.createElement('strong');l.textContent=name;v.textContent=value;d.append(l,v);box.append(d)});
  const table=document.createElement('table'),head=document.createElement('tr');
  const headers=isComparison()?['Residue','Clade 1 count','Clade 1 %','Clade 2 count','Clade 2 %']:['Residue','Count','% of selected group'];
  headers.forEach(t=>{const th=document.createElement('th');th.textContent=t;head.append(th)});table.append(head);
  const counts=isComparison()?[s.clade1,s.clade2]:[s],aas=[...new Set(counts.flatMap(x=>Object.keys(x.counts)))];
  aas.sort((a,b)=>counts.reduce((n,x)=>n+(x.counts[b]||0)-(x.counts[a]||0),0));
  aas.forEach(aa=>{const tr=document.createElement('tr');const values=[aa,...counts.flatMap(x=>[x.counts[aa]||0,x.total?(100*(x.counts[aa]||0)/x.total).toFixed(1)+'%':'—'])];values.forEach(t=>{const td=document.createElement('td');td.textContent=t;tr.append(td)});table.append(tr)});
  const note=document.createElement('p');note.className='hint';note.textContent=isComparison()?'Blue/red = each clade’s most common amino-acid frequency. Green = strongest amino-acid frequency shared by both clades. Neither clade reaches the selected minimum conservation = dark gray. Gaps/X use the coverage cutoff; missing data in either group = beige.':'Below the non-gap cutoff, the score uses the full group size; otherwise it excludes gaps and X.';
  $('detail').replaceChildren(box,table,note);
}
function renderLegend(){
  $('legendTitle').textContent=isComparison()?'Two-clade conservation':'Reference identity';
  if(isComparison()){
    const examples=[['Only clade 1 conserved','#0000ff'],['Only clade 2 conserved','#ff0000'],['Different AAs conserved','#ff00ff'],['Same AA in both clades','#ffffff'],['Neither reaches '+$('threshold').value+'%','#555555'],['Unmapped/ No data','#bfbf94']];
    $('legend').replaceChildren(...examples.map(([label,col])=>{const d=document.createElement('div');d.style.borderColor=col;d.textContent=label;return d}));
    $('legend').title='R = most common AA frequency in clade 2; B = most common AA frequency in clade 1; G = max over amino acids of min(frequency in clade 1, frequency in clade 2). Blue/red stripes indicate the clade channels; actual site colors include the other channels according to their frequencies.';
  }else{$('legend').title='';$('legend').replaceChildren(...bins.map((b,i)=>{const d=document.createElement('div');d.style.borderColor=color(b);d.textContent=i===0?'0–<50':i===8?'100%':b+'–<'+bins[i+1];return d}));const missing=document.createElement('div');missing.style.borderColor='#bfbf94';missing.textContent='Unmapped/ No data';$('legend').append(missing)}
}
function update(){
  if(!state.seqs.length)return;const cutoff=Number($('cutoff').value);if(!Number.isFinite(cutoff)||cutoff<0||cutoff>1){message('Non-gap fraction must be between 0 and 1.');return}
  const threshold=isComparison()?Number($('threshold').value):50;if(!Number.isFinite(threshold)||threshold<0||threshold>100){message('Minimum conservation must be between 0 and 100%.');return}
  $('thresholdControl').hidden=false;$('threshold').disabled=!isComparison();$('threshold').title=isComparison()?'Color sites where at least one clade meets this threshold.':'Used only in Comparison mode.';
  const groups=effectiveGroups();
  state.scores=isComparison()?calculateComparison(state.seqs,state.reference,groups.one,groups.two,cutoff,threshold):calculateSites(state.seqs,state.reference,state.selectedNode.ids,cutoff);
  $('all').textContent=isComparison()?'Select root':'All sequences';$('palette').disabled=isComparison();$('palette').title=isComparison()?'Comparison uses fixed RGB channels. Switch to Identity to use this palette.':'';
  $('cladeAssignment').hidden=!isComparison();$('filterLabel').textContent=isComparison()?'Selected groups + reference only':'Selected group + reference only';
  $('treeHint').textContent=isComparison()?'Select a node, then assign it to clade 1 or 2. Rerooting preserves assigned sequence groups.':'Select a tree node to compare the reference against its clade.';
  $('comparison').textContent=isComparison()?'Clade 1: '+groups.one.length+' · Clade 2: '+groups.two.length:(state.selectedNode===state.tree?'All tree sequences':state.selectedNode.name||'Selected clade')+' · n='+state.selectedNode.ids.length;
  for(const n of [1,2]){
    const clade=state['clade'+n],ids=n===1?groups.one:groups.two,excluded=n===1?groups.excluded1:groups.excluded2;
    $('clade'+n+'Summary').textContent='Clade '+n+': '+(clade?clade.label+' · n='+ids.length+(excluded?' of '+clade.ids.length+' selected':''):'Not assigned');
    $('clade'+n+'Summary').title='Counted species: '+(ids.join(', ')||'None');
  }
  $('overlap').textContent=groups.overlap.length?($('includeOverlap').checked?groups.overlap.length+' overlapping species contribute to both scores.':groups.overlap.length+' overlapping species: '+[groups.excluded1?'excluded from clade 1 ('+groups.excluded1+')':null,groups.excluded2?'excluded from clade 2 ('+groups.excluded2+')':null].filter(Boolean).join('; ')+'.'+(!groups.one.length||!groups.two.length?' One or both comparison groups are empty; choose different clades or include overlapping species.':'')):'';
  renderTree();renderMsa();renderStructure();renderDetail();renderLegend();
}
function renderFileInformation(names,metadata){
  const stats=document.createElement('div');stats.className='statgrid';
  const fields=[['Alignment',names.alignment||'Loaded alignment'],['Alignment length',state.seqs[0].seq.length.toLocaleString()+' columns'],['Alignment taxa',state.seqs.length.toLocaleString()],['Tree',names.tree||'Loaded tree'],['Tree taxa',state.tree.ids.length.toLocaleString()],['Structure',names.structure||'Loaded structure']];
  for(const [label,value] of fields){const d=document.createElement('div'),l=document.createElement('label'),v=document.createElement('span');l.textContent=label;v.textContent=value;d.append(l,v);stats.append(d)}$('fileStats').replaceChildren(stats);
  const content=document.createDocumentFragment();
  const table=document.createElement('table');table.className='metadatafields';
  const info=[['Format',metadata.format],['Entry ID',metadata.entry],['Title',metadata.title],['Classification',metadata.classification],['Structure authors',metadata.authors],['Deposition date',metadata.deposited],['First release',metadata.released],['Latest revision',metadata.revision],['Experimental method',metadata.method],['Resolution',metadata.resolution?metadata.resolution+' Å':''],['Source organisms',metadata.organisms],['Models',metadata.modelCount],['Water residues',metadata.waterCount],['Citation',metadata.citation],['DOI',metadata.doi],['Software',metadata.software]];
  for(const [label,value] of info){if(value===undefined||value==='')continue;const tr=document.createElement('tr'),th=document.createElement('th'),td=document.createElement('td');th.textContent=label;td.textContent=value;tr.append(th,td);table.append(tr)}content.append(table);
  function section(title,headers,rows){const h=document.createElement('h3');h.textContent=title;content.append(h);if(!rows.length){const p=document.createElement('p');p.textContent='None recorded in this file.';content.append(p);return}const t=document.createElement('table'),head=document.createElement('tr');for(const label of headers){const th=document.createElement('th');th.textContent=label;head.append(th)}t.append(head);for(const values of rows){const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value||'—';tr.append(td)}t.append(tr)}content.append(t)}
  section('Polymer chains',['Chain','Molecule / complex component','Observed residues'],(metadata.chains||[]).map(x=>[x.chain,x.description,x.observedResidues]));
  section('Ligands and non-polymer components',['ID','Name','Chains','Residues'],(metadata.ligands||[]).map(x=>[x.id,x.name,x.chains,x.copies]));
  const note=document.createElement('p');note.className='hint';note.textContent=metadata.error?'Some metadata could not be read: '+metadata.error:'Metadata is read from the loaded file. Chain and residue counts describe observed coordinates, not missing residues or inferred assemblies. Unavailable fields are omitted.';content.append(note);$('structureMetadata').replaceChildren(content);
}
function setDataset(seqs,tree,data,format,fileNames={}){const extra=validateTree(tree,seqs);let metadata;try{metadata=StructureMetadata.parse(data,format)}catch(error){metadata={format:format==='cif'?'mmCIF':'PDB',error:error.message,chains:[],ligands:[]}}loadStructure(data,format);state.structureFile=fileNames.structure||('structure.'+(format==='cif'?'cif':'pdb'));state.seqs=seqs;state.tree=tree;state.inputTree=structuredClone(tree);$('ladderize').value='input';state.selectedNode=tree;state.clade1=tree.children.length>1?makeClade(tree.children[0]):null;state.clade2=tree.children.length>1?makeClade(tree.children[1]):null;state.reference=seqs[0].id;state.structureSeq=seqs[0].id;state.selectedSite=null;options('reference',seqs.map(x=>x.id),state.reference);options('structureSeq',seqs.map(x=>x.id),state.structureSeq);message(extra?extra+' alignment sequences are absent from the tree and excluded from comparisons.':'');mapStructure();renderFileInformation(fileNames,metadata);update()}
async function fetchDemoFile(path){const response=await fetch(path);if(!response.ok)throw Error('Could not load demo file: '+path);return response.text()}
async function demo(){
  try{
    message('Loading GFP demo…');const config=JSON.parse(await fetchDemoFile('demo.json'));
    const [alignment,tree,structure]=await Promise.all([fetchDemoFile(config.alignment),fetchDemoFile(config.tree),fetchDemoFile(config.structure)]);
    const seqs=parseFasta(alignment),phylogeny=parseNewick(tree);
    for(const id of [config.referenceSequence,config.structureSequence])if(!seqs.some(x=>x.id===id))throw Error('Demo sequence is not in alignment: '+id);
    setDataset(seqs,phylogeny,structure,config.format,{alignment:config.alignment,tree:config.tree,structure:config.structure});
    if(!chainResidues.has(config.chain))throw Error('Demo chain is not in structure: '+config.chain);
    state.reference=config.referenceSequence;state.structureSeq=config.structureSequence;state.chain=config.chain;
    $('reference').value=state.reference;$('structureSeq').value=state.structureSeq;$('chain').value=state.chain;
    mapStructure();update();
  }catch(e){message(e.message)}
}

function exportCategory(site){
  if(!site)return 'no_data';
  if(!isComparison())return site.score===null?'no_data':site.score<50?'unconserved':'identical';
  const a=site.clade1,b=site.clade2,t=Number($('threshold').value);
  if(a.score===null||b.score===null)return 'no_data';
  const one=a.score>=t,two=b.score>=t;
  if(!one&&!two)return 'unconserved';
  if(one&&!two)return 'clade_1_conserved';
  if(two&&!one)return 'clade_2_conserved';
  return a.dominant.some(aa=>b.dominant.includes(aa))?'identical':'diverging';
}
function exportResidueColors(){
  const mapped=new Map();for(const [c,r] of state.mapping)mapped.set(r.key,{color:siteColor(state.scores[c]),group:exportCategory(state.scores[c])});
  const rows=[];for(const [chain,residues] of chainResidues)for(const r of residues)rows.push({chain,resi:r.resi,icode:r.icode||'',...(chain===state.chain?mapped.get(r.key):null)||{color:'#bfbf94',group:'no_data'}});
  return rows;
}
function exportScriptName(program){return 'cladescope_'+(state.structureFile||'structure').replace(/\.(?:pdb|cif|mmcif)$/i,'').toLowerCase().replace(/[^a-z0-9_-]+/g,'_')+'_'+program+'.py'}
function colorScript(program){
  const rows=exportResidueColors();if(!rows.length)throw Error('Load a structure before exporting.');
  // One Python string per JSON line avoids long lines and safely escapes filenames.
  const payload=JSON.stringify({file:state.structureFile||'structure.cif',residues:rows},null,2).split('\n').map(line=>'    '+JSON.stringify(line+'\n')).join('\n');
  const header='# CladeScope residue colors: '+(isComparison()?'Comparison':'Identity')+' mode\n# Preload the same PDB/CIF file before running this script. Nothing is loaded here.\n# With multiple structures open, set TARGET below to the desired object/model.\n# Selection names use underscores in place of spaces.\n# Groups follow the current conservation threshold and clade consensus amino acids.\n# Identity mode: identical = identity >= 50%; unconserved = identity < 50%.\n# The temporary gold inspection highlight is excluded.\n';
  const common='import json\ndata = json.loads(\n'+payload+'\n)\ngroup_names = ["clade_1_conserved", "clade_2_conserved", "diverging", "identical", "unconserved", "no_data"]\n';
  if(program==='pymol')return header+'# Run with File > Run Script, or run /path/to/'+exportScriptName(program)+'\n'+`from pymol import cmd
cmd.hide("everything", "solvent")
TARGET = ""  # Optional: name of the already loaded PyMOL object.
objects = cmd.get_object_list("all")
if not TARGET:
    if len(objects) != 1:
        raise ValueError("Preload one structure, or set TARGET to its PyMOL object name.")
    TARGET = objects[0]
if TARGET not in objects:
    raise ValueError("TARGET is not an already loaded PyMOL object.")
`+common+`records = {(r["chain"], str(r["resi"]) + r["icode"]): r for r in data["residues"]}
by_color = {}
groups = {name: [] for name in group_names}
for atom in cmd.get_model(TARGET).atom:
    record = records.get((atom.chain, atom.resi))
    if record:
        by_color.setdefault(record["color"], []).append(atom.index)
        groups[record["group"]].append(atom.index)
def atom_selection(indices):
    return "model " + json.dumps(TARGET) + " and index " + "+".join(map(str, indices)) if indices else "none"
for i, (hex_color, indices) in enumerate(by_color.items()):
    name = "cladescope_rgb_" + str(i)
    cmd.set_color(name, [int(hex_color[j:j+2], 16) / 255.0 for j in (1, 3, 5)])
    cmd.color(name, atom_selection(indices))
for name, indices in groups.items():
    cmd.select(name, atom_selection(indices), enable=0)
print("CladeScope colors applied to " + TARGET + ". Selections: " + ", ".join(group_names))
`;
  if(program==='chimerax')return header+'# Run with open /path/to/'+exportScriptName(program)+'\n'+`from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
run(session, "hide solvent atoms")
TARGET = ""  # Optional: model ID without #, for example "1" or "1.1".
models = session.models.list(type=AtomicStructure)
if TARGET:
    models = [m for m in models if m.id_string == TARGET]
if len(models) != 1:
    raise ValueError("Preload one structure, or set TARGET to its ChimeraX model ID.")
model = models[0]
`+common+`records = {(r["chain"], int(r["resi"]), r["icode"]): r for r in data["residues"]}
groups = {name: [] for name in group_names}
for residue in model.residues:
    record = records.get((residue.chain_id, residue.number, residue.insertion_code.strip()))
    if record:
        rgba = [int(record["color"][j:j+2], 16) for j in (1, 3, 5)] + [255]
        residue.atoms.colors = rgba
        residue.ribbon_color = rgba
        groups[record["group"]].append(residue)
# Build frozen named selections from atom selections, avoiding chain-name parsing.
previous = [(m, m.atoms.selected.copy()) for m in session.models.list(type=AtomicStructure)]
try:
    for name, residues in groups.items():
        run(session, "select clear")
        for residue in residues:
            residue.atoms.selected = True
        run(session, "name frozen " + name + " sel")
finally:
    run(session, "select clear")
    for m, selected in previous:
        m.atoms.selected = selected
session.logger.info("CladeScope colors applied. Selections: " + ", ".join(group_names))
`;
  throw Error('Unsupported export program.');
}
function downloadColorScript(){
  try{const program=$('exportProgram').value,script=colorScript(program),url=URL.createObjectURL(new Blob([script],{type:'text/x-python;charset=utf-8'}));const link=document.createElement('a');link.href=url;link.download=exportScriptName(program);document.body.append(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);message('Script exported. Preload '+state.structureFile+' in '+(program==='pymol'?'PyMOL':'ChimeraX')+', then run the script.')}catch(e){message(e.message)}
}

$('exportScript').onclick=downloadColorScript;
$('load').onclick=async()=>{try{const files=['fasta','newick','pdb'].map(id=>$(id).files[0]);if(files.some(x=>!x))throw Error('Select all three files first.');if(files.some(x=>x.size>25*1024*1024))throw Error('Each file must be below 25 MB for this prototype.');const [fa,nw,pdb]=await Promise.all(files.map(f=>f.text()));const seqs=parseFasta(fa);if(seqs.length*seqs[0].seq.length>250000)throw Error('This prototype supports up to 250,000 alignment cells.');const tree=parseNewick(nw);setDataset(seqs,tree,pdb,/\.mm?cif$|\.cif$/i.test(files[2].name)?'cif':'pdb',{alignment:files[0].name,tree:files[1].name,structure:files[2].name});}catch(e){message(e.message)}};
$('remap').onclick=()=>{try{state.structureSeq=$('structureSeq').value;state.chain=$('chain').value;message('');mapStructure();update()}catch(e){message(e.message)}};
$('reference').onchange=()=>{state.reference=$('reference').value;update()};['mode','palette','cutoff','threshold','includeOverlap','filter','representation','showLigands','showWater','showOtherChains'].forEach(id=>$(id).onchange=update);$('all').onclick=()=>{state.selectedNode=state.tree;update()};$('fit').onclick=()=>{const selections=structureSelections();viewer.zoomTo({predicate:a=>selections.protein.predicate(a)||($('showLigands').checked&&selections.ligand.predicate(a))||($('showWater').checked&&selections.water.predicate(a))});viewer.render()};$('demo').onclick=demo;
$('assign1').onclick=()=>assignClade(1);$('assign2').onclick=()=>assignClade(2);
$('ladderize').onchange=()=>{ladderizeTree(state.tree,$('ladderize').value);update()};
$('reroot').onclick=()=>{try{state.tree=rerootTree(state.tree,state.selectedNode.id);ladderizeTree(state.tree,$('ladderize').value);state.selectedNode=state.tree;message(isComparison()?'Root moved; assigned comparison groups retained.':'Root moved. Identity group reset to all tree sequences.');update()}catch(e){message(e.message)}};
$('resetRoot').onclick=()=>{state.tree=structuredClone(state.inputTree);ladderizeTree(state.tree,$('ladderize').value);state.selectedNode=state.tree;message(isComparison()?'Input root restored; assigned comparison groups retained.':'Input root restored. Identity group reset to all tree sequences.');update()};
try{viewer=$3Dmol.createViewer($('viewer'),{backgroundColor:'#0d1d2c',antialias:true});new ResizeObserver(()=>viewer.resize()).observe($('viewer'));demo()}catch(e){message('The 3D viewer could not initialize. A browser with WebGL support is required. '+e.message)}
