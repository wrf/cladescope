const fs=require('fs'),vm=require('vm'),assert=require('assert');let src=fs.readFileSync(require('path').join(__dirname,'../docs/app.js'),'utf8');src=src.slice(0,src.indexOf("$('load').onclick"));const controls={mode:{value:'comparison'},palette:{value:'blue'},showLigands:{checked:true},showWater:{checked:false},showOtherChains:{checked:true}};const ctx={window:{},document:{getElementById:id=>controls[id]||{},querySelector:()=>({})},console,structuredClone,fs,__dirname,require};vm.createContext(ctx);vm.runInContext(src,ctx);vm.runInContext(`
function check(v,m){if(!v)throw Error(m)}
const rooted=parseNewick("[&R] ((A:1[tip],B:2)AB:3,(C:4,D:5)CD:6)Root:0;");check(rooted.name==='Root'&&rooted.children[0].name==='AB','Preserve input root');
check(parseNewick("[&R] ('a[b]':1,'it''s':2)root;").ids.join('|')==="a[b]|it's",'Quoted names and comments');
function distances(root){const adj=new Map(),leaves=new Map();function walk(o,p){adj.set(o,[]);if(!o.children.length)leaves.set(o.name,o);o.children.forEach(x=>{walk(x,o);adj.get(o).push([x,x.length||0]);adj.get(x).push([o,x.length||0])})}walk(root,null);const out={};for(const [name,start] of leaves){function dfs(o,p,d){if(!o.children.length&&o!==start)out[[name,o.name].sort().join('|')]=d;for(const [x,l] of adj.get(o))if(x!==p)dfs(x,o,d+l)}dfs(start,null,0)}return out}
const unary=parseNewick('((A:1,B:2)AB:3)Stem;');const unaryReroot=rerootTree(unary,unary.children[0].children[0].id);check(unaryReroot.ids.slice().sort().join(',')==='A,B','Unary root stem');
const before=distances(rooted);let rerooted=rerootTree(rooted,rooted.children[1].children[0].id);check(rerooted.ids.slice().sort().join(',')==='A,B,C,D','Reroot leaves');check(JSON.stringify(Object.entries(before).sort())===JSON.stringify(Object.entries(distances(rerooted)).sort()),'Reroot preserves distances');
rerooted=rerootTree(rerooted,rerooted.children[1].id);check(JSON.stringify(Object.entries(before).sort())===JSON.stringify(Object.entries(distances(rerooted)).sort()),'Repeated reroot distances');
const uneven=parseNewick('(A,(B,C,D)big)root;');ladderizeTree(uneven,'descending');check(uneven.children[0].name==='big','Descending ladderize');ladderizeTree(uneven,'ascending');check(uneven.children[0].name==='A','Ascending ladderize');ladderizeTree(uneven,'input');check(uneven.ids.join(',')==='A,B,C,D','Restore input ordering');
const seqs=parseFasta('>a\\nAK-X\\n>b\\nAR-X\\n>c\\nA--X'),sites=calculateSites(seqs,'a',['a','b','c'],.5);check(sites[0].score===100&&sites[1].score===50&&sites[2].score===null&&sites[3].score===null,'Identity scores');check(calculateSites(seqs,'a',['a','b','c'],.8)[1].score===100/3,'Coverage rule');check(alignMap('ABCDE','ABDE').map(p=>p.join(':')).join(',')==='0:0,1:1,3:2,4:3','Missing residue mapping');
chainResidues=new Map([['A',[{key:'1:',resi:1}]],['B',[{key:'2:',resi:2}]]]);state.chain='A';
const a={chain:'A',resi:1},b={chain:'B',resi:2},ligA={chain:'A',resi:50},ligB={chain:'B',resi:50};
$('showOtherChains').checked=false;let selections=structureSelections();check(selections.protein.predicate(a)&&!selections.protein.predicate(b),'Hide other chains');check(selections.ligand.predicate(ligA)&&!selections.ligand.predicate(ligB)&&!selections.ligand.predicate(a),'Ligand classification');
const waterA={chain:'A',resi:51,resn:'HOH'},waterB={chain:'B',resi:51,resn:'WAT'};check(selections.water.predicate(waterA)&&!selections.water.predicate(waterB)&&!selections.ligand.predicate(waterA),'Independent water classification');
$('showOtherChains').checked=true;selections=structureSelections();check(selections.protein.predicate(b)&&selections.ligand.predicate(ligB),'Show other chains');
const calls=[];viewer={removeAllShapes(){},setStyle(sel,style){calls.push([sel,style])},render(){}};model={setClickable(){}};
$('representation').value='cartoon';$('showLigands').checked=false;renderStructure();check(calls.length===2,'Hidden ligands receive no style');calls.length=0;$('showLigands').checked=true;renderStructure();check(calls.length===3&&calls[2][0].predicate(ligA)&&!calls[2][0].predicate(waterA),'Visible ligands exclude water');calls.length=0;$('showWater').checked=true;renderStructure();check(calls.length===4&&calls[3][0].predicate(waterA),'Show water independently');calls.length=0;$('showLigands').checked=false;renderStructure();check(calls.length===3&&calls[2][0].predicate(waterA),'Water remains visible when ligands hidden');

const demoSeqs=parseFasta(fs.readFileSync(require('path').join(__dirname,'../docs/demo-alignment.fasta'),'utf8'));
const demoTree=parseNewick(fs.readFileSync(require('path').join(__dirname,'../docs/demo-tree.nwk'),'utf8'));
check(demoTree.ids.length===13&&demoSeqs.length===13,'GFP demo has 13 taxa');validateTree(demoTree,demoSeqs);
check(demoTree.children[0].ids.length===3&&demoTree.children[1].ids.length===10,'GFP original root split');

const cif=fs.readFileSync(require('path').join(__dirname,'../docs/1EMA.cif'),'utf8');
const residues=cif.split('\\n').filter(l=>l.startsWith('ATOM ')).map(l=>l.trim().split(/\\s+/)).filter(v=>v[3]==='CA'&&v[18]==='A'&&AA[v[17]]).map(v=>AA[v[17]]);
const native=demoSeqs.find(x=>x.id==='Aequorea_victoria|AAA27721.1__WT-GFP').seq.replace(/-/g,'');const demoPairs=alignMap(native,residues.join(''));const agreement=demoPairs.filter(([i,j])=>native[i]===residues[j]).length/demoPairs.length;
check(residues.length>200&&agreement>.95,'GFP structure/reference agreement');console.log('GFP demo:',demoTree.ids.length,'taxa;',residues.length,'observed alpha carbons;', (100*agreement).toFixed(1)+'% sequence agreement');


const toy=parseFasta('>reference\\nA-G\\n>one\\nLLG\\n>two\\nVLG\\n>three\\nVLG\\n>four\\nVLG');
const compare=calculateComparison(toy,'reference',['one'],['two','three','four'],.5);
check(compare[0].comparisonColor==='#ff00ff','Fixed L vs V is magenta despite unequal sizes and absent reference AA');
check(compare[0].clade1.dominant[0]==='L'&&compare[0].clade2.dominant[0]==='V','Each clade has its own consensus');
check(compare[1].comparisonColor==='#ffffff','Reference gaps do not mask conserved comparison columns');
check(compare[2].comparisonColor==='#ffffff','Same AA fixed in both is white');
const otherReference=calculateComparison(toy,'two',['one'],['two','three','four'],.5);check(compare.every((x,i)=>x.comparisonColor===otherReference[i].comparisonColor),'Comparison is independent of selected reference');
const ties=parseFasta('>r\\nA\\n>a\\nL\\n>b\\nV\\n>c\\nL\\n>d\\nV');
const tied=calculateComparison(ties,'r',['a','b'],['c','d'],.5)[0];check(tied.comparisonColor==='#808080'&&tied.clade1.dominant.join('/')==='L/V','Exactly 50/50 meets the inclusive default threshold and ties are retained');
const diffuse=parseFasta('>r\\nA\\n>a\\nL\\n>b\\nV\\n>c\\nG\\n>d\\nL\\n>e\\nV\\n>f\\nG');check(calculateComparison(diffuse,'r',['a','b','c'],['d','e','f'],.5)[0].comparisonColor==='#555555','Diffuse clades are gray');
check(calculateComparison(toy,'reference',[],['two'],.5)[0].comparisonColor==='#bfbf94','Missing group is beige');
const gaps=parseFasta('>r\\nA\\n>a\\nL\\n>b\\n-\\n>c\\nX\\n>d\\nL');const gapScore=calculateComparison(gaps,'r',['a','b','c'],['d'],.5)[0];check(gapScore.clade1.score===100/3&&gapScore.clade1.dominant[0]==='L','Coverage cutoff penalizes gaps/X');
const partial=calculateComparison(ties,'r',['a'],['c','d'],.5)[0];check(partial.comparisonColor==='#8080ff','Fixed L vs half L/half V uses shared green');
const swapped=calculateComparison(ties,'r',['c','d'],['a'],.5)[0];check(swapped.comparisonColor==='#ff8080','Swapping groups swaps red and blue');
check(calculateComparison(ties,'r',['a','b'],['c','d'],.5,51)[0].comparisonColor==='#555555','Raising threshold hides 50 percent sites');
check(calculateComparison(ties,'r',['a'],['c','d'],.5,100)[0].comparisonColor==='#8080ff','One fully conserved clade passes 100 percent threshold');
const ninety=parseFasta('>r\\nA\\n>a\\nL\\n>b\\nL\\n>c\\nV\\n>d\\nG');check(calculateComparison(ninety,'r',['a','b','c'],['b','c','d'],.5,100)[0].comparisonColor==='#555555','Neither fully conserved hides at 100 percent');
check(calculateComparison(toy,'reference',['one'],['two','three','four'],.5,100)[0].comparisonColor==='#ff00ff','Fixed different AA passes 100 percent');
const nested=comparisonGroups(['a','b','c'],['a','b','c','d','e','f','g','h','i','j']);check(nested.one.length===3&&nested.two.length===7&&nested.excluded2===3,'Nested 3 vs 10 becomes 3 vs 7');
const nestedSeqs=[{id:'ref',seq:'A'},...nested.one.map(id=>({id,seq:'L'})),...nested.two.map(id=>({id,seq:'V'}))];check(calculateComparison(nestedSeqs,'ref',nested.one,nested.two,.5)[0].comparisonColor==='#ff00ff','Exclusive parent species drive the score, not full parent membership');
const reverseNested=comparisonGroups(['a','b','c','d'],['a','b']);check(reverseNested.one.join(',')==='c,d'&&reverseNested.two.join(',')==='a,b','Reverse nesting preserves smaller clade');
check(comparisonGroups(['a','b'],['a','b','c'],true).two.length===3,'Include toggle restores full parent');
const identical=comparisonGroups(['a','b'],['a','b']);check(!identical.one.length&&!identical.two.length,'Identical clades leave no exclusive species');
const partialOverlap=comparisonGroups(['a','b'],['b','c']);check(partialOverlap.one.join(',')==='a'&&partialOverlap.two.join(',')==='c','Partial overlap excluded from both');
const originalIds=['a','b'];comparisonGroups(originalIds,['a','b','c']);check(originalIds.join(',')==='a,b','Assignments are not mutated');
state.clade1={ids:['one','two']};state.clade2={ids:['two','three']};check(selectedIds().sort().join(',')==='one,three','Comparison filtering uses non-overlapping effective groups');
$('mode').value='identity';state.selectedNode={ids:['four']};check(selectedIds().join(',')==='four','Identity filter uses clicked node');
check(siteColor({score:100})===color(100),'Identity color follows original palette');$('mode').value='comparison';
const frozen=makeClade(rooted.children[0]);const altered=rerootTree(rooted,rooted.children[1].children[0].id);check(frozen.ids.join(',')==='A,B'&&altered.ids.length===4,'Assignments are independent of rerooted node objects');
`,ctx);console.log('Passed: independent clade consensuses, shared conservation, reference invariance, threshold boundaries, missing data, mode-specific filtering, frozen clades; rooted/commented Newick, quoted names, reroot topology and distances, repeated rerooting, ladderizing and original order, chain/ligand visibility, identity scoring, residue mapping.');
