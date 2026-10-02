const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert');
class Element {
  constructor(tag){this.tagName=tag;this.children=[];this.style={};this.dataset={};this.attributes={};this.value='';this.checked=false;this.className='';this.textContent='';this.classList={toggle(){}}}
  append(...nodes){this.children.push(...nodes)}
  replaceChildren(...nodes){this.children=[...nodes]}
  setAttribute(name,value){this.attributes[name]=value;if(name==='class')this.className=value}
}
const html=fs.readFileSync(path.join(__dirname,'../docs/index.html'),'utf8');assert.ok(!/id="showWater"\s+checked/.test(html));
const elements={};const get=id=>elements[id]||(elements[id]=new Element('div'));
get('mode').value='comparison';get('palette').value='blue';get('cutoff').value='0.5';get('threshold').value='50';get('representation').value='cartoon';
const document={getElementById:get,querySelector:()=>new Element('link'),querySelectorAll:()=>[],createElement:tag=>new Element(tag),createElementNS:(ns,tag)=>new Element(tag),createDocumentFragment:()=>new Element('fragment')};
const ctx={document,window:{},console,structuredClone};vm.createContext(ctx);
let src=fs.readFileSync(path.join(__dirname,'../docs/app.js'),'utf8');src=src.slice(0,src.indexOf('try{viewer=$3Dmol.createViewer'));
vm.runInContext(src,ctx);
vm.runInContext(`
state.seqs=parseFasta('>r\\nAK\\n>a\\nAK\\n>b\\nAR\\n>c\\nAR');state.tree=parseNewick('((r,a)A,(b,c)B)Root;');state.inputTree=structuredClone(state.tree);state.reference='r';state.selectedNode=state.tree;state.clade1=makeClade(state.tree.children[0]);state.clade2=makeClade(state.tree.children[1]);update();
`,ctx);
assert.equal(get('palette').disabled,true);assert.equal(get('threshold').disabled,false);assert.equal(get('thresholdControl').hidden,false);assert.equal(get('cladeAssignment').hidden,false);assert.equal(get('legend').children.length,6);
const nodes=get('tree').children[0].children.filter(e=>e.tagName==='g');assert.equal(nodes.filter(e=>e.className.includes('clade1')).length,1);assert.equal(nodes.filter(e=>e.className.includes('clade2')).length,1);
const leaf=nodes.find(n=>n.children[1].textContent==='b');leaf.onclick();get('assign1').onclick();assert.ok(get('clade1Summary').textContent.includes('n=1'));assert.ok(get('overlap').textContent.includes('excluded from clade 2 (1)'));
get('includeOverlap').checked=true;get('includeOverlap').onchange();assert.ok(get('clade2Summary').textContent.includes('n=2'));assert.ok(get('overlap').textContent.includes('contribute to both'));get('includeOverlap').checked=false;get('includeOverlap').onchange();
get('filter').checked=true;get('filter').onchange();let rows=get('msa').children[0].children.filter(e=>e.className.startsWith('msa-row')&&!e.className.includes('ruler'));assert.equal(rows.length,3); // reference + b,c
get('mode').value='identity';get('mode').onchange();assert.equal(get('palette').disabled,false);assert.equal(get('threshold').disabled,true);assert.equal(get('thresholdControl').hidden,false);assert.equal(get('cladeAssignment').hidden,true);assert.equal(get('legend').children.length,10);rows=get('msa').children[0].children.filter(e=>e.className.startsWith('msa-row')&&!e.className.includes('ruler'));assert.equal(rows.length,2); // reference + clicked b
get('mode').value='comparison';get('mode').onchange();vm.runInContext('selectSite(0)',ctx);assert.equal(get('detail').children[0].children.length,6);assert.equal(get('detail').children[1].children[0].children.length,5);
get('reroot').onclick();assert.ok(get('clade1Summary').textContent.includes('n=1'));assert.ok(get('clade2Summary').textContent.includes('n=1 of 2 selected'));
vm.runInContext("renderFileInformation({alignment:'test.fa',tree:'test.nwk',structure:'test.cif'},{format:'mmCIF',entry:'TEST',chains:[{chain:'A',description:'Alpha',observedResidues:2}],ligands:[],waterCount:0,modelCount:1})",ctx);assert.equal(get('fileStats').children[0].children.length,6);assert.ok(get('structureMetadata').children.length);
console.log('Passed interface checks: assigned node colors, click/assign actions, overlap notice, comparison/identity controls and filters, detail table, reroot membership.');
