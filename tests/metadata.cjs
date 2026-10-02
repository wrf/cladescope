const fs=require('fs'),path=require('path'),assert=require('assert');
const api=require('../docs/metadata.js');
const gfp=api.parse(fs.readFileSync(path.join(__dirname,'../docs/1EMA.cif'),'utf8'),'cif');
assert.equal(gfp.entry,'1EMA');assert.ok(gfp.authors.includes('Ormo'));assert.equal(gfp.resolution,'1.9');assert.equal(gfp.chains[0].chain,'A');assert.ok(gfp.waterCount>0);assert.ok(gfp.revision.includes('version'));assert.equal(gfp.ligands.length,0); // chromophore is a modified polymer residue
const text=`data_test
_entry.id TEST
_struct.title
;A multiline
structure title
;
loop_
_audit_author.name
'Author One'
"Author Two"
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Protein alpha'
2 polymer 'Protein beta'
3 non-polymer 'ATP'
4 water 'Water'
loop_
_chem_comp.id
_chem_comp.name
ATP 'ADENOSINE TRIPHOSPHATE'
loop_
_atom_site.group_PDB
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
ATOM A 1 ALA 1 1
ATOM A 1 ALA 1 1
ATOM B 1 GLY 2 1
HETATM A 9 ATP 3 1
HETATM A 9 ATP 3 1
HETATM B 9 ATP 3 1
HETATM A 10 HOH 4 1
`;
const cif=api.parse(text,'cif');assert.equal(cif.authors,'Author One; Author Two');assert.ok(cif.title.includes('multiline'));assert.equal(cif.chains.length,2);assert.equal(cif.chains[1].description,'Protein beta');assert.equal(cif.ligands[0].copies,2);assert.equal(cif.ligands[0].name,'ADENOSINE TRIPHOSPHATE');assert.equal(cif.waterCount,1);
function atom(tag,name,chain,residue){const line=Array(80).fill(' ');function put(i,s){[...s].forEach((c,n)=>line[i+n]=c)}put(0,tag.padEnd(6));put(12,' CA ');put(17,name.padEnd(3));put(21,chain);put(22,String(residue).padStart(4));return line.join('')}
const pdb=['HEADER    TEST COMPLEX                            01-JAN-00   TEST','TITLE     A TEST COMPLEX','AUTHOR    AUTHOR.ONE,AUTHOR.TWO','COMPND    MOL_ID: 1; MOLECULE: ALPHA; CHAIN: A; MOL_ID: 2; MOLECULE: BETA; CHAIN: B;','EXPDTA    X-RAY DIFFRACTION','REMARK   2 RESOLUTION.    2.10 ANGSTROMS.','HETNAM     ATP ADENOSINE TRIPHOSPHATE',atom('ATOM','ALA','A',1),atom('ATOM','GLY','B',1),atom('HETATM','ATP','A',4),atom('HETATM','HOH','A',5)].join('\n');
const p=api.parse(pdb,'pdb');assert.equal(p.chains.length,2);assert.equal(p.chains[0].description,'ALPHA');assert.equal(p.chains[1].description,'BETA');assert.equal(p.ligands[0].id,'ATP');assert.equal(p.waterCount,1);assert.equal(p.resolution,'2.10');assert.equal(p.authors,'AUTHOR.ONE,AUTHOR.TWO');assert.ok(!p.revision);
const bad=()=>api.parse('data_bad\n_struct.title\n;unclosed','cif');assert.throws(bad,/Unclosed/);
console.log('Passed metadata checks: real GFP mmCIF, multiline and quoted fields, complex chains, ligand/water counts, PDB records, missing values, malformed input.');
