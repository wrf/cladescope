/* File-only PDB/mmCIF metadata extraction; no external lookup. */
(function(root){
  const clean=value=>value===undefined||value==='.'||value==='?'?'':String(value).trim();
  const waters=new Set(['HOH','WAT','H2O','DOD','SOL','TIP3','TIP3P']);
  const amino=new Set(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE']);
  function* tokenize(text){
    let i=0;
    while(i<text.length){
      if(/\s/.test(text[i])){i++;continue}
      if(text[i]==='#'){while(i<text.length&&text[i]!=='\n')i++;continue}
      if(text[i]===';'&&(i===0||text[i-1]==='\n')){
        i++;const start=i;const end=text.indexOf('\n;',i);if(end<0)throw Error('Unclosed mmCIF text field.');
        yield {value:text.slice(start,end).trim(),literal:true};i=end+2;continue;
      }
      if(text[i]==="'"||text[i]==='"'){
        const quote=text[i++],start=i;while(i<text.length&&!(text[i]===quote&&(i+1===text.length||/\s/.test(text[i+1]))))i++;
        if(i===text.length)throw Error('Unclosed mmCIF quoted value.');yield {value:text.slice(start,i),literal:true};i++;continue;
      }
      const start=i;while(i<text.length&&!/\s/.test(text[i]))i++;yield {value:text.slice(start,i),literal:false};
    }
  }
  function cifData(text){
    const iterator=tokenize(text),data={};let token=iterator.next().value;
    const advance=()=>{token=iterator.next().value};
    const keyword=t=>t&&!t.literal&&(/^(loop_|stop_|global_|data_|save_)/i.test(t.value)||t.value.startsWith('_'));
    while(token){
      if(!token.literal&&token.value.toLowerCase()==='loop_'){
        advance();const keys=[];while(token&&!token.literal&&token.value.startsWith('_')){keys.push(token.value.toLowerCase());advance()}
        if(!keys.length)throw Error('mmCIF loop has no column names.');keys.forEach(k=>data[k]=[]);let column=0;
        while(token&&!keyword(token)){data[keys[column]].push(token.value);column=(column+1)%keys.length;advance()}
        if(column)throw Error('Incomplete mmCIF loop row.');
      }else if(!token.literal&&token.value.startsWith('_')){
        const key=token.value.toLowerCase();advance();if(!token||keyword(token))throw Error('Missing mmCIF value for '+key);data[key]=[token.value];advance();
      }else advance();
    }
    return data;
  }
  function collectAtoms(rows,descriptions,componentNames){
    const chains=new Map(),ligands=new Map(),waterResidues=new Set(),models=new Set();
    for(const row of rows){
      const ch=clean(row.chain)||'(blank)',residue=ch+'|'+row.residue+'|'+row.name,model=clean(row.model)||'1';models.add(model);
      if(waters.has(row.name)){waterResidues.add(residue);continue}
      if(row.polymer){
        if(!chains.has(ch))chains.set(ch,{chain:ch,entities:new Set(),names:new Set(),residues:new Set()});
        const item=chains.get(ch);item.residues.add(residue);if(row.entity)item.entities.add(row.entity);const name=descriptions[row.entity]||row.description;if(name)item.names.add(name);
      }else{
        if(!ligands.has(row.name))ligands.set(row.name,{id:row.name,name:componentNames[row.name]||'',chains:new Set(),residues:new Set()});
        const item=ligands.get(row.name);item.chains.add(ch);item.residues.add(residue);
      }
    }
    return {chains:[...chains.values()].map(x=>({chain:x.chain,description:[...x.names].join('; '),observedResidues:x.residues.size,entities:[...x.entities].join(', ')})),ligands:[...ligands.values()].map(x=>({id:x.id,name:x.name,chains:[...x.chains].join(', '),copies:x.residues.size})),waterCount:waterResidues.size,modelCount:models.size};
  }
  function parseCif(text){
    const d=cifData(text),values=key=>(d[key]||[]).map(clean).filter(Boolean),first=key=>values(key)[0]||'';
    const types={},descriptions={},componentNames={};(d['_entity.id']||[]).forEach((id,i)=>{types[id]=clean(d['_entity.type']?.[i]);descriptions[id]=clean(d['_entity.pdbx_description']?.[i])});
    (d['_chem_comp.id']||[]).forEach((id,i)=>componentNames[id]=clean(d['_chem_comp.name']?.[i]));
    const atoms=d['_atom_site.group_pdb']||[],rows=[];
    for(let i=0;i<atoms.length;i++){
      const entity=clean(d['_atom_site.label_entity_id']?.[i]),name=clean(d['_atom_site.auth_comp_id']?.[i])||clean(d['_atom_site.label_comp_id']?.[i]);
      rows.push({chain:clean(d['_atom_site.auth_asym_id']?.[i])||clean(d['_atom_site.label_asym_id']?.[i]),entity,name,residue:(clean(d['_atom_site.auth_seq_id']?.[i])||clean(d['_atom_site.label_seq_id']?.[i]))+clean(d['_atom_site.pdbx_pdb_ins_code']?.[i]),model:d['_atom_site.pdbx_pdb_model_num']?.[i],polymer:types[entity]?types[entity]==='polymer':atoms[i]==='ATOM'||amino.has(name)});
    }
    const revisions=(d['_pdbx_audit_revision_history.revision_date']||[]).map((date,i)=>({date:clean(date),version:[clean(d['_pdbx_audit_revision_history.major_revision']?.[i]),clean(d['_pdbx_audit_revision_history.minor_revision']?.[i])].filter(Boolean).join('.')}));
    const latest=revisions.at(-1);
    const organisms=[...new Set([...values('_entity_src_gen.pdbx_gene_src_scientific_name'),...values('_entity_src_nat.pdbx_organism_scientific'),...values('_pdbx_entity_src_syn.organism_scientific')])];
    const citations=(d['_citation.id']||[]).map((id,i)=>({id,title:clean(d['_citation.title']?.[i]),doi:clean(d['_citation.pdbx_database_id_doi']?.[i])}));const citation=citations.find(x=>x.id==='primary')||citations[0];
    const software=[...new Set(values('_software.name'))];
    return {format:'mmCIF',entry:first('_entry.id'),title:first('_struct.title'),classification:first('_struct_keywords.pdbx_keywords'),authors:values('_audit_author.name').join('; ')||values('_citation_author.name').join('; '),deposited:first('_pdbx_database_status.recvd_initial_deposition_date'),released:revisions[0]?.date||'',revision:latest?latest.date+(latest.version?' (version '+latest.version+')':''):'',method:values('_exptl.method').join('; '),resolution:first('_refine.ls_d_res_high')||first('_reflns.d_resolution_high')||first('_em_3d_reconstruction.resolution'),organisms:organisms.join('; '),citation:citation?.title||'',doi:citation?.doi||'',software:software.join('; '),...collectAtoms(rows,descriptions,componentNames)};
  }
  function parsePdb(text){
    const lines=text.split(/\r?\n/),records=kind=>lines.filter(l=>l.slice(0,6).trim()===kind),join=kind=>records(kind).map(l=>l.slice(10).trim()).join(' ').replace(/\s+/g,' ');
    const header=records('HEADER')[0]||'',descriptionByChain={},componentNames={};
    let compound={};function flush(){if(compound.CHAIN)for(const ch of compound.CHAIN.split(',').map(x=>x.trim()))descriptionByChain[ch]=compound.MOLECULE||''}
    for(const field of join('COMPND').split(';')){const match=field.match(/^\s*([^:]+):\s*(.*)$/);if(!match)continue;const key=match[1].trim();if(key==='MOL_ID'){flush();compound={}}compound[key]=match[2].trim()}flush();
    for(const l of records('HETNAM')){const id=l.slice(11,14).trim();componentNames[id]=((componentNames[id]||'')+' '+l.slice(15).trim()).trim()}
    let model='1';const rows=[];for(const l of lines){const tag=l.slice(0,6).trim();if(tag==='MODEL')model=l.slice(10,14).trim();if(tag!=='ATOM'&&tag!=='HETATM')continue;const name=l.slice(17,20).trim(),chain=l.slice(21,22).trim();rows.push({name,chain,residue:l.slice(22,27).trim(),model,polymer:tag==='ATOM'||amino.has(name),description:descriptionByChain[chain]||''})}
    const revisions=records('REVDAT'),latest=revisions[0];const remark2=records('REMARK').find(l=>l.slice(7,10).trim()==='2'&&/RESOLUTION\.\s+[\d.]+\s+ANGSTROMS/.test(l));
    const sources=[...join('SOURCE').matchAll(/ORGANISM_SCIENTIFIC:\s*([^;]+)/g)].map(x=>x[1].trim());
    return {format:'PDB',entry:header.slice(62,66).trim(),title:join('TITLE'),classification:header.slice(10,50).trim(),authors:join('AUTHOR'),deposited:header.slice(50,59).trim(),released:revisions.at(-1)?.slice(13,22).trim()||'',revision:latest?latest.slice(13,22).trim()+' (revision '+latest.slice(7,10).trim()+')':'',method:join('EXPDTA'),resolution:remark2?.match(/RESOLUTION\.\s+([\d.]+)/)?.[1]||'',organisms:[...new Set(sources)].join('; '),citation:join('JRNL'),doi:records('JRNL').filter(l=>l.slice(12,16).trim()==='DOI').map(l=>l.slice(19).trim()).join(''),software:'',...collectAtoms(rows,{},componentNames)};
  }
  const api={parse:(text,format)=>format==='cif'?parseCif(text):parsePdb(text),cifData};root.StructureMetadata=api;if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof globalThis!=='undefined'?globalThis:window);
