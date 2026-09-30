# CladeScope

A browser app for comparing a protein reference sequence against phylogenetic clades and displaying per-site identity on a structure and sequence alignment.

The app can be found at [https://wrf.github.io/cladescope/](https://wrf.github.io/cladescope/) .

## Use the app

- The demo loads Aequorea GFP (1EMA, mmCIF), the supplied 13-sequence GFP alignment, and the supplied tree converted to rooted Newick with bootstrap annotations removed. The reference and structure-associated sequence is `Aequorea_victoria|AAA27721.1__WT-GFP`, chain A.
- Select a tree node to use its descendant sequences as the comparison group.
- Choose the reference sequence, color palette, and minimum non-gap fraction.
- Click alignment cells or mapped structure residues to inspect a site.
- Use **Show ligands**, **Show water**, and **Show other chains** independently. Water and ligands also respect the chain visibility setting.
- Ladderize by small or large clades first, or restore input ordering.
- Select a non-root node, then choose **Root at selected branch**. Rerooting places the root at the branch midpoint and resets comparison to all tree sequences. **Restore input root** restores the uploaded root.

### Your files

Open **Load your own dataset** and select:

1. An aligned protein FASTA with unique IDs and equal sequence lengths.
2. A Newick tree whose leaf IDs match FASTA IDs. Existing roots and `[&R]` comments are accepted. This is plain Newick, not a NEXUS file.
3. A PDB or mmCIF structure.

Choose the structure-associated alignment sequence and chain, then apply mapping. Verify the mapping agreement before interpreting colors. Sequence files are processed in the browser and are not uploaded by the app. Do not place personal datasets in `docs/` unless you intend to publish them with the website.

## Replace the demo

Place your structure, aligned FASTA and rooted Newick tree inside `docs/`, then edit `docs/demo.json` to specify their filenames, structure format (`pdb` or `cif`), reference sequence ID, structure-associated sequence ID and chain. Sequence IDs must match the alignment. The supplied tree keeps its original 3-versus-10 leaf root split. Bootstrap annotations are omitted.

## Files and development

- `docs/index.html`: interface and labels.
- `docs/app.js`: parsing, scoring, mapping, tree operations, and viewer interactions.
- `docs/style.css`: layout and visual styles.
- `docs/branding.js` and `docs/icon.svg`: branding.
- `docs/3Dmol-min.js`: bundled structure viewer; no CDN dependency.
- `docs/1EMA.cif`, `docs/demo-alignment.fasta`, `docs/demo-tree.nwk`: bundled GFP demo data.
- `docs/demo.json`: demo file paths, reference, structure-associated sequence, and chain.
- `tests/check.cjs`: logic checks, run with Node.js.

```bash
node --check docs/app.js
node tests/check.cjs
```

## Scoring and limitations

Reference identity and nine-bin palettes follow the identity/PyMOL-script behavior in [pdbcolor](https://github.com/wrf/pdbcolor). No structure B-factors are modified. Gaps and X are excluded from the denominator unless the non-gap fraction is below the cutoff, in which case the full group size is used. Columns with no valid residues or a reference gap have no score.

The tree is drawn as a cladogram: branch lengths are retained for rerooting but are not represented by horizontal distance. Structure mapping aligns observed chain residues to the structure-associated sequence; repeated regions and low sequence agreement need manual review. This prototype has limits of 25 MB per input file, 250,000 alignment cells, and 8 million mapping cells. It implements reference identity, not the repository's likelihood, heteropecilly, or other specialized analyses. Logic checks are included; automated browser interaction checks are not included.

## Third-party acknowledgements

3Dmol.js is bundled under its BSD license; see `docs/3Dmol-LICENSE.txt` and `docs/3Dmol-min.js.LICENSE.txt`. Scoring and palette definitions are adapted from pdbcolor. The demo structure is RCSB PDB entry 1EMA: https://www.rcsb.org/structure/1EMA. The demo alignment and rooted tree are supplied by the project owner.

This export contains no hosting credentials, private deployment configuration, or Git history. A project-wide license has not been assigned by this export; third-party notices remain applicable.
