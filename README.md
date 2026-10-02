# CladeScope

A browser app for comparing sequences between protein clades, or reference sequence against phylogenetic clades, and displaying per-site conservation on a structure and sequence alignment.

The app can be found at [https://wrf.github.io/cladescope/](https://wrf.github.io/cladescope/) .

![screenshot of cladescope comparison mode](https://github.com/wrf/cladescope/blob/main/docs/cladescope_comparison_screenshot.png)

## Comparison and Identity modes

**Comparison** is the default. Click a tree node, then choose **Assign to comparison clade 1** (blue) or **Assign to comparison clade 2** (red). Assignments are explicit and persist while clicking other nodes, ladderizing, or rerooting. The two root child groups are assigned initially when available.

Comparison scores use the **most frequent amino acid in each clade**, independently of the selected reference. Each clade can have a different dominant amino acid. Ties are retained in the site inspector.

- Blue = clade 1's most frequent amino-acid percentage, scaled to 255.
- Red = clade 2's most frequent amino-acid percentage, scaled to 255.
- Green = strongest shared amino-acid percentage: `max_AA min(frequency in clade 1, frequency in clade 2)`, scaled to 255. This is not pooled reference identity and is independent of clade sizes.
- Neither clade reaches **Minimum conservation (%)** = dark gray (`#555555`). The default is 50%, and the threshold is inclusive. At 100%, at least one clade must be fully conserved; different fixed residues in both are magenta and the same fixed residue in both is white.
- No usable observations in either clade = beige (no data).

**Minimum non-gap fraction** still controls the denominator, independently of the conservation display threshold. Frequencies exclude gaps and X unless non-gap coverage is below the cutoff, when the full group size is the denominator. Reference gaps do not mask Comparison results, although unmapped sites cannot be displayed on the structure.

| Clade 1 | Clade 2 | Color |
|---|---|---|
| 100% L | 100% L | White, `#ffffff` |
| 100% L | 100% V | Bright magenta, `#ff00ff` |
| 100% L | 50% L / 50% V | `#8080ff` |
| 50% L / 50% V | 50% L / 50% V | `#808080` at threshold 50; dark gray at thresholds above 50 |

The 100% L versus 100% V result remains magenta for unequal group sizes or a reference residue absent from both clades. The palette menu is disabled in Comparison mode because RGB channels are fixed.

Choose **Identity** next to Color Palette for the original single-clade workflow. Clicking a node immediately selects the identity group, and the four original palettes are enabled. Existing structure controls, mapping, loading, and site inspection are retained.

Alignment names have subtle blue/red backgrounds for clade membership. **Include overlapping species** defaults to off. If one assigned clade is a proper subclade of the other, the smaller clade stays intact and its species are removed from the larger clade for scoring and name highlights. For example, 3 species versus a parent containing 10 becomes 3 versus 7. Either assignment can be the parent. Turning the toggle on counts all selected species in both groups. Identical groups have no exclusive species; choose distinct clades or enable the toggle. Partial overlaps, possible after rerooting and retaining assignments, are removed from both groups when the toggle is off. The summaries show the effective count and original selection size. **Selected groups + reference only** filters the alignment to the union in Comparison mode, or the clicked clade in Identity mode, always keeping the reference visible.

The original app and original export remain separate; this is a new version.

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
