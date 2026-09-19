# RDKit behaviours that change what you write

Measured against RDKit **2026.03.5**.
Where a behaviour is version-sensitive it says so; when in doubt, read the enum or the docstring off the installed build rather than recalling it:

```python
list(Chem.BondStereo.values.items())
list(Chem.ChiralType.values.items())
rdMolStandardize.Uncharger.__init__.__doc__
[m for m in dir(Chem.RWMol) if 'tereoGroup' in m]
```

Contents:
1. [SMILES syntax traps](#1-smiles-syntax-traps)
2. [Stereo perception](#2-stereo-perception)
3. [CIP labels](#3-cip-labels)
4. [Canonical atom ordering](#4-canonical-atom-ordering)
5. [CXSMILES and enhanced stereo](#5-cxsmiles-and-enhanced-stereo)
6. [Building stereo programmatically](#6-building-stereo-programmatically)
7. [Uncharger](#7-uncharger)
8. [Things that do not exist](#8-things-that-do-not-exist)

---

## 1. SMILES syntax traps

Valence errors print to **stderr** and `MolFromSmiles` returns `None`.
Any loop that doesn't test `mol is None` fails silently - you get a clean run and no molecules.
These four account for nearly all of them:

| Wrong | Why | Right |
|---|---|---|
| `C[NH3+]C` | A bracketed atom's H count is *fixed*, so this N has five bonds. | `C[NH2+]C` in-chain; keep `[NH3+]` terminal |
| `CC(=O)[O-]CCC` | The pendant anion silently continues the chain and bonds `[O-]` to the next carbon. | `CC(C(=O)[O-])CCC` |
| `CCS(=O)(=O)[O-]CC` | Same, for sulfonate. | parenthesise the whole group |
| `CC(C)(C)(C)C` | Four branches plus an incoming bond is five bonds. | count the incoming bond |

`[BH3-]` behaves like `[NH3+]`: it must terminate a chain or a branch.

**Bare `Cl` is hydrogen chloride, not a chlorine atom** - outside brackets it picks up an implicit H.
The atom is `[Cl]`, the anion `[Cl-]`.
Same for every halogen.

**Ring-closure digits must be unique among simultaneously open rings.**
Reusing an open digit parses cleanly and produces a different molecule with the same formula - see the tetraarylborate example in `motifs.md`.

**RDKit silently fills a charged atom's spare valence with an implicit H.**
Don't count substituents by reading the string; inspect the parsed mol.

**Small-ring double bonds silently lose their `/`…`\`.**
`C1CC/C=C\C1` parses without a warning and comes back as `C1=CCCCC1`.
Trans geometry is impossible below ring size 8, so RDKit drops the notation rather than erroring.

---

## 2. Stereo perception

`Chem.FindPotentialStereo(mol)` is the authority for "what stereo elements does this molecule have".
`bond.GetStereo()` and `Chem.FindMolChiralCenters(mol)` with default flags report only elements that already carry a marker, so they cannot answer "did I forget to mark something".

```python
for e in Chem.FindPotentialStereo(mol):
    kind = str(e.type)          # 'Atom_Tetrahedral', 'Bond_Double', ...
    spec = str(e.specified)     # 'Specified' | 'Unspecified' | 'Unknown'
    idx  = e.centeredOn
```

**`centeredOn` is an atom index for `Atom_*` and a bond index for `Bond_*`.**
Split on `'Atom' in kind` before you look anything up, or every index you print is wrong.

`StereoSpecified` has exactly three states - `Unspecified`, `Specified`, `Unknown`.
`Unknown` is what a wavy bond produces; it is a per-atom state, not a kind of stereo group.
Setting it strips the `@` from that atom by definition.

`FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)` is a convenient atoms-only view; unassigned centres come back with CIP `'?'`.

**Non-stereogenic double bonds are not reported at all** - `FindPotentialStereo` simply omits them.
Find them by subtraction:

```python
stereo_bonds = {e.centeredOn for e in Chem.FindPotentialStereo(mol) if 'Bond' in str(e.type)}
nonstereo = [b.GetIdx() for b in mol.GetBonds()
             if b.GetBondTypeAsDouble() == 2.0 and not b.GetIsAromatic()
             and b.GetIdx() not in stereo_bonds]
```

**Ring double bonds become stereogenic at ring size 8.**
Verified: `C1=CCCCC1` (6) and `C1=CCCCCC1` (7) report nothing; `C1=CCCCCCC1` (8) and `C1=CCCCCCCC1` (9) report an unspecified `Bond_Double`.

**C=N is stereogenic too.**
`CC(=N)N` (acetamidine, 4 heavy atoms) has one unspecified stereogenic bond.
Amidines and guanidines inside a bigger molecule are an easy place for an extra stereogenic bond to appear unnoticed.

---

## 3. CIP labels

Use `rdCIPLabeler.AssignCIPLabels(mol)` and read `_CIPCode` for any R/S or E/Z you state.
It implements the full CIP rules; the legacy `Chem.AssignStereochemistry` path does not.

```python
from rdkit.Chem import rdCIPLabeler
rdCIPLabeler.AssignCIPLabels(mol)
atom.GetPropsAsDict().get('_CIPCode')   # 'R' | 'S' | 'r' | 's'
bond.GetPropsAsDict().get('_CIPCode')   # 'E' | 'Z'
```

Measured divergence: the legacy implementation **reports nothing at all** for pseudo-asymmetric centres, where the modern one correctly emits lowercase `r`/`s`.

| Molecule | legacy | rdCIPLabeler |
|---|---|---|
| `O[C@H]1CC[C@@H](O)CC1` | `{}` | `{1: 's', 4: 's'}` |
| `C[C@]1(CC[C@@H](C)CC1)O` | `{}` | `{1: 'r', 4: 'r'}` |
| `C[C@H]1CC[C@H](C(C)C)CC1` | `{}` | `{1: 'r', 4: 'r'}` |

On ordinary alkene geometry `STEREOE`/`STEREOZ` agreed with `rdCIPLabeler` in every case tested here - 15,972 di-, tri- and tetrasubstituted alkenes, zero disagreements - so the legacy path is not simply wrong about E/Z.
There is still no reason to use it: it goes silent on pseudo-asymmetric centres, and the next point is worse.

**The cis/trans family is a different thing from E/Z, and reporting it as E/Z is wrong half the time.**
`STEREOCIS`/`STEREOTRANS`, and the `descriptor` field of `FindPotentialStereo` (`Bond_Cis`/`Bond_Trans`), are defined relative to `bond.GetStereoAtoms()` - whichever neighbours RDKit nominated - not to CIP priority.
Measured over 7,290 trisubstituted alkenes: `Bond_Trans` corresponded to CIP **E** in 3,645 cases and to CIP **Z** in the other 3,645.
It is a coin flip.
`CC/C=C(/C)C(=O)O` reports `Bond_Trans` and is CIP **Z**.

`descriptor` tells you the geometry relative to two named atoms, which is useful when you are setting geometry programmatically.
For anything you state to a reader as E or Z, read `_CIPCode` after `rdCIPLabeler.AssignCIPLabels`.

**Never compare raw `@`/`@@` or `/`/`\` characters between two strings.**
Parity tags are relative to neighbour write order, so a different traversal inverts the tag without changing the configuration.
Verified:

```
NC(=O)[C@H](F)C[NH3+]   CIP R at atom 3
NC[C@@H](F)C(N)=O       CIP R at atom 2      (same molecule, uncharged and re-canonicalised)
```

`@` became `@@` and the index moved.
Nothing was inverted.
Compare CIP labels.

---

## 4. Canonical atom ordering

**`MolFromSmiles` preserves the input token order.**
Atom index 0 is the first atom in the string, index 1 the second, and so on, through branches, ring closures and `.` separators.
Verified on branched and cyclic inputs.
So indices you compute from a parsed mol *are* indices into the string you parsed.

What does move is the **output** order when you write a canonical SMILES:

```python
Chem.MolToSmiles(mol)                                  # canonicalises
order = json.loads(mol.GetProp('_smilesAtomOutputOrder').replace(',]', ']'))
```

`_smilesAtomOutputOrder` is a *string* property - `GetPropsAsDict()` does not return it, so read it with `GetProp` and parse it.

Two rules, both measured:

- **Changing a formal charge genuinely reorders the canonical output**, moving atoms across symmetry classes.
`NC(=O)CCCN` writes as `NCCCC(N)=O`; protonate the terminal amine and `NC(=O)CCC[NH3+]` writes in input order.
Six of the seven atoms land in a position previously held by an atom of a different symmetry class.
- **Adding or changing stereo tags never does.**
Across every stereoisomer of six test molecules, no atom ever moved to a position held by an atom of a different symmetry class.
The output order can permute, but only among symmetry-equivalent atoms, which leaves the constitutional skeleton of the string unchanged.

The test that distinguishes these is not a naive diff of `_smilesAtomOutputOrder` but:

```python
ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))
crossings = [(i, j) for i, j in zip(order_a, order_b) if ranks[i] != ranks[j]]
```

A claim sometimes made - that any charged atom seizes the canonical root - **is false**.
Of seven singly-charged test molecules, only two began their canonical string at the charged atom; `CCCCCC[N+](C)(C)C` and `CCCCCCC(=O)[O-]` both start at a neutral carbon.
Charge changes the ordering; it does not determine the root.

**Practical rule:** align two forms of a molecule by mol atom index or by `_smilesAtomOutputOrder`, never by position in the canonical string.
Canonicalise *both* sides before comparing - diffing a hand-written input against a canonical output tells you nothing.

---

## 5. CXSMILES and enhanced stereo

Parse with CXSMILES enabled:

```python
ps = Chem.SmilesParserParams(); ps.allowCXSMILES = True
mol = Chem.MolFromSmiles(cxsmiles, ps)
```

**Verify the groups after parsing; a successful parse proves nothing.**
An index pointing at an atom that is not a stereocentre causes RDKit to **drop that group silently, with no warning and no error**.
Verified: a string with five `o` groups where one index was wrong parsed fine and returned four groups.

```python
assert len(mol.GetStereoGroups()) == expected
for g in mol.GetStereoGroups():
    print(g.GetGroupType(), [a.GetIdx() for a in g.GetAtoms()])
```

### Attaching groups programmatically

`Chem.CreateStereoGroup(...)` alone does not attach anything.
Go through an `RWMol` and call `SetStereoGroups` - there is no `Mol.SetStereoGroups`.

```python
rw = Chem.RWMol(mol)
g = Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_OR, rw, [3])
rw.SetStereoGroups([g])
mol = rw.GetMol()
```

`StereoGroupType` has exactly three members: `STEREO_ABSOLUTE`, `STEREO_AND`, `STEREO_OR`.
There is no "unknown" group type - unknown is a per-atom wavy-bond state (below).

### Index frames

CXSMILES indices are 0-based over every atom left to right, including bracket atoms, skipping bonds, parentheses and dots, and running **continuously across `.` fragment separators**.
Branches break any "+2 per stereocentre" pattern, so never hand-count them; enumerate `atom.GetIdx()` instead.

An index block is only meaningful against the exact string it is attached to.
**`MolToCXSmiles` recomputes the indices for whatever ordering it writes**, so you don't have to.
Verified round-trip:

```
in : C[C@H](O)[C@@H](N)[C@H](F)C(Cl)CC |o1:1,&1:3,w:5.5|
out: CCC(Cl)C(F)[C@H](N)[C@H](C)O      |w:4.4,o1:8,&1:6|
```

If the reader needs *their* atom ordering preserved, write with `canonical=False` and the indices are recomputed against that ordering too:

```python
wp = Chem.SmilesWriteParams(); wp.canonical = False
Chem.MolToCXSmiles(mol, wp)
```

Sanity-check by confirming the output differs from the default canonical string and still round-trips to the same canonical form.

### Group numbering, and what counts

`o1`, `o2`, … and `&1`, `&2`, … are per-type sequential counters, and RDKit renumbers them on output: `|o5:1,o9:4,&3:7|` is accepted on input and comes back as `|o1:…,o2:…,&1:…|`.
Don't read meaning into the numbers you get back, and don't try to preserve a particular numbering across a round-trip.

Which tokens actually work:

| Token | Behaviour |
|---|---|
| `o1:`, `&1:` | OR / AND groups. Parsed, stored, re-emitted. |
| `a:` | Absolute group. **Valid** - parsed as `STEREO_ABSOLUTE` and re-emitted. |
| `w:atom.bond` | Wavy bond → `StereoSpecified.Unknown`. Re-emitted for atoms. |
| `u:`, `?:` | **Silently ignored.** Parse succeeds, group count is unchanged. |

Stereocentres in no group are implicitly ABS, so `a:` is usually unnecessary rather than invalid - but it is available if a reader's parser wants it explicit.

**A single-member OR group is not a no-op.**
It means "this centre is R or S and we don't know which", and RDKit treats it that way: `C[C@H](O)CC |o1:1|` enumerates to two stereoisomers and round-trips intact.
Don't drop one on the grounds that it does nothing.

### Wavy bonds (unknown stereo)

`w:` is how unknown stereo is expressed - there is no `?` group type.
The syntax is `w:atom.bond`: the atom the wavy bond starts from, then the bond it runs along.

```
CCC(O)C |w:2.1|   → atom 2 is Atom_Tetrahedral, Unknown
```

For **atoms** this is fully supported: it parses, sets `StereoSpecified.Unknown`, coexists with `o`/`&` groups in the same block, is re-emitted by the writer with recomputed indices, and counts as unassigned for stereoisomer enumeration.
Applying it to an atom that carried a `@` strips the tag, which is the definition rather than a bug.

For **bonds** it is lossy.
`CC=CCC |w:2.1|` parses and reports the double bond as `Unknown`, but the writer emits plain `CC=CCC` with no `|w:|` block, and re-parsing that gives `Unspecified` instead.
If a spec needs an explicitly unknown *double bond* to survive a round-trip, say that it can't be expressed here rather than shipping a string that quietly downgrades.

### The one real gap

**The writer serialises only the *atom* members of a stereo group.**
`CreateStereoGroup` accepts a `bondIds` argument and the mol object retains those bonds, but neither the default writer nor `CXSmilesFields.CX_ALL` emits them:

```python
g = Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_AND, rw, [3], [1])
# mol keeps bond 1 in the group; MolToCXSmiles writes  C/C=C/[C@H](O)CC |&1:3|
```

So a double bond cannot be placed in an enhanced stereo group in a way that survives serialisation.
If a spec asks for that, say plainly that it can't be expressed rather than shipping something lossy.

---

## 6. Building stereo programmatically

Setting bond stereo has a strict ordering requirement, and getting it wrong is worse than an exception:

```python
bond.SetStereoAtoms(a, b)                    # each end atom's OTHER neighbour
bond.SetStereo(Chem.BondStereo.STEREOTRANS)  # or STEREOCIS / STEREOE / STEREOZ
Chem.SetDoubleBondNeighborDirections(mol)    # harmless once stereo atoms are set
```

**Always call `SetStereoAtoms` first.**
Two different things happen if you don't:

- `SetStereo(STEREOCIS)` or `STEREOTRANS` raises a pre-condition violation ("Stereo atoms should be specified before specifying CIS/TRANS") - noisy, but it stops you.
- `SetStereo(STEREOE)` or `STEREOZ` does *not* raise.
It writes bare `CC=CC`, and a subsequent `Chem.SetDoubleBondNeighborDirections(mol)` **segfaults the interpreter** - no exception, no traceback, the process dies.
If a build script vanishes without output, this is why.

With stereo atoms set, `SetStereo` writes the directional bonds immediately: `SetStereoAtoms(0, 3)` then `SetStereo(STEREOE)` on `CC=CC` gives `C/C=C/C` without needing anything further.
`SetDoubleBondNeighborDirections` is then a no-op, and calling it is a cheap way to cover mols built atom-by-atom with `RWMol.AddAtom`/`AddBond`, which do need it.

**Clear double-bond directions locally, never globally.**
Setting every bond's direction to `NONE` to unspecify one double bond destroys E/Z everywhere else - verified: a molecule with two specified E bonds ended with both unspecified.
Clear only the directional single bonds adjacent to the target bond, re-run perception, and audit the rest.

---

## 7. Uncharger

```python
from rdkit.Chem.MolStandardize import rdMolStandardize
rdMolStandardize.Uncharger(canonicalOrder=True, force=False, protonationOnly=False)
```

Those are the real parameter names and the real positional order - note it is `canonicalOrder`, not `canonicalOrdering`, and that the first positional argument is *not* `protonationOnly`.
Always pass them by keyword.

- **`protonationOnly=True`** restricts the operation to proton transfer.
This is what makes it a test for "removable charge".
Without it the Uncharger will do things proton transfer cannot - it reduces `CC[BH3-]` to neutral `BCC`.
- **`force=True`** makes it a per-atom test.
Without it the Uncharger preserves overall charge balance: a betaine `C[N+](C)(C)CCC(=O)[O-]` (net 0) comes back untouched, so its perfectly removable carboxylate looks permanent.
With `force=True` it protonates to net +1.
- **Parse with `removeHs=True`** (the default).
The Uncharger edits hydrogen *counts* and cannot delete a hydrogen *node*, so with `removeHs=False` it is a silent no-op on explicit H: `CC([N+]([H])([H])[H])C` comes back at +1.

Selected measured results with `protonationOnly=True, force=True`:

| Input | Output | Verdict |
|---|---|---|
| `CS(=O)(=O)[O-]` | `CS(=O)(=O)O` | removable |
| `[O-]S(=O)(=O)c1ccccc1` | `O=S(=O)(O)c1ccccc1` | removable |
| `CC(=O)[O-]` | `CC(=O)O` | removable |
| `CC[NH3+]` | `CCN` | removable |
| `[NH4+]` | `N` | removable |
| `[Cl-]` | `Cl` | removable - **and now hydrogen chloride** |
| `C[N+](C)(C)C` | unchanged | permanent |
| `C[P+](C)(C)C` | unchanged | permanent |
| `C[B-](C)(C)C` | unchanged | permanent |
| `CC[B-](F)(F)F` | unchanged | permanent |
| `[Na+]`, `[Ca+2]` | unchanged | permanent |
| `F[P-](F)(F)(F)(F)F` | neutral P | **removable** |

Bare halide counter-ions turning into hydrogen halides is a trap for any "neutralise then compare" workflow - the fragment count is unchanged but the formula isn't.

Note `F[P-](F)(F)(F)(F)F` in that table: the Uncharger will happily produce chemically ridiculous neutral forms (here a `[PH]` with six fluorines).
The test answers "can proton transfer remove this charge in RDKit's model", which is a good proxy for removable-vs-permanent and not a chemical judgment.
Say which one you are asserting.

**Check the uncharged form's stereo too, not just the submitted one.**
A `@` tag written on an atom that is not chiral in the state you submitted is silently stripped (correct behaviour, easily misread as a bug), and conversely the neutralised form can gain or lose an element.
It costs one extra call to run `FindPotentialStereo` on both.

**The default (non-`force`) Uncharger can invent a stereocentre.**
Balancing charge, it protonates only as many groups as it needs to, and on two symmetric removable arms it protonates exactly one - which breaks the symmetry that was keeping the shared carbon achiral:

```
C[N+](C)(C)CCC(CC(=O)[O-])CC(=O)[O-]      as written:  0 stereo elements
  Uncharger(protonationOnly=True)      →  net 0, one arm protonated, 1 unspecified centre
  Uncharger(protonationOnly=True, force=True) → net +1, both arms protonated, 0 centres
```

If a reader's pipeline runs the default configuration, that phantom centre is what they will see.
Walk the protonation states that matter, not only the one you submitted.

**Re-parse before you deliver.**
After any in-place edit - `uncharge`, setting tags, clearing directions - the mol can carry a chiral tag on an atom that is no longer a stereocentre, and the writer emits it:

```
C[C@@H](C(=O)O)C(=O)[O-]  --uncharge-->  writes "C[C@@H](C(=O)O)C(=O)O"
                                          but FindPotentialStereo says 0 elements
re-parsed                                 canonicalises to "CC(C(=O)O)C(=O)O"
```

The string you paste and the mol you measured are only the same object if you round-tripped it.

---

## 8. Things that do not exist

- **`rdkit.Chem.Atropisomers`** - `ImportError`.
`BondStereo` does have `STEREOATROPCW`/`STEREOATROPCCW`, but atropisomerism cannot be perceived from bare SMILES at all: it needs 2D coordinates with wedge/hash bonds (a molfile), and only then round-trips through CXSMILES.
`EnumerateStereoisomers` does not handle it.
- **`Mol.SetStereoGroups`** - it's on `RWMol` only.
- **A `?` stereo-group type** - unknown stereo is the per-atom `w:` wavy-bond state, not a group.

Limits worth knowing before you promise something:

- **Allene and cumulene axial chirality cannot be specified in SMILES here.**
It is a bond property, not an `@AL` atom tag - but RDKit 2026.03.5 discards the directional markers entirely on parse.
Every one of `F/C(F)=C=C(\Cl)Br`, `F/C=C=C/F`, `C/C=C=C/C` and `F/C=C=C=C/F` comes back with the slashes stripped and each cumulated bond reported `Unspecified`; no atom carries a chiral tag.
So an allene is a reliable way to *add unspecified stereogenic bonds*, and not a way to add a specified one.
If a spec wants specified axial chirality, say it needs a molfile with wedge bonds.
- **Non-tetrahedral stereo essentially requires metals.**
`@SP`/`@TB`/`@OH` are coordination-compound notation.
A sulfoxide `C[S@@](=O)CC` is ordinary tetrahedral stereo with the lone pair as the fourth substituent - not a non-tetrahedral counterexample.
RDKit does perceive `Atom_SquarePlanar` on `Cl[Pt@SP3](Cl)([NH3])[NH3]`, but normalises the molecule to dative form, `[NH3]->[Pt@SP3](<-[NH3])([Cl])[Cl]`, and `->`/`<-` is an OpenSMILES extension many parsers reject.
- **`EmbedMolecule` is a screening heuristic, not a proof of chemical possibility.**
Its verdict depends on configuration: `useBasicKnowledge=True` adds planarity and ring physics and rejects strained-but-real molecules; without a `randomSeed` it is not even deterministic, and `maxAttempts` is a ceiling paid only on rejection, with each attempt under a fixed seed being a distinct deterministic draw.
A single-attempt failure is never a valid negative verdict.
If you use it, state the exact settings and call the result a screen.

When a reader names a downstream consumer - a depicter, their own parser - an RDKit-clean string is necessary but not sufficient.
Say which tool you actually tested.
