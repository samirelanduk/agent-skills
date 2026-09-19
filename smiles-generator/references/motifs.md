# Construction motifs

Building blocks whose behaviour follows from symmetry.
Every SMILES here has been parsed and its stated property measured in RDKit.

**This is a starting set, not a taxonomy of what can be asked for.**
It covers the axes that come often - charges, stereo, multi-fragment, tautomer, meso - each with a trap that costs a build cycle to rediscover.
Plenty of specs ask for none of them - a ring system, an element count, a molecular-weight window, a particular functional group, an isotope label.
Those need no catalogue: write the group, measure it, attach it.

To get a property P reliably:

1. Write the smallest fragment you believe has P.
2. Measure P on that fragment alone, in RDKit, before attaching it to anything.
3. Prefer a construction where P follows from symmetry or connectivity rather than from a marker you wrote - markers get stripped, symmetry doesn't.
4. Re-measure after assembly, because attachment changes the neighbourhood that P depended on.

Section 9 covers everything else.

Contents:
1. [Charges](#1-charges)
2. [Stereo elements on demand](#2-stereo-elements-on-demand)
3. [Killing stereocentres you don't want](#3-killing-stereocentres-you-dont-want)
4. [Charge-dependent stereogenicity (twin-arm gadgets)](#4-charge-dependent-stereogenicity)
5. [Meso forms](#5-meso-forms)
6. [Tautomer pairs](#6-tautomer-pairs)
7. [Multi-fragment specs](#7-multi-fragment-specs)
8. [Hitting a target stereoisomer count](#8-hitting-a-target-stereoisomer-count)
9. [Everything else: composition, rings, substructure, size](#9-everything-else)

---

## 1. Charges

The classification below is what `Uncharger(protonationOnly=True, force=True)` does to each group - the operational test for removable vs permanent.

### Permanent (survive the Uncharger)

| Motif | Example | Note |
|---|---|---|
| Quaternary ammonium | `CCCC[N+](C)(C)C` | The workhorse. Three identical methyls add no stereo. |
| Quaternary phosphonium | `CCC[P+](C)(C)C` | |
| Sulfonium | `C[S+](C)CCC` | |
| N-alkyl pyridinium | `CCCC[n+]1ccccc1` | Compact, aromatic, adds no stereo. |
| Tetraalkylborate | `CCC[B-](C)(C)C` | The standard permanent negative. |
| Organotrifluoroborate | `CCC[B-](F)(F)F` | Very compact. |
| Tetrafluoroborate | `F[B-](F)(F)F` | Counter-ion form. |
| Tetraarylborate | `c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1` | Large - poor for depiction. |
| Alkali/alkaline-earth cation | `[Na+]`, `[K+]`, `[Ca+2]` | Survives the Uncharger, but it's an ion, not a functionalised group. |

Note the ring-closure digits in the tetraarylborate.
Writing `c1` for the outer ring and reusing it while it is still open - `c1ccc([B-](c1ccccc1)(c1ccccc1)c1ccccc1)cc1` - **parses without error and gives a completely different molecule**, a fused polycyclic cage, because a digit that is still open gets matched to the next occurrence rather than starting a new ring.
Same atom count, same formula, so neither an atom tally nor a parse check catches it.
Use a fresh digit for each simultaneously-open ring, and read the canonical output back.

**Compact template for two permanent negatives on one fragment:**
`C[B-](C)(C)CC[B-](C)(C)C` - net −2, zero stereo elements.

### Removable (the Uncharger neutralises them)

| Motif | Example |
|---|---|
| Protonated primary amine | `CCCC[NH3+]` |
| Protonated secondary amine | `CCC[NH2+]CC` |
| Ammonium | `[NH4+]` |
| Carboxylate | `CCCC(=O)[O-]` |
| Phenoxide / alkoxide | `c1ccccc1[O-]` |
| **Sulfonate** | `CCCS(=O)(=O)[O-]` |

Sulfonate belongs in this list.
Its pKa is extremely low, so it is tempting to call it permanent, but the Uncharger protonates it (`CS(=O)(=O)[O-]` → `CS(=O)(=O)O`), which means any pipeline using that test will classify it as removable.
If a spec asks for a permanent negative, use borate.

### Two traps

- `[BH3-]` survives `protonationOnly=True` but **plain `Uncharger()` reduces it to neutral boron** (`CC[BH3-]` → `CCB`).
If you don't control which Uncharger configuration the reader runs, prefer tetraalkyl- or trifluoroborate.
- Hexafluorophosphate is **not** permanent under the Uncharger - `F[P-](F)(F)(F)(F)F` comes back neutral, and RDKit also perceives an unspecified stereo element on the phosphorus.
Use tetrafluoroborate `[B-](F)(F)(F)F` instead when you want an inert permanent anion.

### Mixed on one fragment

Betaine style: `[O-]C(=O)C[N+](C)(C)CC(=O)[O-]` - one permanent positive, two removable negatives, fragment charge −1.

---

## 2. Stereo elements on demand

| Want | Motif | Verified |
|---|---|---|
| Specified stereocentre | `CC[C@H](F)CCO` | 1 specified atom |
| Unspecified stereocentre | `CCC(F)CCO` | 1 unspecified atom |
| Specified stereogenic C=C | `CC/C=C/CC` | 1 specified bond, CIP E |
| Unspecified stereogenic C=C | `CCC=CCC` | 1 unspecified bond |
| Non-stereogenic C=C | `C(=C)CCO` | terminal vinyl: two H on one end |
| Non-stereogenic C=O | `CCC(=O)CC` | any ketone, ester or carboxyl |
| Unknown (wavy) stereo | `CCC(O)C |w:2.1|` | `StereoSpecified.Unknown` |

**Stereogenic ring double bonds need a ring of 8 or more.**
`C1=CCCCCC1` (7) is not stereogenic; `C1=CCCCCCC1` (8) is.
RDKit silently strips `/`…`\` written on a smaller ring double bond - `C1CC/C=C\C1` parses without complaint and comes back as plain `C1=CCCCC1`.
Keep stereogenic C=C acyclic unless the spec demands otherwise.

**C=N counts.**
Acetamidine `CC(=N)N` has one unspecified stereogenic bond.
A guanidinium or amidine buried in a larger molecule is an easy place for an extra stereogenic bond to appear unnoticed.
Note that a `=N` carrying only an implicit H has no explicit single bond for `/` or `\` to attach to, so it cannot be *specified* from bare SMILES - substitute the nitrogen (`CC(=NC)N`) if you need to pin it.

**Packing many elements:** alternate a directional double bond with a tagged carbon along a chain, and cap with an ester or terminal vinyl to supply the non-stereogenic double bonds:

```
CC/C=C\[C@@H](F)/C=C\[C@@H](F)C(=O)OC
→ 2 specified centres, 2 specified bonds, 1 non-stereogenic C=O, 14 heavy atoms
```

**Conjugated double bonds work, but check them.**
`C/C=C/C=C/C` gives two specified E bonds.
The hazard is that the single bond between them carries one direction serving *both* double bonds, so a slash you change to fix one bond flips the other.
RDKit will also rewrite your slashes: `C/C=C/C=C\C` comes back as `C/C=C\C=C\C` - same molecule, different notation.
**Read CIP labels from `rdCIPLabeler`, never the slash characters**, to confirm what you built.

---

## 3. Killing stereocentres you don't want

The reliable way is symmetry, not saturation.
A carbon with four attachments is fine as long as two of them are identical.

- **Geminal identical substituents**: `CCC(C)(C)CO` - four attachments, no stereocentre.
This is how you park an odd number of charged groups without creating stereo: `CC([NH3+])([NH3+])[NH3+]` has three cations and zero stereo elements.
- **Aromatic hub + CH2 linkers**: aromatic ring carbons are sp2, and a benzylic CH2 carries two hydrogens.
Neither can be a stereocentre no matter what you hang off it.
- **Identical methyls on the charged atom**: `[N+](C)(C)C`, `[P+](C)(C)C`, `[B-](C)(C)C`, `[B-](F)(F)F` - charge with no stereo.

Verified zero-stereocentre scaffold carrying mixed charges:

```
c1cc(CCO)cc(CC(=O)[O-])c1CC([NH3+])([NH3+])[NH3+]
→ 0 stereo atoms, 0 stereo bonds, net +2, 1 non-stereogenic C=O
```

**A quaternary ammonium with four *different* groups is itself a stereocentre.**
`CC[N+](C)(CCC)CCCC` has one unspecified centre; `CC[N+](C)(CC)CCC` (two ethyls) has none.
Useful either way - as a stereocentre that cannot be neutralised away, or as an accident to avoid.

**Move charged pendant groups onto a CH2 linker** rather than attaching them directly to a backbone CH, when direct attachment would create an unwanted centre.

---

## 4. Charge-dependent stereogenicity

Make an element stereogenic *only* in one charge state by giving an atom two arms that are constitutionally identical except for the charge.

**Removable-charge-dependent stereocentre** - one arm `-COOH`, one arm `-COO⁻`:

```
OC(=O)C[C@@H](C)CC(=O)[O-]   → 1 specified stereocentre, net −1
OC(=O)CC(C)CC(=O)O           → 0 stereocentres   (arms now identical)
```

Minimal version: `C[C@@H](C(=O)O)C(=O)[O-]` (methylmalonate mono-anion).

**Permanent-charge-dependent stereocentre** - twin arms differing by a charge-enabled methyl, so removing the charge would force the extra substituent off and the arms would become identical:

```
CN(C)CC[C@@H](C)CC[N+](C)(C)C   → 1 specified stereocentre, net +1, survives Uncharger
CN(C)CCC(C)CCN(C)C              → 0 stereocentres   (the hypothetical neutral form)
```

**Removable-charge-dependent double bond**: put `[NH3+]` and `N` on the same sp2 carbon; protonation state is the only thing distinguishing the two ends.

**Inert charge ballast** that adds no stereo, for pushing net charge off zero when your working gadgets happen to cancel: `C[N+](C)(C)C`.

Whichever you use, **run the counterfactual** - build the other charge state and confirm exactly one element appears or disappears and nothing else changes.

For a removable-charge gadget the Uncharger gives you the counterfactual for free.
For a permanent-charge gadget it cannot: the whole point is that proton transfer doesn't touch it, so `CN(C)CC[C@@H](C)CC[N+](C)(C)C` comes back from the Uncharger unchanged, still +1 and still one centre.
Write the hypothetical neutral analogue by hand and count that instead.

One caution when two removable-charge arms sit on the same atom: if they are identical, the centre exists only in the charged form, which is the point - but make sure that's what was asked.
If the spec wants a centre that is real in *every* protonation state, make the arms structurally different (acetate arm vs propionate arm):

```
C[C@@H](CC(=O)[O-])CCC(=O)[O-]   → 1 specified stereocentre, net −2
C[C@H](CCC(=O)O)CC(=O)O          → 1 specified stereocentre, net  0   (after Uncharger)
```

---

## 5. Meso forms

**Writing the same chiral tag on both centres of a symmetric chain gives the meso form; opposite tags give the chiral one.**
Reasoning from R/S to tag characters gets this backwards.
Verified on tartaric acid:

| SMILES | CIP | Form |
|---|---|---|
| `OC(=O)[C@H](O)[C@H](O)C(=O)O` | 3R, 5S | **meso** |
| `OC(=O)[C@H](O)[C@@H](O)C(=O)O` | 3R, 5R | (R,R), chiral |

Same pattern in 2,3-butanediol: `C[C@H](O)[C@H](O)C` is meso (S,R); `C[C@H](O)[C@@H](O)C` is (S,S).

**Testing for meso.**
Invert every chiral tag and every bond stereo, then compare.
InChI equality and canonical-SMILES equality both gave correct answers on every case tested here, but InChI is the safer primary test; use canonical SMILES as a cheap pre-check.

```python
def invert(m):
    m2 = Chem.Mol(m)
    for a in m2.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return m2
# meso  <=>  MolToInchi(m) == MolToInchi(invert(m))
```

Cross-check with CIP: a meso pair shows opposite descriptors (R and S) on the two centres.
Note that `rdCIPLabeler` also emits **lowercase `r`/`s`** for pseudo-asymmetric centres - e.g. cis-cyclohexane-1,4-diol `O[C@H]1CC[C@@H](O)CC1` is (1s,4s).
The legacy `AssignStereochemistry` path reports *nothing at all* for those atoms.

**Symmetry test set** (assignments → unique structures, via `EnumerateStereoisomers(onlyUnassigned=True, unique=False, maxIsomers=0)`):

| Molecule | SMILES | assignments | unique |
|---|---|---|---|
| tartaric acid | `OC(=O)C(O)C(O)C(=O)O` | 4 | 3 |
| 2,3-butanediol | `CC(O)C(O)C` | 4 | 3 |
| cyclohexane-1,4-diol | `OC1CCC(O)CC1` | 4 | 2 |
| 2,3,4-pentanetriol | `CC(O)C(O)C(O)C` | 8 | 4 |
| inositol | `OC1C(O)C(O)C(O)C(O)C1O` | 64 | 9 |

Inositol breaks code that assumes 2-to-1 collapses: its multiplicities are not all 2.

---

## 6. Tautomer pairs

The toolkit's rules define the answer.
Test the pair with the hash the reader will use:

```python
from rdkit.Chem import RegistrationHash
from rdkit.Chem.RegistrationHash import HashLayer
RegistrationHash.GetMolLayers(mol)[HashLayer.TAUTOMER_HASH]
```

`GetMolLayers` returns a plain **dict** keyed by the `HashLayer` enum - `layers.tautomer_hash` raises `AttributeError`.
Available layers: `CANONICAL_SMILES`, `ESCAPE`, `FORMULA`, `NO_STEREO_SMILES`, `SGROUP_DATA`, `TAUTOMER_HASH`, `NO_STEREO_TAUTOMER_HASH`.
The tautomer hash retains stereo, so enantiomers hash differently; `NO_STEREO_TAUTOMER_HASH` is the stereo-free variant.

Cheap pre-check first: if the two canonical SMILES are already identical you haven't got a pair at all, you've got one molecule written twice.

**Verified as sharing a tautomer hash:**

| Pair | Type |
|---|---|
| `C[C@@H](F)C(N)=O` / `C[C@@H](F)C(=N)O` | amide/iminol - also carries a stereocentre |
| `O=c1cc[nH]cc1` / `Oc1ccncc1` | 4-pyridone |
| `O=c1cccc[nH]1` / `Oc1ccccn1` | 2-pyridone |
| `O=c1[nH]cnc2[nH]cnc12` / `Oc1ncnc2[nH]cnc12` | hypoxanthine |
| `NC1=NC=CC=N1` / `N=C1N=CC=CN1` | aminopyrimidine |
| `O=C1NC=CC(=O)N1` / `OC1=NC=CC(=O)N1` | uracil |

**Verified as NOT sharing a tautomer hash** - do not offer these:

- **Keto/enol of any kind.**
`CC(=O)C` / `CC(O)=C` hash differently.
RDKit's heteroatom tautomer rules only collapse certain heteroatom tautomers, and carbon-to-oxygen proton shifts are not among them.
- **2-pyridone written in Kekulé form.**
`O=C1C=CC=CN1` / `OC1=CC=CN=C1` fail, while the aromatic lowercase spelling of the same pair succeeds.
For ring tautomers, always write the aromatic form.

`TautomerEnumerator()` destroys stereo by default - `CleanupParameters` has `removeSp3Stereo` and `removeBondStereo` both set to `True`.

---

## 7. Multi-fragment specs

- Keep the interesting chemistry on one fragment and use simple counter-ions for the rest: a stereogenic dianion plus two `[Na+]` reads naturally.
- Write metal salts fully disconnected - `[Na+].CC(=O)[O-]`, not Na–O bonded.
- Keep every fragment charged when the set should read as one salt; a neutral fragment sitting alongside charged ones invites a "mixture" reading.
(Solvates are the deliberate exception.)
- Reach the target net charge by adjusting counter-ion *multiplicity* rather than restructuring a fragment.
- When fragments must be "different", differ in heteroatom set, chain length and branching topology.
Two fragments with the same formula that differ by a CH2 and a positional swap will be called basically the same thing.
- For a repeated-fragment spec, repeat something substantive - a ligand or a real counter-ion (`[Pt+2].[Cl-].[Cl-].NCCN.NCCN`) rather than three waters.
- Atom indices run **continuously across `.` separators** and do not reset per fragment, so a stereo group can legitimately span fragments: `[NH3+][C@@H](F)[C@H](Cl)CC.[O-][C@@H](Br)[C@H](I)CC |o1:1,8|` parses to one OR group holding atoms 1 and 8.
- Isolate the largest fragment before running CIP or stereo analysis if you only care about the organic part - otherwise counter-ions pollute the counts: `max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())`.

See `salts-solvates.md` for validated counter-ion and solvate structures.

---

## 8. Hitting a target stereoisomer count

For an **asymmetric** molecule the count is exactly

```
2 ^ (unspecified stereo elements + enhanced stereo groups)
```

Each AND/OR group contributes one factor of 2 regardless of how many atoms it holds, and an atom marked `Unknown` by a wavy bond counts as unassigned.
Verified:

| SMILES | unspec | groups | count |
|---|---|---|---|
| `CC(F)C(Cl)C` | 2 | 0 | 4 |
| `CC=CC(F)C` | 2 | 0 | 4 |
| `C[C@H](O)[C@@H](N)CC \|o1:1,3\|` | 0 | 1 | 2 |
| `C[C@H](O)[C@@H](N)CC \|o1:1,&1:3\|` | 0 | 2 | 4 |
| `C[C@H](O)[C@@H](N)C(F)CC \|o1:1,&1:3\|` | 1 | 2 | 8 |

So to hit 2^14, supply fourteen independent unspecified centres, or mix in groups.
A chain of `C(X)` carbons each with a distinct substituent and no `@` gives you unspecified centres cheaply: `NC(F)C(F)C(F)...O`.

**The formula only holds when the molecule has no internal symmetry.**
With `unique=True` (the default) a symmetric skeleton collapses: a symmetric ten-centre chain gives 528 unique isomers, not 1024.
If the spec is a count, either keep the termini different, or measure rather than predict.

**Always pass `maxIsomers=0`.**
The default silently truncates at exactly 1024, and 1024 is itself a plausible true answer, so a capped run is indistinguishable from a real one.

```python
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=0)
n = len(list(EnumerateStereoisomers(mol, opts)))
```

To emit one fully-specified isomer from an unspecified skeleton, use `StereoEnumerationOptions(maxIsomers=1, onlyUnassigned=True)`, take the first, re-parse it and re-run `FindPotentialStereo` to confirm nothing is left unassigned.

---

## 9. Everything else

These axes need no motif catalogue - write the thing and measure it.
The notes below are where measurement is less obvious than construction.

### Substructure and functional-group requirements

Express the requirement as SMARTS and count matches; don't eyeball the string.

```python
q = Chem.MolFromSmarts('[SX4](=O)(=O)[NX3]')   # sulfonamide
len(mol.GetSubstructMatches(q, uniquify=True))
```

`Chem.MolFromSmarts` returns `None` on a bad pattern exactly as `MolFromSmiles` does, so a typo'd query silently matches nothing and the requirement looks satisfied at zero.
Check the query parsed before trusting a count of 0.

Two counting points:

- **`uniquify=True` deduplicates by atom set, not by chemistry.**
A symmetric query on a symmetric molecule still collapses matches you might have meant to count separately.
- **Aromatic vs aliphatic is a real distinction in SMARTS.**
`C` matches only aliphatic carbon and `c` only aromatic; `[#6]` matches either.
A pattern written with `C` that "should" match a benzene carbon is the usual cause of an unexpected zero.

### Ring systems

`mol.GetRingInfo()` gives `NumRings()`, `AtomRings()` and `BondRings()`; the count is the smallest set of smallest rings, which is what "how many rings" normally means but is not the same as the number of ring-closure digits you typed.
For fused systems, `rdMolDescriptors.CalcNumAromaticRings`, `CalcNumAliphaticRings` and `CalcNumRings` split it up.

Building a ring correctly is mostly about ring-closure digits: **a digit is only reusable once its ring has closed.**
Reusing one that is still open parses cleanly and gives a different molecule with the same formula, so neither the parse nor an atom count catches it:

```
c1ccc(-c2ccccc2)cc1   → biphenyl, C12H10, two separate rings
c1ccc(-c1ccccc1)cc1   → a fused tricycle, also C12H10 - not biphenyl
```

Read the canonical output back whenever a structure has more than one ring.

### Size and property targets

`Descriptors` and `rdMolDescriptors` cover the usual asks - `MolWt`, `NumRotatableBonds`, `NumHDonors`, `NumHAcceptors`, `MolLogP`, `TPSA`, `CalcNumAromaticRings`.
Two things to keep straight:

- **`MolWt` includes implicit hydrogens; `GetNumHeavyAtoms` does not count them at all.**
A heavy-atom budget and a MW budget are different constraints.
- **Multi-fragment molecules are measured as a whole.**
For "the organic part must be under 400", isolate it first: `max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())`.

### Isotopes, atom maps and dummy atoms

| Want | Write | Read back |
|---|---|---|
| Isotope label | `[13CH3]CO` | `atom.GetIsotope()` |
| Atom map number | `[CH3:1]CO` | `atom.GetAtomMapNum()` |
| Attachment point / R-group | `*CCO` or `[*:1]CCO` | `atom.GetAtomicNum() == 0` |

All three survive a canonical round-trip.
Note that atom maps and dummy atoms change canonical output ordering, and a dummy atom counts as a heavy atom in `GetNumHeavyAtoms()` - relevant if a size requirement and an R-group requirement are in the same spec.

### Radicals and unusual valences

`[CH2]` and `[N]` are radicals, not typos, and RDKit accepts them silently: `sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())` is the check.
If a spec does *not* want radicals, assert that sum is zero - writing a bracket atom to control an H count is an easy way to create one by accident, because bracket notation fixes the H count and any shortfall becomes an unpaired electron.
