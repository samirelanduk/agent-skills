---
name: smiles-generator
description: >-
  Construct a SMILES or CXSMILES string that provably satisfies a stated structural specification, verified in RDKit before it is presented.
  Use whenever someone asks for a molecule described by its properties rather than its name - a number of chiral centres, specified vs unspecified stereo, stereogenic vs non-stereogenic double bonds, permanent vs removable charges, a target net charge, N fragments, a salt, solvate, tautomer pair, meso form, enhanced stereo (AND/OR/ABS) groups, a wanted stereoisomer count, or stereo
  that depends on protonation state.
  Use it even when RDKit is never mentioned: verification is part of the deliverable.
---

# Building a SMILES to a structural specification

Someone hands you a list of structural properties, and you return one SMILES string that demonstrably has all of them. That is the whole job.

**Verify before you present.**
A string that looks right routinely has an extra stereocentre, a stereo group RDKit dropped without complaint, or a charge that evaporates.
None of that shows up by reading.
So the check has to happen before the answer exists, not after someone asks for it - never 'offer' to verify, just verify.

**Search in private, deliver in public.**
Build, test, fix and re-test inside a script.
What you show is one string plus a checklist.
A visible tour of the candidates you rejected is unneccessary.

## Step 0 - get RDKit running

```bash
python3 -c "import rdkit; print(rdkit.__version__)" || pip install rdkit --break-system-packages -q
```

Install `rdkit`. Not `rdkit-pypi` - that name resolves to nothing ("Could not find a version that satisfies the requirement").

Read the install result before doing anything else. Firing the install and then reasoning in prose while it fails in the background is how an entire answer ends up unverified.

## Step 1 - turn the spec into predicates

Before touching a SMILES, rewrite the request as one testable predicate per requirement.
You will re-run all of them after every edit, so they need to exist as executable data rather than as a memory of what was asked.

Concretely, that means a `checks` list - one entry per requirement, each a label and a boolean, ending in a single overall verdict.
Step 3 has the template.
The list is short, it is written fresh for each spec, and it is the thing you re-run rather than re-reading the request and trusting your memory of it.

Anything RDKit can compute is a legitimate predicate, so the spec is never limited by a fixed vocabulary of properties:

```python
len(Chem.GetMolFrags(mol)) == 2
Descriptors.MolWt(mol) < 400
Descriptors.NumRotatableBonds(mol) == 0
len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SX4](=O)(=O)[NX3]'))) == 1
sorted(len(r) for r in mol.GetRingInfo().AtomRings()) == [5, 6, 6]
any(a.GetIsotope() == 13 for a in mol.GetAtoms())
```

Write the predicate for what was actually asked, not for the nearest thing that
is convenient to measure.
And write it fresh: a checker carried over from an earlier turn, still asserting that turn's thresholds, reports `False` on a molecule that is completely correct, and then you are explaining your own tooling instead of the answer.

When a requirement resists expression as a predicate, that is usually a signal that you haven't pinned down what it means yet - which is Step 1's real work.
"Three fragments" is unambiguous; "drug-like", "a realistic salt", "not too strained" are not, and deciding what they mean is cheaper now than after you've built to the wrong reading.

Two important points:

**Follow-up turns are additive.**
"Now add a stereogenic double bond" means *in addition to* everything already agreed.
Edit the existing molecule and re-run the full list.
Rebuilding from scratch is how a satisfied requirement gets quietly dropped while the new one is satisfied - which is the most irritating possible way to fail, because the reader has to re-check work they had already accepted.

**Don't widen the spec.**
If net charge was never mentioned, it is not a constraint.
Inventing one, then reporting at length on why the invented version is hard, and spends the reader's time on a problem they don't have.

### Terms that carry more weight than they look like they do

Most of a spec reads literally.
The failures cluster on a handful of words that sound descriptive but are being used as precise technical requirements, where the everyday reading and the intended one produce different molecules.
When one of these appears, the molecule can be invalid even though the string is perfectly valid and every other requirement is met.

The list below is the set that comes up often, not a complete inventory - the general habit is what matters: when a word could be read two ways and the readings differ structurally, settle it before building.

- **Removable charge** - neutralised by proton transfer alone: the neutral form differs by one H and nothing else. `[NH3+]`→`NH2`, `[O-]`→`OH`.
- **Permanent charge** - no acid–base equilibrium at all.
Not "an extreme pKa": neutralising it would require making or breaking a bond to something other than hydrogen.
Quaternary ammonium and tetrasubstituted borate qualify.
**Sulfonate and carboxylate do not**, however low the pKa - and a caveat saying sulfonate is "effectively permanent at any reasonable pH" does not rescue it, it just delays the rejection.
This is machine-checkable; see Step 3.
- **Unspecified stereo** - a genuine stereogenic element carrying no marker.
Stereogenicity is a property of the graph, not of the string, so deleting `@` from a stereocentre makes the molecule *ambiguous*, not achiral.
If someone asks for a molecule that "completely specifies a structure" or a "non-isomeric" one, they mean zero unspecified elements - which is a much stronger requirement than zero `@` characters.
- **Non-stereogenic double bond** - a C=C or C=O that *cannot* express E/Z because one end carries two identical substituents.
Removing `/` and `\` from a stereogenic bond produces an unspecified one, not a non-stereogenic one.
- **"Stereocentres"** is often used loosely for the whole stereo inventory, atoms *and* bonds.
Report both, and say which is which.
- **Tautomer pair** - a charge-conserving proton shift.
Deprotonation (`CC(=O)O` → `CC(=O)[O-]`) is acid–base chemistry, not tautomerism.
- **Single fragment** means no `.` anywhere.
Don't reach for a two-part salt to dodge a valence or stereo problem.
- **Charged** is ambiguous between "contains charged atoms" and "has non-zero net charge" - a zwitterion satisfies the first and not the second.
Check which is meant before building, because the two lead to different molecules.

Where two readings would produce materially different molecules and the cost of guessing wrong is a wasted build, ask.
Otherwise decide, state the convention in one line, and build - don't ask about something you could have resolved yourself.

## Step 2 - build from motifs that behave predictably

Pick pieces whose behaviour follows from symmetry.
`references/motifs.md` is a catalogue of verified building blocks keyed by the property each one supplies.
It covers the axes that come often - charges, stereo elements, substructure and ring requirements, size and property targets - but it is a starting set, not a menu of what can be asked for.
When a spec wants something not in it, the method is the same: build the smallest fragment that has the property, confirm it in RDKit on its own, then attach it.

The assembly pattern:

1. **Settle the global quantities before building anything local.**
Anything that must sum to a target across the whole molecule - net charge, fragment count, an element tally, a heavy-atom or MW budget - is cheap to solve on paper and expensive to fix by editing a built structure.
Charge is the usual case: permanent total must equal the post-neutralisation charge, and removable total = target − permanent.
2. **Verify each piece standalone**, batched into one script rather than one at a time.
Testing ten candidates costs the same as testing one, and the batch is what turns "which of these behaves?" from a guess into a table.
3. **Hang them off an inert hub.**
A benzene ring with CH2 linkers is the reliable default: aromatic ring carbons are sp2 and a benzylic CH2 carries two H, so neither becomes an accidental stereocentre no matter what you attach.
The general requirement for a hub is that it contributes nothing the spec counts - pick one whose own properties are zero on every axis in play.
4. **Verify the assembled whole.**
Assembly creates properties the pieces didn't have: the carbon you attached three different arms to is now a stereocentre, two fragments that were fine apart may now share a formula, and a ring you closed for one requirement counts against a ring-count target.

Aim for the least contrived molecule that meets the letter of the spec.
The degenerate solution to a counting requirement is usually available and usually
unwelcome: three water molecules for "three fragments", `[NH4+].[Cl-]` for "a
salt", a methyl for "a substituent".
Satisfying the count while ignoring what the count was for reads as not having engaged with the request.
Unless the spec sets a floor, keep every non-solvent fragment substantive - around 9 heavy atoms is a comfortable minimum - and when several pieces must be "different", make them differ in kind (heteroatom set, branching topology, ring system), not just in where one substituent sits.

Keep it small enough to depict.
These strings usually get pasted into a structure renderer, so a compact alkyl solution beats a bulky fused-ring one carrying the same properties.

## Step 3 - verify

Write a short script per candidate.
This template carries the measurements whose obvious API call gives a confidently wrong answer - swap the `checks` list for your spec and leave the rest:

```python
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

smi = '<candidate>'

mol = Chem.MolFromSmiles(smi)          # default removeHs=True: the Uncharger needs it
assert mol is not None, 'parse failed - RDKit wrote the reason to stderr'
rdCIPLabeler.AssignCIPLabels(mol)      # before reading any _CIPCode

# Stereo inventory. FindPotentialStereo, not GetStereo/FindMolChiralCenters:
# those report only what is already marked. centeredOn is an ATOM index for
# Atom_*, a BOND index for Bond_* - never look one up in the other's table.
elems  = list(Chem.FindPotentialStereo(mol))
s_atom = [e for e in elems if 'Atom' in str(e.type)]
s_bond = [e for e in elems if 'Bond' in str(e.type)]
n = lambda seq, state: sum(1 for e in seq if str(e.specified) == state)

# Removable vs permanent charge. protonationOnly restricts it to proton transfer;
# force makes it per-atom (without it, charges are left in place to balance).
# Re-parse: in-place edits leave stale @ tags that the writer will emit.
unc       = rdMolStandardize.Uncharger(protonationOnly=True, force=True)
neutral   = Chem.MolFromSmiles(Chem.MolToSmiles(unc.uncharge(Chem.Mol(mol))))
permanent = [a.GetIdx() for a in neutral.GetAtoms() if a.GetFormalCharge()]

# Non-stereogenic double bonds have no API - subtract the stereogenic ones.
marked    = {e.centeredOn for e in s_bond}
nonstereo = [b.GetIdx() for b in mol.GetBonds()
             if b.GetBondTypeAsDouble() == 2.0 and not b.GetIsAromatic()
             and b.GetIdx() not in marked]

# maxIsomers=0 or the count silently truncates at 1024 - a plausible true answer.
isomers = len(list(EnumerateStereoisomers(
    mol, StereoEnumerationOptions(maxIsomers=0, unique=True))))

checks = [                                     # one line per stated requirement
    ('single fragment',            len(Chem.GetMolFrags(mol)) == 1),
    ('net charge 0',               Chem.GetFormalCharge(mol) == 0),
    ('1 permanent charge',         len(permanent) == 1),
    ('1 unspecified stereocentre', n(s_atom, 'Unspecified') == 1),
    ('no specified stereo',        n(s_atom, 'Specified') + n(s_bond, 'Specified') == 0),
    ('MW under 400',               Descriptors.MolWt(mol) < 400),
]
for label, ok in checks:
    print(('[OK] ' if ok else '[X]  ') + label)
print('ALL REQUIREMENTS MET:', all(ok for _, ok in checks))
```

Print the inventory alongside it - `elems`, the fragment charges, the neutralised form - even when every check passes.

Four of those measurements are worth understanding rather than copying:

**Parse, and test for `None`.**
Valence errors go to *stderr* while `MolFromSmiles` returns `None`.
A generation loop that doesn't test for `None` produces a clean run and no molecules, which looks like success.

**`FindPotentialStereo` for the stereo inventory** - not `bond.GetStereo()`, not bare `FindMolChiralCenters`.
Those report only elements that already carry a marker, so they can confirm what you wrote but cannot tell you about the centre you forgot.
`FindPotentialStereo` returns potential elements whether specified or not, which is exactly the question "does this molecule have stereo I'm not marking?"
It also catches the elements nobody looks for: an amidine or guanidinium C=N is stereogenic (`CC(=N)N`, four heavy atoms, is the minimal case), and a ring double bond becomes stereogenic once the ring reaches 8 atoms.

**The Uncharger for removable vs permanent.**

```python
from rdkit.Chem.MolStandardize import rdMolStandardize
unc = rdMolStandardize.Uncharger(protonationOnly=True, force=True)
neutral = unc.uncharge(Chem.Mol(mol))
```

`protonationOnly=True` restricts it to proton transfer, which is precisely the
removable/permanent distinction - without it the Uncharger does things proton
transfer cannot, reducing `CC[BH3-]` to neutral boron.
`force=True` makes it a *per-atom* test: by default it preserves overall charge balance, so a betaine `C[N+](C)(C)CCC(=O)[O-]` comes back untouched and its perfectly removable carboxylate reads as permanent.

Parse with the default `removeHs=True`.
The Uncharger edits hydrogen *counts* and cannot delete a hydrogen *node*, so on a molecule parsed with `removeHs=False` it is a silent no-op and every charge looks permanent.

**The counterfactual, for any "X is stereogenic only because of Y" claim.**
Build the variant without Y, re-count, and confirm exactly one element disappeared and nothing else moved.
An argument is not evidence:

```
OC(=O)CC(CC(=O)[O-])C   → 1 stereocentre   (as written)
CC(CC(=O)O)CC(=O)O      → 0 stereocentres  (protonated: the two arms become identical)
```

For a *removable*-charge dependency the `neutral` mol from the template is already the counterfactual - run `FindPotentialStereo` on it as well as on `mol` and compare.
For a *permanent*-charge dependency the Uncharger by definition won't touch it, so write the hypothetical neutral analogue by hand and count that instead.

### Final checks

**Re-parse the string you are about to deliver.**
After any in-place edit - uncharging, setting tags, clearing directions - the writer will happily emit a `@` on an atom that is no longer a stereocentre.
Round-tripping the final string through `MolFromSmiles` and re-running the checks is the only way to know that what you paste is what you measured.

**Check the neutralised form's stereo too, not just the submitted one.**
A `@` written on an atom that isn't chiral in the state you submitted is silently stripped, and neutralisation can create an element that wasn't there: the default Uncharger protonates just *one* of two symmetric removable arms, turning a symmetric carbon into an unspecified stereocentre that the molecule as written did not have.

## Step 4 - present it

One string, then a compact numbered checklist mapping each stated requirement to the observed value that confirms it, with explicit OK marks.
Report the numbers RDKit produced - atom indices, CIP labels, per-fragment charges.

Separate what the toolkit measured from what you are asserting chemically.
RDKit confirms that a quaternary ammonium survives the Uncharger; that it is chemically non-ionisable is your reasoning, and saying which is which is welcomed rather than treated as hedging.
Likewise a COO⁻/COOH twin-arm stereocentre is real in the static connectivity graph while the proton exchanges in solution - worth one line.

Three things to keep out of the visible answer:

- **Alternatives and abandoned candidates.**
Deciding between eleven strings in prose is a worse deliverable than one string, even if the eleventh is right.
- **Mid-answer self-correction.**
If you change your mind, redo it silently and present the clean version.
- **Hedges that mean the answer is wrong.**
If a requirement can't be met, say so plainly and say why.
A nearby structure presented as though it were the requested one - an ABS group where an unknown group was asked for, a sulfonate where a permanent negative was asked for - is a failure.

If part of the spec turns out to be unsatisfiable, deliver everything else and isolate the blocked part.
A working partial answer with a clear note beats a long silence followed by nothing.

And if a check contradicts something you said earlier, re-run it rather than folding to whichever position was stated most recently.
"I haven't established that; here is what I did establish" is a legitimate answer.
Reversing under pushback without new evidence is not, and neither is over-correcting back.

## Reference files

- **`references/motifs.md`** - verified construction motifs keyed by the requirement each satisfies. Read in Step 2, before building.
- **`references/rdkit-notes.md`** - API behaviours that affect the output: stereo perception, CIP labelling, canonical ordering, CXSMILES and enhanced stereo groups, stereoisomer enumeration, the Uncharger's parameters, and the SMILES syntax mistakes that cause valence errors.
Read it before writing any CXSMILES with stereo groups, or when a check behaves unexpectedly.
- **`references/salts-solvates.md`** - counter-ion and solvate SMILES with measured charge, heavy-atom count and stereo status.
Use it for salt/solvate specs instead of transcribing a structure from memory.
