# Counter-ions and solvates

Counter-ion and solvate SMILES with measured charge, heavy-atom count and stereo status.
Use these instead of writing one from memory.

**Charge** is the formal charge as written.
**Unc** is after `Uncharger(protonationOnly=True, force=True)` - `→0` means removable, unchanged means permanent.

Same structure can sit in both lists: neutral form as solvate, charged form as counter-ion (acetic acid / acetate, water / hydroxide).

---

## Solvates (neutral, use with `.` )

| Name | SMILES | Heavy | Note |
|---|---|---|---|
| water | `O` | 1 | the only fragment normally exempt from a "no tiny fragments" rule |
| methanol | `CO` | 2 | |
| ethanol | `CCO` | 3 | |
| 1-propanol | `CCCO` | 4 | |
| 2-propanol | `CC(C)O` | 4 | |
| 1-butanol | `CCCCO` | 5 | |
| *tert*-butanol | `CC(C)(C)O` | 5 | |
| ethylene glycol | `OCCO` | 4 | |
| acetone | `CC(C)=O` | 4 | 1 non-stereogenic C=O |
| methyl ethyl ketone | `CCC(C)=O` | 5 | |
| ethyl acetate | `CCOC(C)=O` | 6 | |
| methyl acetate | `COC(C)=O` | 5 | |
| isopropyl acetate | `CC(C)OC(C)=O` | 7 | |
| diethyl ether | `CCOCC` | 5 | |
| MTBE | `COC(C)(C)C` | 6 | |
| THF | `C1CCOC1` | 5 | |
| 2-methyl-THF | `CC1CCCO1` | 6 | **1 unspecified stereocentre** |
| 1,4-dioxane | `C1COCCO1` | 6 | |
| anisole | `COc1ccccc1` | 8 | |
| dichloromethane | `ClCCl` | 3 | |
| chloroform | `ClC(Cl)Cl` | 4 | |
| 1,2-dichloroethane | `ClCCCl` | 4 | |
| carbon tetrachloride | `ClC(Cl)(Cl)Cl` | 5 | |
| DMF | `CN(C)C=O` | 5 | |
| dimethylacetamide | `CC(=O)N(C)C` | 6 | |
| NMP | `CN1CCCC1=O` | 7 | |
| DMSO | `CS(C)=O` | 4 | |
| sulfolane | `O=S1(=O)CCCC1` | 7 | |
| acetonitrile | `CC#N` | 3 | |
| nitromethane | `C[N+](=O)[O-]` | 4 | charges present, net 0 |
| pyridine | `c1ccncc1` | 6 | |
| triethylamine | `CCN(CC)CC` | 7 | |
| benzene | `c1ccccc1` | 6 | |
| toluene | `Cc1ccccc1` | 7 | |
| *o*-xylene | `Cc1ccccc1C` | 8 | |
| chlorobenzene | `Clc1ccccc1` | 7 | |
| *n*-pentane | `CCCCC` | 5 | |
| *n*-hexane | `CCCCCC` | 6 | |
| *n*-heptane | `CCCCCCC` | 7 | |
| cyclohexane | `C1CCCCC1` | 6 | |
| acetic acid | `CC(=O)O` | 4 | |
| formic acid | `O=CO` | 3 | |
| trifluoroacetic acid | `OC(=O)C(F)(F)F` | 7 | |

---

## Cationic counter-ions

| Name | SMILES | Charge | Unc | Heavy |
|---|---|---|---|---|
| sodium | `[Na+]` | +1 | +1 | 1 |
| potassium | `[K+]` | +1 | +1 | 1 |
| lithium | `[Li+]` | +1 | +1 | 1 |
| calcium | `[Ca+2]` | +2 | +2 | 1 |
| magnesium | `[Mg+2]` | +2 | +2 | 1 |
| zinc | `[Zn+2]` | +2 | +2 | 1 |
| aluminium | `[Al+3]` | +3 | +3 | 1 |
| ammonium | `[NH4+]` | +1 | →0 | 1 |
| choline | `C[N+](C)(C)CCO` | +1 | +1 | 7 |
| tetrabutylammonium | `CCCC[N+](CCCC)(CCCC)CCCC` | +1 | +1 | 17 |

Metal cations survive the Uncharger and look permanent under that test, but they are bare ions.
For a permanent charge on a molecule, use quaternary ammonium or borate (`motifs.md`).

**Amines as cations.**
Write them protonated:

| Name | Free base | Cation form | Charge |
|---|---|---|---|
| diethylamine | `CCNCC` | `CC[NH2+]CC` | +1 |
| tromethamine | `OCC(N)(CO)CO` | `[NH3+]C(CO)(CO)CO` | +1 |
| meglumine | `CNCC(O)C(O)C(O)C(O)CO` | `C[NH2+]CC(O)C(O)C(O)C(O)CO` | +1 |
| benzathine | `C(CNCc1ccccc1)NCc1ccccc1` | protonate either N | +1 / +2 |

Meglumine as written has **4 unspecified stereocentres**.
The real compound comes from D-glucose - add the tags if the spec needs them.

---

## Anionic counter-ions

| Name | SMILES | Charge | Unc | Heavy | Stereo |
|---|---|---|---|---|---|
| chloride | `[Cl-]` | −1 | →0 | 1 | |
| bromide | `[Br-]` | −1 | →0 | 1 | |
| iodide | `[I-]` | −1 | →0 | 1 | |
| fluoride | `[F-]` | −1 | →0 | 1 | |
| sulfate | `[O-]S(=O)(=O)[O-]` | −2 | →0 | 5 | |
| hydrogen sulfate | `OS(=O)(=O)[O-]` | −1 | →0 | 5 | |
| phosphate | `[O-]P(=O)([O-])[O-]` | −3 | →0 | 5 | |
| dihydrogen phosphate | `OP(=O)(O)[O-]` | −1 | →0 | 5 | |
| nitrate | `[O-][N+](=O)[O-]` | −1 | →0 | 4 | |
| carbonate | `[O-]C(=O)[O-]` | −2 | →0 | 4 | |
| bicarbonate | `OC(=O)[O-]` | −1 | →0 | 4 | |
| mesylate | `CS(=O)(=O)[O-]` | −1 | →0 | 5 | |
| esylate | `CCS(=O)(=O)[O-]` | −1 | →0 | 6 | |
| besylate | `[O-]S(=O)(=O)c1ccccc1` | −1 | →0 | 10 | |
| tosylate | `Cc1ccc(cc1)S(=O)(=O)[O-]` | −1 | →0 | 11 | |
| napsylate | `[O-]S(=O)(=O)c1ccc2ccccc2c1` | −1 | →0 | 14 | |
| isethionate | `OCCS(=O)(=O)[O-]` | −1 | →0 | 7 | |
| edisylate | `[O-]S(=O)(=O)CCS(=O)(=O)[O-]` | −2 | →0 | 10 | |
| camsylate | `CC1(C)[C@@H]2CC[C@]1(CS(=O)(=O)[O-])C(=O)C2` | −1 | →0 | 15 | 2 specified (RDKit calls this pattern 3R,6R) |
| triflate | `[O-]S(=O)(=O)C(F)(F)F` | −1 | →0 | 8 | |
| acetate | `CC(=O)[O-]` | −1 | →0 | 4 | |
| formate | `[O-]C=O` | −1 | →0 | 3 | |
| trifluoroacetate | `[O-]C(=O)C(F)(F)F` | −1 | →0 | 7 | |
| benzoate | `[O-]C(=O)c1ccccc1` | −1 | →0 | 9 | |
| salicylate | `Oc1ccccc1C(=O)[O-]` | −1 | →0 | 10 | |
| (S)-lactate | `C[C@H](O)C(=O)[O-]` | −1 | →0 | 6 | 1 specified |
| gluconate | `OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]` | −1 | →0 | 13 | 4 specified |
| stearate | `CCCCCCCCCCCCCCCCCC(=O)[O-]` | −1 | →0 | 20 | |
| oxalate | `[O-]C(=O)C(=O)[O-]` | −2 | →0 | 6 | |
| malonate | `[O-]C(=O)CC(=O)[O-]` | −2 | →0 | 7 | |
| succinate | `[O-]C(=O)CCC(=O)[O-]` | −2 | →0 | 8 | |
| adipate | `[O-]C(=O)CCCCC(=O)[O-]` | −2 | →0 | 10 | |
| fumarate | `[O-]C(=O)/C=C/C(=O)[O-]` | −2 | →0 | 8 | 1 specified bond (E) |
| maleate | `[O-]C(=O)/C=C\C(=O)[O-]` | −2 | →0 | 8 | 1 specified bond (Z) |
| (S)-malate | `[O-]C(=O)C[C@H](O)C(=O)[O-]` | −2 | →0 | 9 | 1 specified (S) |
| L-(+)-tartrate | `[O-]C(=O)[C@H](O)[C@@H](O)C(=O)[O-]` | −2 | →0 | 10 | 2 specified, (R,R) |
| citrate | `[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]` | −3 | →0 | 13 | none - the central C has two identical arms |
| embonate (pamoate) | `Oc1c(Cc2c(O)c3ccccc3cc2C(=O)[O-])cc2ccccc2c1C(=O)[O-]` | −2 | →0 | 29 | |
| tetrafluoroborate | `[B-](F)(F)(F)F` | −1 | **−1** | 5 | **permanent** |
| tetraphenylborate | `[B-](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1` | −1 | **−1** | 25 | **permanent** |
| hexafluorophosphate | `F[P-](F)(F)(F)(F)F` | −1 | →0 | 7 | Uncharger neutralises it; RDKit also flags an unspecified stereo element on P |

Every anion here is **removable** under the Uncharger except the two borates.
That includes the sulfonates: sulfonate does not satisfy a permanent-negative-charge spec, however low its pKa.

Tartrate: `[C@H](O)[C@@H](O)` - opposite tags give (R,R); matching tags give meso.
See meso in `motifs.md`.

Named configs (L-tartrate, (S)-malate, (S)-lactate, D-gluconate, camsylate) are the conventional ones.
RDKit checks the CIP on the tags you wrote; it does not check that the name maps to that enantiomer.
