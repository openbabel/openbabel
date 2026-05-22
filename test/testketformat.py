r"""Test Ketcher KET (JSON) format using the OpenBabel Python bindings.

Mirrors the structure of testcdjsonformat.py. Covers:
- reading a single-molecule KET file
- reading multi-molecule KET files
- reading a reaction KET (arrow + plus -> OBReaction)
- writing valid JSON output
- minified -xm output option
- SMILES <-> KET round-trip
- KET -> KET round-trip preserving atom/bond counts and unknown sections

Run directly while developing:
    python ../../test/testketformat.py
"""

import json
import os
import unittest

from testbindings import pybel, PybelWrapper

import faulthandler
faulthandler.enable()

filedir = os.path.join(os.path.dirname(__file__), "files", "ket")


def _readket(name):
    """Helper: read a KET file from test/ket/ as a list of pybel.Molecule."""
    return list(pybel.readfile("ket", os.path.join(filedir, name)))


def _readketstring(text, opt=None):
    """Helper: read a KET document from a string."""
    return pybel.readstring("ket", text, opt=opt)


class TestKetFormat(PybelWrapper):
    """Read/write/round-trip tests for the KET format."""

    # ------------------------------------------------------------------
    #  Read
    # ------------------------------------------------------------------

    def test_read_single_molecule(self):
        """Reads ethanol (C-C-O) and verifies atom/bond counts and structure."""
        mols = _readket("ethanol.ket")
        self.assertEqual(len(mols), 1)
        mol = mols[0]
        self.assertEqual(mol.OBMol.NumAtoms(), 3)
        self.assertEqual(mol.OBMol.NumBonds(), 2)
        self.assertEqual(
            sorted(a.atomicnum for a in mol.atoms), [6, 6, 8]
        )
        bond_elements = sorted(
            tuple(
                sorted(
                    [
                        bond.GetBeginAtom().GetAtomicNum(),
                        bond.GetEndAtom().GetAtomicNum(),
                    ]
                )
            )
            for bond in pybel.ob.OBMolBondIter(mol.OBMol)
        )
        self.assertEqual(bond_elements, [(6, 6), (6, 8)])

    def test_rejects_unsupported_future_major_version(self):
        """Future KET major versions fail explicitly instead of partial parsing."""
        text = json.dumps(
            {
                "ket_version": "3.0.0",
                "root": {"nodes": []},
            }
        )
        with self.assertRaises(OSError):
            _readketstring(text)

    def test_read_multi_component(self):
        """Reads methoxide + sodium as a single OBMol with two components."""
        mols = _readket("multi.ket")
        self.assertEqual(len(mols), 1)
        mol = mols[0]
        # 1 carbon + 1 oxygen + 1 sodium = 3 atoms.
        self.assertEqual(mol.OBMol.NumAtoms(), 3)
        atoms = list(mol.atoms)
        charges = sorted(a.formalcharge for a in atoms)
        self.assertEqual(charges, [-1, 0, 1])
        self.assertIn(11, [a.atomicnum for a in atoms])  # sodium

    def test_read_reaction(self):
        """Reads ethene + Br2 -> 1,2-dibromoethane as an OBReaction."""
        mols = _readket("reaction.ket")
        self.assertEqual(len(mols), 1)
        mol = mols[0]
        self.assertTrue(mol.OBMol.IsReaction())

        # OBReactionFacade exposes per-role component counts.
        facade = pybel.ob.OBReactionFacade(mol.OBMol)
        # 2 reactants (ethene + Br2), 1 product (1,2-dibromoethane)
        self.assertEqual(facade.NumComponents(pybel.ob.REACTANT), 2)
        self.assertEqual(facade.NumComponents(pybel.ob.PRODUCT), 1)

        # Sanity-check the product fragment has the expected element makeup.
        product = pybel.ob.OBMol()
        facade.GetComponent(product, pybel.ob.PRODUCT, 0)
        atomic_nums = sorted(
            product.GetAtom(i).GetAtomicNum()
            for i in range(1, product.NumAtoms() + 1)
        )
        self.assertEqual(atomic_nums, [6, 6, 35, 35])

    def test_read_charges(self):
        """Formal charges round-trip through the reader."""
        mols = _readket("multi.ket")
        mol = mols[0]
        # Oxygen is -1, sodium is +1, carbon is 0.
        by_elem = {a.atomicnum: a.formalcharge for a in mol.atoms}
        self.assertEqual(by_elem[8], -1)
        self.assertEqual(by_elem[11], 1)
        self.assertEqual(by_elem[6], 0)

    # ------------------------------------------------------------------
    #  Write
    # ------------------------------------------------------------------

    def test_write_is_valid_json(self):
        """Writer emits a parseable JSON document with the expected structure."""
        mols = _readket("ethanol.ket")
        text = mols[0].write("ket")
        doc = json.loads(text)
        self.assertIn("root", doc)
        self.assertIn("nodes", doc["root"])
        self.assertIn("mol0", doc)
        self.assertEqual(doc["mol0"]["type"], "molecule")
        self.assertEqual(len(doc["mol0"]["atoms"]), 3)
        self.assertEqual(len(doc["mol0"]["bonds"]), 2)

    def test_write_atom_coordinates(self):
        """All emitted atom locations are 3-element [x, y, z] arrays."""
        mols = _readket("ethanol.ket")
        doc = json.loads(mols[0].write("ket"))
        for atom in doc["mol0"]["atoms"]:
            self.assertIn("location", atom)
            self.assertEqual(len(atom["location"]), 3)

    def test_write_minified(self):
        """-xm produces no embedded newlines; default is pretty-printed."""
        mols = _readket("ethanol.ket")
        minified = mols[0].write("ket", opt={"m": None})
        self.assertNotIn("\n", minified)
        pretty = mols[0].write("ket")
        self.assertIn("\n", pretty)

    def test_write_reaction_has_arrow(self):
        """Reaction output contains exactly one arrow meta object."""
        mols = _readket("reaction.ket")
        doc = json.loads(mols[0].write("ket"))
        nodes = doc["root"]["nodes"]
        arrows = [n for n in nodes if n.get("type") == "arrow"]
        self.assertEqual(len(arrows), 1)
        # Arrow head should be to the right of arrow tail (positive x delta).
        tail, head = arrows[0]["data"]["pos"]
        self.assertGreater(head["x"], tail["x"])

    def test_write_reaction_has_pluses(self):
        """Reaction with 2 reactants emits at least one plus between them."""
        mols = _readket("reaction.ket")
        doc = json.loads(mols[0].write("ket"))
        nodes = doc["root"]["nodes"]
        pluses = [n for n in nodes if n.get("type") == "plus"]
        # 2 reactants and 1 product => at least 1 plus (between reactants).
        self.assertGreaterEqual(len(pluses), 1)

    # ------------------------------------------------------------------
    #  Round-trip
    # ------------------------------------------------------------------

    def test_smiles_to_ket_to_smiles(self):
        """SMILES -> KET -> SMILES preserves the canonical structure."""
        mol = pybel.readstring("smi", "C(=O)CO")
        mol.make2D()
        text = mol.write("ket")
        roundtrip = pybel.readstring("ket", text)
        self.assertEqual(roundtrip.write("smi").strip(), "C(=O)CO")

    def test_ket_to_ket_preserves_atom_bond_counts(self):
        """KET -> KET -> KET preserves atom/bond counts after each pass."""
        mol = _readket("ethanol.ket")[0]
        first = mol.write("ket")
        m2 = pybel.readstring("ket", first)
        second = m2.write("ket")
        m3 = pybel.readstring("ket", second)
        self.assertEqual(m2.OBMol.NumAtoms(), mol.OBMol.NumAtoms())
        self.assertEqual(m2.OBMol.NumBonds(), mol.OBMol.NumBonds())
        self.assertEqual(m3.OBMol.NumAtoms(), mol.OBMol.NumAtoms())
        self.assertEqual(m3.OBMol.NumBonds(), mol.OBMol.NumBonds())
        # JSON structure remains stable across re-serialization.
        self.assertEqual(json.loads(first)["mol0"]["atoms"],
                         json.loads(second)["mol0"]["atoms"])

    def test_round_trip_preserves_charges(self):
        """Round-trip through KET preserves per-atom formal charges."""
        mols = _readket("multi.ket")
        text = mols[0].write("ket")
        rt = pybel.readstring("ket", text)
        rt_charges = sorted(a.formalcharge for a in rt.atoms)
        self.assertEqual(rt_charges, [-1, 0, 1])

    def test_round_trip_preserves_explicit_zero_implicit_h_count(self):
        """Explicit implicitHCount: 0 survives a KET round-trip."""
        text = json.dumps(
            {
                "root": {"nodes": [{"$ref": "mol0"}]},
                "mol0": {
                    "type": "molecule",
                    "atoms": [
                        {
                            "label": "N",
                            "location": [0.0, 0.0, 0.0],
                            "implicitHCount": 0,
                        }
                    ],
                    "bonds": [],
                },
            }
        )
        doc = json.loads(_readketstring(text).write("ket"))
        self.assertEqual(doc["mol0"]["atoms"][0]["implicitHCount"], 0)

    def test_atom_alias_is_available_as_alias_data(self):
        """KET atom aliases are exposed through Open Babel AliasData."""
        text = json.dumps(
            {
                "root": {"nodes": [{"$ref": "mol0"}]},
                "mol0": {
                    "type": "molecule",
                    "atoms": [
                        {
                            "label": "*",
                            "alias": "COOH",
                            "location": [0.0, 0.0, 0.0],
                        }
                    ],
                    "bonds": [],
                },
            }
        )
        atom = _readketstring(text).atoms[0].OBAtom
        data = atom.GetData(pybel.ob.AliasDataType)
        self.assertTrue(data)
        self.assertEqual(pybel.ob.toAliasData(data).GetAlias(), "COOH")

    def test_atom_alias_expands_with_read_option(self):
        """The -ia read option expands chemically meaningful KET atom aliases."""
        text = json.dumps(
            {
                "root": {"nodes": [{"$ref": "mol0"}]},
                "mol0": {
                    "type": "molecule",
                    "atoms": [
                        {
                            "label": "C",
                            "location": [0.0, 0.0, 0.0],
                        },
                        {
                            "label": "*",
                            "alias": "COOH",
                            "location": [1.0, 0.0, 0.0],
                        }
                    ],
                    "bonds": [{"type": 1, "atoms": [0, 1]}],
                },
            }
        )
        mol = _readketstring(text, opt={"a": None})
        self.assertEqual(mol.OBMol.NumAtoms(), 4)
        self.assertEqual(mol.OBMol.NumBonds(), 3)
        self.assertEqual(sorted(a.atomicnum for a in mol.atoms), [6, 6, 8, 8])
        data = mol.atoms[1].OBAtom.GetData(pybel.ob.AliasDataType)
        self.assertTrue(data)
        alias = pybel.ob.toAliasData(data)
        self.assertTrue(alias.IsExpanded())
        self.assertEqual(alias.GetAlias(), "COOH")

    def test_round_trip_preserves_reaction(self):
        """A reaction KET stays a reaction after a round-trip."""
        mol = _readket("reaction.ket")[0]
        text = mol.write("ket")
        rt = pybel.readstring("ket", text)
        self.assertTrue(rt.OBMol.IsReaction())
        facade = pybel.ob.OBReactionFacade(rt.OBMol)
        self.assertEqual(facade.NumComponents(pybel.ob.REACTANT), 2)
        self.assertEqual(facade.NumComponents(pybel.ob.PRODUCT), 1)

    def test_round_trip_preserves_connections(self):
        """root.connections is preserved verbatim across a round-trip."""
        mol = _readket("with_connections.ket")[0]
        text = mol.write("ket")
        doc = json.loads(text)
        self.assertIn("connections", doc["root"])
        self.assertEqual(len(doc["root"]["connections"]), 1)
        conn = doc["root"]["connections"][0]
        self.assertEqual(conn["connectionType"], "single")
        self.assertEqual(conn["endpoint1"]["monomerId"], "monomer0")
        self.assertEqual(conn["endpoint1"]["attachmentPointId"], "R1")
        self.assertEqual(conn["endpoint2"]["moleculeId"], "mol5")

    def test_round_trip_preserves_ref_name(self):
        """Original $ref name (mol5, not mol0) survives a round-trip."""
        mol = _readket("with_connections.ket")[0]
        doc = json.loads(mol.write("ket"))
        refs = [n["$ref"] for n in doc["root"]["nodes"] if "$ref" in n]
        self.assertIn("mol5", refs)
        self.assertIn("mol5", doc)

    def test_round_trip_preserves_node_order(self):
        """The text meta object stays between mol5 and monomer0."""
        mol = _readket("with_connections.ket")[0]
        doc = json.loads(mol.write("ket"))
        nodes = doc["root"]["nodes"]
        kinds = [
            ("ref", n["$ref"]) if "$ref" in n else ("type", n.get("type"))
            for n in nodes
        ]
        self.assertEqual(
            kinds, [("ref", "mol5"), ("type", "text"), ("ref", "monomer0")]
        )

    def test_sgroup_bonds_remain_local_per_molecule(self):
        """S-group bond indices are local to their owning molecule."""
        mol = _readket("full_fidelity.ket")[0]
        doc = json.loads(mol.write("ket"))
        # mol0: 3 bonds, S-group references bonds 0 and 2.
        self.assertIn("sgroups", doc["mol0"])
        gen0 = next(s for s in doc["mol0"]["sgroups"] if s["type"] == "GEN")
        self.assertEqual(gen0["bonds"], [0, 2])
        # mol1: 1 bond, S-group references bond 0 (local!), NOT 3 (global).
        gen1 = next(s for s in doc["mol1"]["sgroups"] if s["type"] == "GEN")
        self.assertEqual(gen1["bonds"], [0])

    def test_per_component_highlight_and_properties_scope(self):
        """Highlight and properties don't leak across components."""
        mol = _readket("full_fidelity.ket")[0]
        doc = json.loads(mol.write("ket"))
        self.assertIn("highlight", doc["mol0"])
        self.assertIn("properties", doc["mol0"])
        self.assertNotIn("highlight", doc["mol1"])
        self.assertNotIn("properties", doc["mol1"])

    def test_no_duplicate_arrow_on_reaction_round_trip(self):
        """Original arrow is preserved; no extra synthesized one."""
        mol = _readket("reaction.ket")[0]
        doc = json.loads(mol.write("ket"))
        arrows = [n for n in doc["root"]["nodes"] if n.get("type") == "arrow"]
        self.assertEqual(len(arrows), 1)
        # Original arrow tail at x=2, head at x=5 (not the synthesized values).
        tail, head = arrows[0]["data"]["pos"]
        self.assertAlmostEqual(tail["x"], 2.0, places=2)
        self.assertAlmostEqual(head["x"], 5.0, places=2)

    def test_empty_obmol_emits_valid_empty_document(self):
        """Writing an empty OBMol does not produce dangling refs."""
        mol = pybel.ob.OBMol()
        text = pybel.Molecule(mol).write("ket")
        doc = json.loads(text)
        self.assertEqual(doc["root"]["nodes"], [])
        # No mol* members.
        mol_members = [k for k in doc if k.startswith("mol")]
        self.assertEqual(mol_members, [])

    def test_round_trip_preserves_passthrough(self):
        """Unknown sections (monomer / monomerTemplate / image) survive a round-trip."""
        mol = _readket("with_passthrough.ket")[0]
        text = mol.write("ket")
        doc = json.loads(text)

        # The monomer ref and its top-level object must reappear.
        node_refs = [n.get("$ref") for n in doc["root"]["nodes"] if "$ref" in n]
        self.assertIn("monomer0", node_refs)
        self.assertIn("monomer0", doc)
        self.assertEqual(doc["monomer0"]["type"], "monomer")

        # The template ref + top-level object must reappear under root.templates.
        self.assertIn("templates", doc["root"])
        template_refs = [t.get("$ref") for t in doc["root"]["templates"]]
        self.assertIn("monomerTemplate-A___Adenine", template_refs)
        self.assertIn("monomerTemplate-A___Adenine", doc)

        # The inline image meta object survives.
        images = [n for n in doc["root"]["nodes"] if n.get("type") == "image"]
        self.assertEqual(len(images), 1)
        self.assertEqual(images[0]["format"], "image/png")

        # ket_version 2.0.0 round-trips since the input requested it.
        self.assertEqual(doc.get("ket_version"), "2.0.0")


if __name__ == "__main__":
    unittest.main()
