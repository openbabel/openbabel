/*
 * AromaticityResolver - Port of OpenChemLib's RingCollection aromaticity detection
 *
 * This class provides Java-compatible aromaticity detection for molecular structures.
 * It uses the same pattern-matching algorithm as OpenChemLib's RingCollection class
 * to ensure consistent CLogP predictions between Java and C++ implementations.
 *
 * Original Java source: com.actelion.research.chem.RingCollection
 * Copyright (c) 1997 - 2016 Actelion Pharmaceuticals Ltd.
 *
 * C++ port maintains the same BSD-style license.
 */

#ifndef AROMATICITY_RESOLVER_H
#define AROMATICITY_RESOLVER_H

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <vector>
#include <unordered_set>

namespace OpenChemLib {

/**
 * AromaticityResolver provides Java-compatible aromaticity detection.
 *
 * This class analyzes a molecule's ring systems and determines aromaticity
 * using the same pattern-matching rules as OpenChemLib, ensuring consistent
 * fingerprint generation and CLogP predictions.
 *
 * Usage:
 *   AromaticityResolver resolver(mol);
 *   if (resolver.isAromaticBond(bond)) {
 *       // Bond is in an aromatic ring
 *   }
 */
class AromaticityResolver {
public:
    // Maximum ring size considered a "small ring" (matches Java)
    static constexpr int MAX_SMALL_RING_SIZE = 7;
    static constexpr int MAX_SMALL_RING_COUNT = 1024;

    /**
     * Constructor - analyzes molecule and computes aromaticity.
     * This performs the full ring detection and aromaticity analysis.
     *
     * @param mol OpenBabel molecule to analyze
     */
    explicit AromaticityResolver(OpenBabel::OBMol& mol);

    /**
     * Destructor
     */
    ~AromaticityResolver();

    // ========================================================================
    // Bond-level queries (primary API for fingerprinting)
    // ========================================================================

    /**
     * Check if a bond is aromatic (member of an aromatic ring).
     * Uses Java's pattern-matching algorithm.
     *
     * @param bond OpenBabel bond pointer
     * @return true if bond is in an aromatic ring
     */
    bool isAromaticBond(const OpenBabel::OBBond* bond) const;

    /**
     * Check if a bond is delocalized.
     * Delocalized bonds are aromatic bonds in 6-membered rings.
     * This is a subset of aromatic bonds.
     *
     * @param bond OpenBabel bond pointer
     * @return true if bond is delocalized
     */
    bool isDelocalizedBond(const OpenBabel::OBBond* bond) const;

    /**
     * Check if a bond is heteroaromatic.
     *
     * @param bond OpenBabel bond pointer
     * @return true if bond is in a heteroaromatic ring
     */
    bool isHeteroAromaticBond(const OpenBabel::OBBond* bond) const;

    // ========================================================================
    // Atom-level queries
    // ========================================================================

    /**
     * Check if an atom is aromatic (member of an aromatic ring).
     *
     * @param atom OpenBabel atom pointer
     * @return true if atom is in an aromatic ring
     */
    bool isAromaticAtom(const OpenBabel::OBAtom* atom) const;

    /**
     * Check if an atom is delocalized.
     *
     * @param atom OpenBabel atom pointer
     * @return true if atom is delocalized
     */
    bool isDelocalizedAtom(const OpenBabel::OBAtom* atom) const;

    /**
     * Check if an atom is heteroaromatic.
     *
     * @param atom OpenBabel atom pointer
     * @return true if atom is heteroaromatic
     */
    bool isHeteroAromaticAtom(const OpenBabel::OBAtom* atom) const;

    // ========================================================================
    // Ring-level queries (for debugging and validation)
    // ========================================================================

    /**
     * Get the number of small rings detected.
     *
     * @return number of rings (up to MAX_SMALL_RING_SIZE)
     */
    int getRingCount() const;

    /**
     * Check if a specific ring is aromatic.
     *
     * @param ringNo ring index (0-based)
     * @return true if ring is aromatic
     */
    bool isAromaticRing(int ringNo) const;

    /**
     * Check if a specific ring is delocalized.
     *
     * @param ringNo ring index (0-based)
     * @return true if ring is delocalized
     */
    bool isDelocalizedRing(int ringNo) const;

    /**
     * Get the size (number of atoms) of a ring.
     *
     * @param ringNo ring index (0-based)
     * @return number of atoms in ring
     */
    int getRingSize(int ringNo) const;

    /**
     * Get the atoms in a specific ring.
     *
     * @param ringNo ring index (0-based)
     * @return vector of atom indices (0-based)
     */
    const std::vector<int>& getRingAtoms(int ringNo) const;

    /**
     * Get the bonds in a specific ring.
     *
     * @param ringNo ring index (0-based)
     * @return vector of bond indices (0-based)
     */
    const std::vector<int>& getRingBonds(int ringNo) const;

private:
    // Reference to the molecule being analyzed
    OpenBabel::OBMol& mMol;

    // Ring data structures
    std::vector<std::vector<int>> mRingAtoms;        // Atoms in each ring (0-based indices)
    std::vector<std::vector<int>> mRingBonds;        // Bonds in each ring (0-based indices)
    std::vector<std::vector<int>> mAnnelatedRing;    // Annelated ring index for each bond position (-1 if none)
    std::vector<bool> mAromaticityHandled;           // Has ring aromaticity been determined?
    std::vector<bool> mIsAromaticRing;               // Is ring aromatic?
    std::vector<bool> mIsDelocalizedRing;            // Is ring delocalized?
    mutable std::vector<int> mHeteroPosition;        // Heteroatom position (-1 if none, mutable for lazy eval)
    std::unordered_set<int> mAzuleneBridgeBonds;     // Special azulene bridge bonds

    // Per-atom/bond feature flags (for O(1) lookup)
    std::vector<int> mAtomRingFeatures;              // Bit flags per atom
    std::vector<int> mBondRingFeatures;              // Bit flags per bond

    // Feature bit flags (matches Java's FEATURES_* constants)
    static constexpr int FEATURES_AROMATIC = 0x00010000;
    static constexpr int FEATURES_DELOCALIZED = 0x00020000;
    static constexpr int FEATURES_HETERO_AROMATIC = 0x00040000;

    // ========================================================================
    // Initialization methods
    // ========================================================================

    /**
     * Find all small rings in the molecule using OpenBabel's SSSR.
     */
    void findSmallRings();

    /**
     * Determine aromaticity for all rings.
     */
    void determineAromaticity();

    /**
     * Build annelation map - identify which bonds are shared between rings.
     */
    void buildAnnelationMap();

    /**
     * Determine aromaticity for a specific ring.
     * This is the core pattern-matching logic ported from Java.
     *
     * @param ringNo ring index
     * @return true if aromaticity was successfully determined
     */
    bool determineAromaticityForRing(int ringNo);

    /**
     * Propagate ring features (aromatic, delocalized, heteroaromatic)
     * to atoms and bonds for fast lookup.
     */
    void propagateRingFeatures();

    // ========================================================================
    // Helper methods (ported from Java's RingCollection)
    // ========================================================================

    /**
     * Check if an atom qualifies to be part of an aromatic system.
     * Ported from RingCollection.qualifiesAsAromaticAtom()
     *
     * @param atomIdx atom index (0-based)
     * @return true if atom can participate in aromatic system
     */
    bool qualifiesAsAromaticAtom(int atomIdx) const;

    /**
     * Check if a bond qualifies as a pi bond (double or triple).
     * Ported from RingCollection.qualifiesAsPiBond()
     *
     * @param bondIdx bond index (0-based)
     * @return true if bond is a pi bond
     */
    bool qualifiesAsPiBond(int bondIdx) const;

    /**
     * Build bond sequence for pattern matching.
     * Creates a binary pattern representing pi vs sigma bonds in the ring.
     * Handles annelated rings by treating shared aromatic ring bonds as pi bonds.
     *
     * @param ringNo ring index
     * @param aromaticButNotDelocalizedSequence output parameter for non-delocalized pattern
     * @return bond sequence as integer (binary pattern)
     */
    int buildBondSequence(int ringNo, int& aromaticButNotDelocalizedSequence) const;

    /**
     * Check if bond sequence matches aromatic pattern for given ring size.
     * Implements pattern matching for 3, 5, 6, and 7-membered rings.
     *
     * @param bondSequence binary pattern of pi bonds
     * @param ringSize number of atoms in ring
     * @param ringNo ring index (for atom access)
     * @param aromaticButNotDelocalizedSequence pattern for non-delocalized aromatics
     * @return true if pattern indicates aromaticity
     */
    bool checkAromaticPattern(int bondSequence, int ringSize, int ringNo,
                             int aromaticButNotDelocalizedSequence,
                             bool& hasDelocalizationLeak) const;

    // ========================================================================
    // Utility methods
    // ========================================================================

    /**
     * Convert OpenBabel 1-based index to 0-based.
     */
    inline int toBondIndex(const OpenBabel::OBBond* bond) const {
        return bond->GetIdx();
    }

    /**
     * Convert OpenBabel 1-based index to 0-based.
     */
    inline int toAtomIndex(const OpenBabel::OBAtom* atom) const {
        return atom->GetIdx() - 1;
    }
};

} // namespace OpenChemLib

#endif // AROMATICITY_RESOLVER_H

