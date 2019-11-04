/**********************************************************************
Spectrophore.h - Spectrophore(TM) calculator.
                 Declarations of OBSpectrophore

Copyright (C) 2005-2010 by Silicos NV

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

The algorithm in this software has been covered by patent WO2009146735.
However, Silicos NV and the inventors of the above mentioned patent assure
that no patent infringment claims will be issued against individuals or
institutions that use this software under the GNU General Public License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


#ifndef OB_SPECTROPHORE_H
#define OB_SPECTROPHORE_H

#include <openbabel/babelconfig.h>
#include <vector>

/** \file  spectrophore.h
 *  \brief Class to compute Spectrophores&tm;
 *
 *  \author Hans De Winter
 *  \date   2010-06-30
 *  \version 1
 **/




#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define N_PROPERTIES 4



namespace OpenBabel
{
  class OBMol;

/**
\class  OBSpectrophore spectrophore.h <openbabel/spectrophore.h>
\brief  Class to compute Spectrophores&tm;
\since  Version 2.3
<P>
<H3>Introduction</H3>
Spectrophores&tm; are one-dimensional descriptors generated from the property
fields surrounding the molecules. The technology allows the accurate description
of molecules in terms of their surface properties or fields. Comparison of
molecules’ property fields provides a robust structure-independent method of
aligning actives from different chemical classes. When applied to molecules such
as ligands and drugs, Spectrophores&tm; can be used as powerful molecular
descriptors in the fields of chemoinformatics, virtual screening, and QSAR
modeling.
</P><P>
The computation of Spectrophores&tm; is independent of the position and
orientation of the molecule and this enables easy and fast comparison of
Spectrophores&tm; between different molecules. Molecules having similar
three-dimensional properties and shapes always yield similar Spectrophores&tm;.
</P><P>
A Spectrophore&tm; is calculated by surrounding the three-dimensional
conformation of the molecule by a three-dimensional arrangement of points,
followed by calculating the interaction between each of the atom properties and
the surrounding the points. The three-dimensional arrangement of the points
surrounding the molecule can be regarded as an ‘artificial’ cage or receptor,
and the interaction calculated between the molecule and the cage can be regarded
as an artificial representation of an affinity value between molecule and cage.
Because the calculated interaction is dependent on the relative orientation of
the molecule within the cage, the molecule is rotated in discrete angles and the
most favorable interaction value is kept as final result. The angular stepsize
at which the molecule is rotated along its three axis can be specified by the
user and influences the accuracy of the method.
</P><P>
<H3>Atomic properties</H3>
The calculation of a Spectrophore&tm; starts by calculating the atomic
contributions of each property from which one wants to calculate a
Spectrophore&tm; from. In the current implementation, four atomic properties are
converted into a Spectrophore&tm;; these four properties include the atomic
partial charges, the atomic lipohilicities, the atomic shape deviations and the
atomic electrophilicities.
</P><P>
The atomic partial charges and atomic electrophilicity properties are calculated
using the electronegativity equalisation method (EEM) as described by Bultinck
and coworkers (<I>J. Phys. Chem.</I> 2002, <B>A106</B>, 7895-7901; <I>J. Chem.
Inf. Comput. Sci.</I> 2003, <B>43</B>, 422-428). The following table lists the
atomic electronegativity and hardness parameters that are used in the current
implementation of the EEM method:
<TABLE>
<TR><TD><B>Atom symbol</B></TD><TD><B>Atomic electronegativity</B></TD>
<TD><B>Atomic hardness</B></TD></TR>
<TR><TD>H</TD><TD><CENTER>0.20606</CENTER></TD><TD><CENTER>0.65971</CENTER></TD></TR>
<TR><TD>C</TD><TD><CENTER>0.36237</CENTER></TD><TD><CENTER>0.32966</CENTER></TD></TR>
<TR><TD>N</TD><TD><CENTER>0.49279</CENTER></TD><TD><CENTER>0.34519</CENTER></TD></TR>
<TR><TD>O</TD><TD><CENTER>0.73013</CENTER></TD><TD><CENTER>0.54428</CENTER></TD></TR>
<TR><TD>F</TD><TD><CENTER>0.72052</CENTER></TD><TD><CENTER>0.72664</CENTER></TD></TR>
<TR><TD>S</TD><TD><CENTER>0.62020</CENTER></TD><TD><CENTER>0.20640</CENTER></TD></TR>
<TR><TD>Cl</TD><TD><CENTER>0.36237</CENTER></TD><TD><CENTER>0.32966</CENTER></TD></TR>
<TR><TD>Br</TD><TD><CENTER>0.70052</CENTER></TD><TD><CENTER>0.54554</CENTER></TD></TR>
<TR><TD>I</TD><TD><CENTER>0.68052</CENTER></TD><TD><CENTER>0.30664</CENTER></TD></TR>
<TR><TD>Default</TD><TD><CENTER>0.36237</CENTER></TD><TD><CENTER>0.32966</CENTER></TD></TR>
</TABLE>
</P><P>
Atomic lipophilic potential parameters are calculated using a rule-based method
using parameters from the following table. These parameters were obtained by
fitting against the logP values of 10,881 molecules.
<TABLE>
<TR><TD><B>Atom</B></TD><TD><B>Lipophilicity parameter</B></TD></TR>
<TR><TD>H bound to C</TD><TD><CENTER>-0.018</CENTER></TD></TR>
<TR><TD>H bound to heteroatom</TD><TD><CENTER>-0.374</CENTER></TD></TR>
<TR><TD>C</TD><TD><CENTER>+0.271</CENTER></TD></TR>
<TR><TD>N</TD><TD><CENTER>-0.137</CENTER></TD></TR>
<TR><TD>O</TD><TD><CENTER>-0.321</CENTER></TD></TR>
<TR><TD>F</TD><TD><CENTER>+0.217</CENTER></TD></TR>
<TR><TD>S</TD><TD><CENTER>+0.385</CENTER></TD></TR>
<TR><TD>Cl</TD><TD><CENTER>+0.632</CENTER></TD></TR>
<TR><TD>Br</TD><TD><CENTER>+0.815</CENTER></TD></TR>
<TR><TD>I</TD><TD><CENTER>+0.198</CENTER></TD></TR>
<TR><TD>Default</TD><TD><CENTER>-0.175</CENTER></TD></TR>
</TABLE>
</P><P>
Finally, the atomic shape deviation is generated by calculating, for each atom,
the atom’s deviation from the average molecular radius. This is done in a
four steps process:
- The molecular center of geometry (COG) is calculated;
- The distances between each atom and the molecular COG are calculated;
- The average molecular radius is calculated by averaging all the atomic
distances.
- The distances between each atom and the COG are then divided by the average
molecular radius and centered on zero.
</P><P>
<H3>Interaction between the atoms and cage points</H3>
Following the calculation of all required atomic properties, the next step in
the calculation of a Spectrophore&tm; consists of determining the total
interaction value <I>V(c,p)</I> between each of the atomic contributions of
property <I>p</I> with a set of interaction points on an artificial cage
<I>c</I> surrounding the molecular conformation. For this purpose, each of these
interaction points <I>i</I> on cage <I>c</I> is assigned a value <I>P(c,i)</I>
which is either +1 or -1, with the constraint that the sum of all interaction
points on a particular cage should be zero. In a typical Spectrophore&tm;
calculation, a cage is represented as a rectangular box encompassing the
molecular conformation in all three dimensions, with the centers of the box
edges being the interaction points. Such a configuration gives twelve
interaction points per cage, and, in the case of a non-stereospecific
distribution of the interaction points, leads to 12 different cages. Although
there are no particular requirements as to the dimensions of the rectangular
cage, the distance between the interaction points and the geometrical extremes
of the molecule should be such that a meaningful interaction value between each
cage point and the molecular entity can be calculated. In this respect, the
default dimensions of the cage are constantly adjusted to enclose the molecule
at a minimum distance of 3 &Aring; along all dimensions. This cage size can be
modified by the user and influences the resolution of the Spectrophore&tm;.
</P><P>
\image html spectrophore_cage.png "Schematic representation of a molecule surrounded by the artifical cage"
The total interaction value <I>V(c,p)</I> between the atomic contribution values
<I>A(j,p)</I> of property <I>p</I> for a given molecular conformation and the
cage interaction values <I>P(c,i)</I> for a given cage <I>c</I> is calculated
according a standard interaction energy equation. It takes into account the
Euclidean distance between each atom and each cage point. This total interaction
<I>V(c,p)</I> for a given property <I>p</I> and cage <I>c</I> for a given
molecular conformation is minimized by sampling the molecular orientation along
the three axis in angular steps and the calculation of the interaction value for
each orientation within the cage. The final total interaction <I>V(c,p)</I> for
a given cage <I>c</I> and property <I>p</I> corresponds to the lowest
interaction value obtained this way, and corresponds to the <I>c</I>’th value in
the one-dimensional Spectrophore&tm; vector calculated for molecular property
<I>p</I>. As a result, a Spectrophore&tm; is organized as a vector of minimized
interaction values <I>V</I>, each of these organized in order of cages and
property values. Since for a typical Spectrophore&tm; implementation twelve
different cages are used, the total length of a Spectrophore&tm; vector equals
to 12 times the number of properties. Since four different properties are used
in the current implementation (electrostatic, lipophilic, electrophilic
potentials, and an additional shape index as described before), this leads to a
total Spectrophore&tm; length of 48 real values per molecular conformation.
</P><P>
Since Spectrophore&tm; descriptors are dependent on the actual
three-dimensional conformation of the molecule, a typical analysis includes the
calculation of Spectrophores&tm; from a reasonable set of different
conformations. It is then up to the user to decide on the most optimal strategy
for processing the different Spectrophore&tm; vectors. In a typical virtual
screening application, calculating the average Spectrophore&tm; vector from all
conformations of a single molecule may be a good strategy; other applications
have benefit from calculating a weighted average or the minimal values.
</P><P>
<H3>Accuracy</H3>
As already mentioned, the total interaction between cage and molecule for a
given property is minimized by sampling the molecular orientation in angular
steps of a certain magnitude. As a typical angular step size, 20º was found to
be the best compromise between accuracy and computer speed. Larger steps sizes
are faster to calculate but have the risk of missing the global interaction
energy minimum, while smaller angular steps sizes do sample the rotational space
more thoroughly but at a significant computational cost. The accuracy can be
specified by the user using the
OBSpectrophore::SetAccuracy(const OBSpectrophore::AccuracyOption) method.
</P><P>
<H3>Resolution</H3>
Spectrophores&tm; capture information about the property fields surrounding the
molecule, and the amount of detail that needs to be captured can be regulated by
the user. This is done by altering the minimal distance between the molecule and
the surrounding cage. The resolution can be specified by the user with the
OBSpectrophore::SetResolution(const double) method. The default distance along
all dimensions is 3.0 &Aring;. The larger the distance, the lower the resolution.
With a higher resolution, more details of the property fields surrounding the
molecule are contained by the Spectrophore&tm;. On the contrary, low resolution
settings may lead to a more general representation of the property fields, with
little or no emphasis on small local variations within the fields. Using a low
resolution can be the method of choice during the initial virtual screening
experiments in order to get an initial, but not so discriminative, first
selection. This initial selection can then further be refined during subsequent
virtual screening steps using a higher resolution. In this setting, small local
differences in the fields between pairs of molecules will be picked up much more
easily.
</P><P>
The absolute values of the individual Spectrophore&tm; data points are dependent
on the used resolution. Low resolution values lead to small values of the
calculated individual Spectrophore&tm; data points, while high resolutions will
lead to larger data values. It is therefore only meaningful to compare only
Spectrophores&tm; that have been generated using the same resolution settings or
after some kind of normalization is performed.
</P><P>
Computation time is not influenced by the specified resolution, hence the
computation time is identical for all different resolution settings.
</P><P>
<H3>Stereospecificity</H3>
Some of the cages that are used to calculated Spectrophores&tm; have a
stereospecific distribution of the interaction points. The resulting
interaction valus resulting from these cages are therefore sensitive to the
enantiomeric configuration of the molecule within the cage. The fact that both
stereoselective as well as stereo non-selective cages can be used makes it
possible to include or exclude stereospecificity in the virtual screening
search. Depending on the desired output, the stereospecificity of
Spectrophores&tm; can be specified by the user:
- No stereospecificity (default). Spectrophores&tm; are generated using cages that are not
stereospecific. For most applications, these Spectrophores&tm; will suffice.
- Unique stereospecificity. Spectrophores&tm; are generated using unique
stereospecific cages.
- Mirror stereospecificity. Mirror stereospecific Spectrophores&tm; are
Spectrophores&tm; resulting from the mirror enantiomeric form of the input
molecules.
</P><P>
The stereospecificity can be specified by the user using the
OBSpectrophore::SetStereo(const OBSpectrophore::StereoOption) method.
</P><P>
The differences between the corresponding data points of unique and mirror
stereospecific Spectrophores&tm; are very small and require very long
calculation times to obtain a sufficiently high quality level. This increased
quality level is triggered by the accuracy setting and will result in
calculation times being increased by at least a factor 100. As a consequence, it
is recommended to apply this increased accuracy only in combination with a
limited number of molecules, and when the small differences between the
stereospecific Spectrophores&tm; are really critical. However, for the vast
majority of virtual screening applications, this increased accuracy is not
required as long as it is not the intention to draw conclusions about
differences in the underlying molecular stereoselectivity. Non-stereospecific
Spectrophores&tm; will therefore suffice for most applications.
</P><P>
<H3>Normalisation</H3>
It may sometimes be desired to focus on the relative differences between the
Spectrophore&tm; data points rather than focussing on the absolute differences.
In these cases, normalization of Spectrophores&tm; may be required. The current
implementation offers with the
OBSpectrophore::SetNormalization(const OBSpectrophore::NormalizationOption)
method the possibility to normalize in four different ways:
- No normalization (default);
- Normalization towards zero mean;
- Normalization towards standard deviation;
- Normalization towards zero mean and unit standard deviation.
</P><P>
In all these cases, normalization is performed on a ‘per-property’ basis, which
means that the data points belonging to the same property set are treated as a
single set and that normalization is only performed on the data points within
each of these sets and not across all data points.
</P><P>
Normalization may be important when comparing the Spectrophores&tm; of charged
molecules with those of neutral molecules. For molecules carrying a global
positive charge, the resulting Spectrophore&tm; data points of the charge and
electrophilicity properties will both be shifted in absolute value compared to
the corresponding data points of the respective neutral species. Normalization
of the Spectrophores&tm; removes the original magnitude differences for the data
points corresponding to the charge and electrophilicity properties of charged
and neutral species. Therefore, if the emphasis of the virtual screening
consists of the identification of molecules with similar property fields without
taking into account differences in absolute charge, then Spectrophores&tm;
should be normalized towards zero mean. However, if absolute charge differences
should be taken into account to differentiate between molecules, unnormalized
Spectrophores&tm; are recommended.
</P>

@since version 2.3
*/
   class OBAPI OBSpectrophore
   {
      public:

         /** Accuracy options
         */
         enum AccuracyOption {AngStepSize1,
                              AngStepSize2,
                              AngStepSize5,
                              AngStepSize10,
                              AngStepSize15,
                              AngStepSize20,
                              AngStepSize30,
                              AngStepSize36,
                              AngStepSize45,
                              AngStepSize60};

         /** Normalisation options
         */
         enum NormalizationOption {NoNormalization,
                                   NormalizationTowardsZeroMean,
                                   NormalizationTowardsUnitStd,
                                   NormalizationTowardsZeroMeanAndUnitStd};

         /** Stereo options
         */
         enum StereoOption {NoStereoSpecificProbes,
                            UniqueStereoSpecificProbes,
                            MirrorStereoSpecificProbes,
                            AllStereoSpecificProbes};

         //@}

         /** \name Structors */
         //@{
         /** Default constructor to create an OBSpectrophore object.
         */
         OBSpectrophore(void);

         /** Copy constructor to create an OBSpectrophore object by taking
         a copy from another OBSpectrophore object.

         \param   sphore A reference to the source OBSpectrophore object
         */
         OBSpectrophore(const OBSpectrophore& sphore);

         /** Default destructor object
         */
         virtual ~OBSpectrophore(void);
         //@}

         /** \name Overloaded operators */
         //@{
         /** Assignment operator that assigns a source OBSpectrophore object to
         the target OBSpectrophore object.

         \param   sphore A reference to the source OBSpectrophore object
         \return  A reference to the new target OBSpectrophore object
         */
         OBSpectrophore& operator=(const OBSpectrophore& sphore);
         //@}

         /** \name Set methods */
         //@{
         /** Method to set the resolution at which Spectrophores&tm; will be
         calculated.

         Sets the resolution at which Spectrophores&tm; are calculated. The
         resolution is an arbitrary positive number larger than 0.0, and is used
         to increase the box size for the calculation of Spectrophores&tm; in
         all directions. For example, a resolution of 3 means that the box size
         will be increased by a value of 3 &Aring; in all directions. Smaller
         values for the resolution have the effect that small differences in the
         atomic properties will result in more enhanced differences in the
         resulting Spectrophores&tm;, while larger resolution values will result
         in the smoothing of local property differences. All values &lt;= 0.0
         are automatically reset to the default value of 3.0.

         \param   r The desired resolution expressed as a double number
         */
         void SetResolution(const double r = 3.0);

         /** Method to set the accuracy at which Spectrophores&tm; will be
         calculated.

         Sets the accuracy by which Spectrophores&tm; are calculated. The
         accuracy is linked to the angular step increment that is used to
         calculate Spectrophores&tm;. The accuracy parameter is an enumeration
         type defined by OBSpectrophore::AccuracyOption. should be a number ranging
         between 0 and 9 (inclusive), where a value of 9 means a very high
         accuracy, and 0 means no accuracy at all. Values between 2-6 are a good
         compromise between accuracy and speed. The following table provides the
         link between the parameter value and the angular step size:
         <TABLE>
         <TR><TD><B>Accuracy</B></TD><TD><B>OBSpectrophore::AccuracyOption parameter</B></TD><TD><B>Minimal angular steps size</B></TD></TR>
         <TR><TD>Extremely high</TD><TD><CENTER>AngStepSize1</CENTER></TD><TD><CENTER>1º</CENTER></TD></TR>
         <TR><TD>Extremely high</TD><TD><CENTER>AngStepSize2</CENTER></TD><TD><CENTER>2º</CENTER></TD></TR>
         <TR><TD>Very High</TD><TD><CENTER>AngStepSize5</CENTER></TD><TD><CENTER>5º</CENTER></TD></TR>
         <TR><TD>Very high</TD><TD><CENTER>AngStepSize10</CENTER></TD><TD><CENTER>10º</CENTER></TD></TR>
         <TR><TD>High</TD><TD><CENTER>AngStepSize15</CENTER></TD><TD><CENTER>15º</CENTER></TD></TR>
         <TR><TD>Default</TD><TD><CENTER>AngStepSize20</CENTER></TD><TD><CENTER>20º</CENTER></TD></TR>
         <TR><TD>Low</TD><TD><CENTER>AngStepSize30</CENTER></TD><TD><CENTER>30º</CENTER></TD></TR>
         <TR><TD>Very low</TD><TD><CENTER>AngStepSize36</CENTER></TD><TD><CENTER>36º</CENTER></TD></TR>
         <TR><TD>Extremely low</TD><TD><CENTER>AngStepSize45</CENTER></TD><TD><CENTER>45º</CENTER></TD></TR>
         <TR><TD>Extremely low</TD><TD><CENTER>AngStepSize60</CENTER></TD><TD><CENTER>60º</CENTER></TD></TR>
         </TABLE>

         \param  a The desired accuracy expressed as an
         OBSpectrophore::AccuracyOption enumeration type
         */
         void SetAccuracy(const AccuracyOption a = AngStepSize20);

         /** Method to set the required stereoselectivity of the resulting
         Spectrophores&tm;.

         Sets the stereoselectivity of the generated Spectrophores&tm;. This is
         achieved by selecting the appropriate probes for the calculation. In
         the default case that no stereoselectivity is required, only
         non-stereospecific probes are used in the calculation. Stereospecific
         differences between enantiomers are only captured in combination with a
         very high accuracy, hence leading to very long computation times. The
         following table provides the link between the parameter value and the
         desired stereochemical treatment:
         <TABLE>
         <TR><TD><B>Stereo treatment</B></TD><TD><B>OBSpectrophore::StereoOption parameter</B></TD></TR>
         <TR><TD>Non-stereospecific probes (default)</TD><TD><CENTER>NoStereoSpecificProbes</CENTER></TD></TR>
         <TR><TD>Only the mirror stereospecific probes</TD><TD><CENTER>UniqueStereoSpecificProbes</CENTER></TD></TR>
         <TR><TD>Only the unique stereospecific probes</TD><TD><CENTER>MirrorStereoSpecificProbes</CENTER></TD></TR>
         <TR><TD>All stereospecific probes (unique and mirror)</TD><TD><CENTER>AllStereoSpecificProbes</CENTER></TD></TR>
         </TABLE>

         \param   s The desired stereospecificity treatment expressed as an
         OBSpectrophore::StereoOption enumeration type
         */
         void SetStereo(const StereoOption s = NoStereoSpecificProbes);

         /** Method to set the desired normalization treatment of the calculated
         Spectrophores&tm;.

         The following table provides the link between the parameter value
         and the desired normalization treatment:
         <TABLE>
         <TR><TD><B>Normalization treatment</B></TD><TD><B>OBSpectrophore::NormalizationOption parameter</B></TD></TR>
         <TR><TD>No normalisation (default)</TD><TD><CENTER>NoNormalization</CENTER></TD></TR>
         <TR><TD>Normalization towards zero mean</TD><TD><CENTER>NormalizationTowardsZeroMean</CENTER></TD></TR>
         <TR><TD>Normalization towards unit standard deviation</TD><TD><CENTER>NormalizationTowardsUnitStd</CENTER></TD></TR>
         <TR><TD>Normalization towards zero mean and unit standard deviation</TD><TD><CENTER>NormalizationTowardsZeroMeanAndUnitStd</CENTER></TD></TR>
         </TABLE>
         </P>

         \param   n The desired normalization treatment expressed as an
         OBSpectrophore::NormalizationOption enumeration type
         */
         void SetNormalization(const NormalizationOption n = NoNormalization);
         //@}

         /** \name Get methods */
         //@{
         /** Returns the accuracy at which Spectrophores&tm; will be calculated.

         \return  The accuracy expressed as an OBSpectrophore::AccuracyOption
                  enumeration type. For a link between the returned value and
                  the angular step size, see the
                  OBSpectrophore::SetAccuracy(const OBSpectrophore::AccuracyOption)
                  method.
         */
         AccuracyOption GetAccuracy(void) const;

         /** Returns the resolution at which Spectrophores&tm; will be calculated.

         \return  The resolution in &Aring; units
         */
         double GetResolution(void) const;

         /** Returns the stereoselectivity setting at which Spectrophores&tm;
         will be calculated.

         \return  The level of stereospecificity expressed as an
                  OBSpectrophore::StereoOption enumeration type. For a link
                  between the returned value and the stereoselectivity
                  treatment, see the
                  OBSpectrophore::SetStereo(const OBSpectrophore::StereoOption)
                  method.
         */
         StereoOption GetStereo(void) const;

         /** Returns the normalization settings at which Spectrophores&tm; will
         be calculated.

         \return  The normalization treatmet expressed as an
                  OBSpectrophore::NormalizationOption enumeration type.
                  For a link between the returned value and the normalization
                  treatment, see the
                  OBSpectrophore::SetNormalization(const OBSpectrophore::NormalizationOption)
                  method.
         */
         NormalizationOption GetNormalization(void) const;

         /** Calling this method starts the calculation of the Spectrophore&tm;.
         After succesful calculation, the Spectrophore&tm; is returned as a
         standard vector of 48 doubles. The 48 doubles are organised into 4 sets
         of 12 doubles each:-
         - numbers 01-11: Spectrophore&tm; values calculated from the atomic partial charges;
         - numbers 13-24: Spectrophore&tm; values calculated from the atomic lipophilicity properties;
         - numbers 25-36: Spectrophore&tm; values calculated from the atomic shape deviations;
         - numbers 37-48: Spectrophore&tm; values calculated from the atomic electrophilicity properties;

\code
OpenBabel::OBConversion obconversion;
obconversion.SetInFormat("sdf");
OpenBabel::OBMol mol;
OpenBabel::OBSpectrophore s;
s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);          // Default
s.SetResolution(3.0);                                             // Default
s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);   // Default
s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);   // Default
std::ifstream ifs;
ifs.open(argv[1]);
while (obconversion.Read(&mol, &ifs))
{
   std::vector<double> result = s.GetSpectrophore(&mol);
   for (unsigned int i(0); i < result.size(); ++i)
   {
      if (!(i%12)) std::cerr << std::endl;
      std::cerr << result[i] << "  ";
   }
   std::cerr << std::endl;
   mol.Clear();
}
\endcode

         \param   mol Pointer to the OBMol object from which to calculate a
                  Spectrophore&tm;. For proper functioning, the input molecule
                  should have all explicit hydrogens assigned and the molecule
                  should have a three-dimensional conformation assigned. It is
                  the responsability of the programmer to make sure that the
                  molecule is in the desired protonation state.

         \return  The calculated Spectrophore&tm; as a standard vector of 48
                  doubles. An empty vector of doubles is returned in case an
                  error occurred during the calculation. For example, this could
                  happen in the case that the molecule contains less than three
                  atoms. It is therefore up to the user to check the size of the
                  returned vector to capture errors.
         */
         std::vector<double> GetSpectrophore(OpenBabel::OBMol* mol);
         //@}


      protected:

         //@{
         double               _resolution;
         AccuracyOption       _accuracy;
         StereoOption         _stereoFlag;
         NormalizationOption  _normalization;
         std::vector<int>     _rotationStepList;
         unsigned int         _nAtoms;
         double**             _property;
         double*              _radii;
         double**             _oricoor;
         double**             _coor;
         unsigned int         _beginProbe;
         unsigned int         _endProbe;
         unsigned int         _numberOfProbes;
         std::vector<double>  _spectro;
         struct
         {
            int value[12];
         }                    _probe[48];
         struct
         {
            double x;
            double y;
            double z;
            double v[N_PROPERTIES];
         }                    _boxPoint[12];
         //@}

         //@{
         void _getMoleculeData(OpenBabel::OBMol*);
         void _orient(void);
         void _getBox(double**);
         void _setBox(void);
         void _getEnergies(double**, double*);
         void _initiateSpectrophore(double*, double*);
         void _rotateX(double**, double**, const double, const double);
         void _rotateY(double**, double**, const double, const double);
         void _rotateZ(double**, double**, const double, const double);
         void _updateSpectrophore(double*, double*);
         void _calculateProperties(OpenBabel::OBMol*);
         void _solveMatrix(double**, double*, unsigned int);
         void _luDecompose(double**, std::vector<int>&, unsigned int);
         void _luSolve(double**, std::vector<int>&, double*, unsigned int);
         void _swapRows(double*, unsigned int, unsigned int);
         void _swapRows(double**, unsigned int, unsigned int, unsigned int);
         //@}

   };

} // end namespace OpenBabel


#endif   //OB_SPECTROPHORE_H
