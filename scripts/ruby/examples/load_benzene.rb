require 'openbabel'

conv = OpenBabel::OBConversion.new
conv.set_in_and_out_formats("smi", "mdl")

mol = OpenBabel::OBMol.new
conv.read_string(mol, 'c1ccccc1')

puts "Benzene has #{mol.num_atoms} heavy atoms." 

mol.add_hydrogens

puts "Benzene has #{mol.num_atoms} total atoms."
puts "The molecular weight of benzene is #{ mol.get_mol_wt}."
