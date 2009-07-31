require 'openbabel'

class Convertor
  def initialize
    @conversion = OpenBabel::OBConversion.new
    @conversion.set_in_and_out_formats 'smi', 'can'
  end

  def convert smiles
    mol = OpenBabel::OBMol.new

    @conversion.read_string mol, smiles
    @conversion.write_string mol
  end
end
 
c=Convertor.new
aminopterin_smiles = "C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=NC(=N3)N)N" # from PubChem CID 2154

puts "The canonical SMILES for aminopterin is:\n#{c.convert(aminopterin_smiles)}"
