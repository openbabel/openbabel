import org.openbabel._

object OBTest {

  def main(args: Array[String]) {
    println("Loading openbabel_java...")
    System.loadLibrary("openbabel_java")

    println("Running OBTest...")
    run
  }

  def run() {
    val c = new OBConversion
    val mol = new OBMol

    c.SetInFormat("smi")
    c.ReadString(mol, "c1ccccc1")

    println("Benzene has " + mol.NumAtoms + " heavy atoms.")

    mol.AddHydrogens

    println("Benzene has " + mol.NumAtoms + " atoms in total.")
    println("The molecular weight of benzene is " + mol.GetMolWt)
  }
}
