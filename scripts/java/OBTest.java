public class OBTest
{
  public OBTest()
  {
    System.loadLibrary("openbabeljni");
  }

  public void run()
  {
    OBConversion conv = new OBConversion();
    conv.SetInFormat("smi");

    OBMol mol = new OBMol();
    conv.ReadString(mol, "CCCCC");
    System.out.println(mol.NumAtoms());
  }

  public static void main(String[] args)
  {
    System.out.println("Running OBTest...");

    OBTest test = new OBTest();

    test.run();
  }
}
