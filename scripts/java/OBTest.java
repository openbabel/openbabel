import org.openbabel.*

public class OBTest
{
  public OBTest()
  {
    System.loadLibrary("openbabel");
  }

  public void run()
  {
    OBConversion c = new OBConversion();
  }

  public static void main(String[] args)
  {
    System.out.println("Running OBTest...");

    OBTest test = new OBTest();

    test.run();
  }
}
