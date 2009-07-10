"%JAVA_HOME%\bin\javac" org\openbabel\*.java
del /Q org\openbabel\*.java
"%JAVA_HOME%\bin\jar" cf openbabel.jar org
del /Q org\openbabel\*.class